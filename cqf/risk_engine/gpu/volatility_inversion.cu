#include <math.h>
#include <stdio.h>

#include "cuda_profiler_api.h"

#include "volatility_inversion.cuh"


#ifndef PI 
#define PI 3.141592653589793238462643f
#endif 

const int serial_threads_per_block = 448;

const int parallel_guesses_per_option = 4;
const int parallel_threads_per_block = 192;

const int vega_parallel_guesses_per_option = 4;
const int vega_parallel_threads_per_block = 192;


#define HANDLE_CUDA_ERROR(err)	if (err) { printf("%s", cudaGetErrorString(err)); return; }

__device__ __host__ float normal_pdf(const float z) 
{
	return (1.0f / sqrt(2.0f * PI)) * exp(-0.5f * z);
}


__device__ __host__ float normal_cdf(const float z)
{
    const float b1 =  0.31938153f; 
    const float b2 = -0.356563782f; 
    const float b3 =  1.781477937f;
    const float b4 = -1.821255978f;
    const float b5 =  1.330274429f; 
    const float p  =  0.2316419f; 
    const float c2 =  0.3989423f; 

    if (z >  6.0f)
	{
		return 1.0f;
	}

    if (z < -6.0f)
	{
		return 0.0f;
	}

    const float a = abs(z); 
    const float t = 1.0f / (1.0f + a * p); 
    const float b = c2 * exp((-z) * (z / 2.0f)); 
    
	float n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t; 
    n = 1.0f - b * n; 
    if ( z < 0.0f )
	{
		n = 1.0f - n;
	}
    return n; 
}


__device__ __host__ float call_price(const float s, const float r, const float v, const float t, const float k)
{
	const float sqrt_t = sqrt(t);
	const float d1 = (1.0f / (v * sqrt_t)) * (log(s / k) + (r + v * v * 0.5f) * t);
	const float d2 = d1 - v * sqrt_t;

	return (s * normal_cdf(d1)) - (k * exp(-r * t) * normal_cdf(d2));
}


__device__ __host__ float call_vega(const float s, const float r, const float v, const float t, const float k)
{
	const float sqrt_t = sqrt(t);
	const float d1 = (1.0f / (v * sqrt_t)) * (log(s / k) + (r + v * v * 0.5f) * t);
	return s * normal_pdf(d1) * sqrt_t;
}


__device__ __host__ float put_price(const float s, const float r, const float v, const float t, const float k)
{
	const float sqrt_t = sqrt(t);
	const float d1 = (1.0f / (v * sqrt_t)) * (log(s / k) + (r + v * v * 0.5f) * t);
	const float d2 = d1 - v * sqrt_t;

	return (k * exp(-r * t) * normal_cdf(-d2)) - (s * normal_cdf(-d1));
}


__device__ __host__ float put_vega(const float s, const float r, const float v, const float t, const float k)
{
	return call_vega(s, r, v, t, k);
}


__global__ void volatility_inversion_device(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num)
{
	const int option = (blockIdx.x * blockDim.x) + threadIdx.x;
	if (option >= num)
	{
		return;
	}

	float v_local = v[option];
	float error;
	int i = 0;
	do
	{
		const float price = call_price(s, r, v_local, t[option], k[option]);
		const float vega = call_vega(s, r, v_local, t[option], k[option]);
		
		error = (p[option] - price);
		v_local += error / vega;
	} while ((abs(error) > tol) && (i++ < iter));

	v[option] = v_local;
}

template<int GUESSES, int BLOCK_SIZE>
__device__ int find_best_guess(const float *const ferr_ladder, const int tid, const int lower_guess)
{
	__shared__ int smallest_idx[BLOCK_SIZE];
	__shared__ float smallest_error[BLOCK_SIZE];

	/* Find minimum error */
	smallest_idx[tid] = tid;
	smallest_error[tid] = ferr_ladder[tid];
	if ((GUESSES > 16) && (smallest_error[tid + 16] < smallest_error[tid]))
	{
		smallest_error[tid] = smallest_error[tid + 16];
		smallest_idx[tid] = tid + 16;
	}

	if ((GUESSES > 8) && (smallest_error[tid + 8] < smallest_error[tid]))
	{
		smallest_error[tid] = smallest_error[tid + 8];
		smallest_idx[tid] = tid + 8;
	}

	if ((GUESSES > 4) && (smallest_error[tid + 4] < smallest_error[tid]))
	{
		smallest_error[tid] = smallest_error[tid + 4];
		smallest_idx[tid] = tid + 4;
	}

	if ((GUESSES > 2) && (smallest_error[tid + 2] < smallest_error[tid]))
	{
		smallest_error[tid] = smallest_error[tid + 2];
		smallest_idx[tid] = tid + 2;
	}

	if (smallest_error[tid + 1] < smallest_error[tid])
	{
		smallest_error[tid] = smallest_error[tid + 1];
		smallest_idx[tid] = tid + 1;
	}
	
	return smallest_idx[lower_guess];
}


template<int GUESSES>
__global__ void parallel_volatility_inversion_device(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter)
{
	const int tid = threadIdx.x;
	const int guess = tid & (GUESSES - 1);
	const int lower_guess = tid & ~(GUESSES - 1);
	const int upper_guess = lower_guess + GUESSES - 1;
	const int log_guesses = 31 - __clz(GUESSES);
	const int option = (blockIdx.x * (parallel_threads_per_block >> log_guesses)) + (tid >> log_guesses);

	/* Build guesses */
	__shared__ float v_ladder[parallel_threads_per_block];
	__shared__ float err_ladder[parallel_threads_per_block];
	__shared__ float ferr_ladder[parallel_threads_per_block];

	float ladder_span = 0.04f;	/* Span the guesses over 4% */
	float v_min = v[option] - (ladder_span * 0.5f);

	int i = 0;
	const float guess_fraction = guess / static_cast<float>(GUESSES - 1);
	do
	{
		/* Work out guess */
		v_ladder[tid] = v_min + (ladder_span * guess_fraction);

		/* Price */
		err_ladder[tid] = call_price(s, r, v_ladder[tid], t[option], k[option]) - p[option];
		ferr_ladder[tid] = abs(err_ladder[tid]);
		
		/* Find minimum error */
		const int min_err_pos = find_best_guess<GUESSES, parallel_threads_per_block>(ferr_ladder, tid, lower_guess);
		if (ferr_ladder[min_err_pos] < tol)
		{
			if (tid == lower_guess)
			{
				v[option] = v_ladder[min_err_pos];
			}
			break;
		}

		/* Pick the span for the next ladder */
		/* Doesnt matter if v_min is actually higher than v_max so long as 0 is crossed */
		if ((err_ladder[lower_guess] * err_ladder[upper_guess]) >= 0.0f) /* Root not found (or two roots found) */
		{
			ladder_span *= 2.0f;
			if (abs(err_ladder[lower_guess] - err_ladder[upper_guess]) < tol) /* Ladder is very flat so no direction */
			{
				v_min -= ladder_span * 0.5f;
				ladder_span *= 2.0f;
			}
			else if (ferr_ladder[lower_guess] < ferr_ladder[upper_guess]) /* Lower end is closer to root so expand it */
			{
				v_min -= ladder_span;
			}
			else /* Upper end is closer to root so expand it */
			{
				v_min = v_ladder[upper_guess];
			}
		}
		else if (min_err_pos == lower_guess) /* Root found at lower extreme of ladder */
		{
			if ((err_ladder[lower_guess] * err_ladder[lower_guess + 1]) >= 0.0f)
			{
				ladder_span *= 2.0f;
				v_min = (v_ladder[min_err_pos] - ladder_span);
			}
			else
			{
				ladder_span *= 0.5f;
				v_min = v_ladder[min_err_pos];
			}
		}
		else if (min_err_pos == upper_guess) /* Root found at upper extreme of ladder */
		{
			if ((err_ladder[upper_guess] * err_ladder[upper_guess - 1]) >= 0.0f)
			{
				ladder_span *= 2.0f;
				v_min = v_ladder[min_err_pos];
			}
			else
			{
				ladder_span *= 0.5f;
				v_min = v_ladder[upper_guess - 1];
			}
		}
		else if ((err_ladder[min_err_pos] * err_ladder[min_err_pos - 1]) < 0.0f) /* Root in bin below min error */
		{
			ladder_span *= 1.0f / GUESSES;
			v_min = v_ladder[min_err_pos - 1];
		}
		else /* Root in bin above min error */
		{
			ladder_span *= 1.0f / GUESSES;
			v_min = v_ladder[min_err_pos];
		}
	} while (i++ < iter);
}


template<int GUESSES>
__global__ void vega_guided_parallel_volatility_inversion_device(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter)
{
	const int tid = threadIdx.x;
	const int guess = tid & (GUESSES - 1);
	const int lower_guess = tid & ~(GUESSES - 1);
	const int upper_guess = lower_guess + GUESSES - 1;
	const int log_guesses = 31 - __clz(GUESSES);
	const int option = (blockIdx.x * (parallel_threads_per_block >> log_guesses)) + (tid >> log_guesses);

	/* Build guesses */
	__shared__ float v_ladder[vega_parallel_threads_per_block];
	__shared__ float vega_ladder[vega_parallel_threads_per_block];
	__shared__ float err_ladder[vega_parallel_threads_per_block];
	__shared__ float ferr_ladder[vega_parallel_threads_per_block];

	float ladder_span = 0.04f;	/* Span the guesses over 4% */
	float v_mid = v[option];

	int i = 0;
	int min_err_pos;
	const float guess_fraction = guess / static_cast<float>(GUESSES);
	do
	{
		/* Work out guess */
		v_ladder[tid] = v_mid + (ladder_span * guess_fraction) - (ladder_span * 0.5f) + (ladder_span / static_cast<float>(GUESSES)) * 0.5f;

		/* Price */
		err_ladder[tid] = p[option] - call_price(s, r, v_ladder[tid], t[option], k[option]);
		vega_ladder[tid] = call_vega(s, r, v_ladder[tid], t[option], k[option]);
		ferr_ladder[tid] = abs(err_ladder[tid]);
		
		/* Find minimum error */
		min_err_pos = find_best_guess<GUESSES, vega_parallel_threads_per_block>(ferr_ladder, tid, lower_guess);

		/* Pick the span for the next ladder */
		v_mid = v_ladder[min_err_pos] + (err_ladder[min_err_pos] / vega_ladder[min_err_pos]);
		if ((v_mid < v_ladder[upper_guess]) && (v_mid > v_ladder[lower_guess]))
		{
			ladder_span *= (1.0f / 1024.0f);
		}
		else
		{
			ladder_span *= 2.0f;
		}
	} while ((ferr_ladder[min_err_pos] > tol) && (i++ < iter));

	if (tid == lower_guess)
	{
		v[option] = v_ladder[min_err_pos];
	}
}


void volatility_inversion(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num, const kernel_t kernel)
{
	/* Start the profiler */
	//cudaProfilerStart();

	/* Get device */
	HANDLE_CUDA_ERROR(cudaSetDevice(0));

	/* Allocate and copy memory */
	float *dev_v;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_v, num * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_v, v, num * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_t;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_t, num * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_t, t, num * sizeof(float), cudaMemcpyHostToDevice));
	
	float *dev_k;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_k, num * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_k, k, num * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_p;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_p, num * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_p, p, num * sizeof(float), cudaMemcpyHostToDevice));

	/* Run kernels */
	const int serial_threads = min(num, serial_threads_per_block);
	const int serial_blocks = (int)ceil(num / static_cast<float>(serial_threads));

	const int parallel_threads = parallel_threads_per_block;
	const int parallel_blocks = static_cast<int>(ceil((parallel_guesses_per_option * num) / static_cast<float>(parallel_threads_per_block)));

	const int vega_parallel_threads = vega_parallel_threads_per_block;
	const int vega_parallel_blocks = static_cast<int>(ceil((vega_parallel_guesses_per_option * num) /  static_cast<float>(vega_parallel_threads_per_block)));
	switch (kernel)
	{
		case serial :
			volatility_inversion_device<<<serial_blocks, serial_threads>>>(s, r, dev_v, dev_t, dev_k, dev_p, tol, iter, num);
			break;
		case parallel :
			parallel_volatility_inversion_device<parallel_guesses_per_option><<<parallel_blocks, parallel_threads>>>(s, r, dev_v, dev_t, dev_k, dev_p, tol, iter);
			break;
		case vega_parallel :
			vega_guided_parallel_volatility_inversion_device<vega_parallel_guesses_per_option><<<vega_parallel_blocks, vega_parallel_threads>>>(s, r, dev_v, dev_t, dev_k, dev_p, tol, iter);
			break;
	}

	/* Copy data back */
	HANDLE_CUDA_ERROR(cudaDeviceSynchronize());
	HANDLE_CUDA_ERROR(cudaMemcpy(v, dev_v, num * sizeof(float), cudaMemcpyDeviceToHost));

	/* Clean up */
	HANDLE_CUDA_ERROR(cudaFree(dev_v));
	HANDLE_CUDA_ERROR(cudaFree(dev_t));
	HANDLE_CUDA_ERROR(cudaFree(dev_k));
	HANDLE_CUDA_ERROR(cudaFree(dev_p));

	/* Flush profiling info */
	HANDLE_CUDA_ERROR(cudaDeviceReset());

	/* Stop the profiler */
	//cudaProfilerStop();
}