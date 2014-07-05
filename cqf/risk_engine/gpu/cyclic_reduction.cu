#include <stdio.h>

const int max_threads_per_block = 512;

#define HANDLE_CUDA_ERROR(err)	if (err) { printf("%s", cudaGetErrorString(err)); return; }



__device__ int cyclic_reduction_forward_reduction(float *lower, float *diagonal, float *upper, float *equal, const int dim, int step, const int to)
{
	/* Forward reduction */
	for (; (step * to * 3) <= dim; step <<= 1)
	{
		const int addr = (threadIdx.x * (step << 1)) + (step << 1) - 1;
		if (addr < dim)
		{
			if (addr - step >= 0)
			{
				const float alpha = -lower[addr] / diagonal[addr - step];
				equal[addr]    += (alpha * equal[addr - step]);
				diagonal[addr] += (alpha * upper[addr - step]);
				lower[addr]		= alpha * lower[addr - step];
			}

			if (addr + step < dim)
			{
				const float gamma = -upper[addr] / diagonal[addr + step];
				equal[addr]	   += (gamma * equal[addr + step]);
				diagonal[addr] += (gamma * lower[addr + step]);
				upper[addr]		= gamma * upper[addr + step];
			}
		}
		__syncthreads();
	}

	return step;
}


__device__ void cyclic_reduction_back_substitution(float *lower, float *diagonal, float *upper, float *equal, const int dim, int step, const int to)
{
	/* Backward substitution */
	for (; step > to; step >>= 1)
	{
		const int addr = (threadIdx.x * (step << 1)) + step - 1;
		if (addr < dim)
		{
			if (addr - step >= 0)
			{
				equal[addr] -= (lower[addr] * equal[addr - step]);
			}

			if (addr + step < dim)
			{
				equal[addr] -= (upper[addr] * equal[addr + step]);
			}

			equal[addr] = equal[addr] / diagonal[addr];
		}
		__syncthreads();
	}
}


__global__ void cyclic_reduction_device(float *lower_glb, float *diagonal_glb, float *upper_glb, float *equal_glb, const int dim)
{
	__shared__ float lower[512];
	__shared__ float diagonal[512];
	__shared__ float upper[512];
	__shared__ float equal[512];

	lower[threadIdx.x] = lower_glb[threadIdx.x];
	diagonal[threadIdx.x] = diagonal_glb[threadIdx.x];
	upper[threadIdx.x] = upper_glb[threadIdx.x];
	equal[threadIdx.x] = equal_glb[threadIdx.x];
	__syncthreads();

	/* Forward reduction */
	int step = cyclic_reduction_forward_reduction(lower, diagonal, upper, equal, dim, 1, 1);

	/* Solve base system */
	if (threadIdx.x == 0)
	{
		if ((dim / step) == 2) /* Solve simultaneous equations */
		{
			const int equal_addr = (step << 1) - 1;
			const float a0 = diagonal[equal_addr - step];
			const float a1 = lower[equal_addr];
			const float b0 = upper[equal_addr - step];
			const float b1 = diagonal[equal_addr];
			const float c0 = equal[equal_addr - step];
			const float c1 = equal[equal_addr];

			equal[equal_addr] = (c0 * a1 - a0 * c1) / (a1 * b0 - a0 * b1);
			equal[equal_addr - step] = (c0 - b0 * equal[equal_addr]) / a0;
		}
		else /* blk_size == 1, equations are already solved */
		{
			const int equal_addr = step - 1;
			equal[equal_addr] = equal[equal_addr] / diagonal[equal_addr];
		}
	}

	__syncthreads();
	step >>= 1;

	/* Backward substitution */
	cyclic_reduction_back_substitution(lower, diagonal, upper, equal, dim, step, 0);
	equal_glb[threadIdx.x] = equal[threadIdx.x];
}


void cyclic_reduction(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	const int log_dim = static_cast<int>(ceil(log(static_cast<float>(dim)) / log(2.0f)));

	/* Get device */
	HANDLE_CUDA_ERROR(cudaSetDevice(0));

	/* Allocate and copy memory */
	float *dev_equal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_equal, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_equal, equal, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_lower;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_lower, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_lower, lower, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_diagonal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_diagonal, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_diagonal, diagonal, dim * sizeof(float), cudaMemcpyHostToDevice));
	
	float *dev_upper;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_upper, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_upper, upper, dim * sizeof(float), cudaMemcpyHostToDevice));

	/* Run kernel */
	if (dim > max_threads_per_block)
	{
		printf("Thead count (%i) exceeds maximum", dim);
		return;
	}
	cyclic_reduction_device<<<1, dim>>>(dev_lower, dev_diagonal, dev_upper, dev_equal, dim);

	/* Copy data back */
	HANDLE_CUDA_ERROR(cudaDeviceSynchronize());
	HANDLE_CUDA_ERROR(cudaMemcpy(equal, dev_equal, dim * sizeof(float), cudaMemcpyDeviceToHost));

	/* Clean up */
	HANDLE_CUDA_ERROR(cudaFree(dev_equal));
	HANDLE_CUDA_ERROR(cudaFree(dev_lower));
	HANDLE_CUDA_ERROR(cudaFree(dev_diagonal));
	HANDLE_CUDA_ERROR(cudaFree(dev_upper));

	/* Flush profiling info */
	HANDLE_CUDA_ERROR(cudaDeviceReset());
}


__global__ void parallel_cyclic_reduction_device(float *lower_glb, float *diagonal_glb, float *upper_glb, float *equal_glb, const int dim)
{
	const int rank = threadIdx.x;
	__shared__ float lower[512];
	__shared__ float diagonal[512];
	__shared__ float upper[512];
	__shared__ float equal[512];

	lower[rank] = lower_glb[rank];
	diagonal[rank] = diagonal_glb[rank];
	upper[rank] = upper_glb[rank];
	equal[rank] = equal_glb[rank];
	__syncthreads();

	float lower_tmp;
	float upper_tmp;
	float result_tmp;
	float diag_tmp;
	for (int span = 1 ; span < dim; span <<= 1)
	{
		if (rank < dim)
		{
			result_tmp = equal[rank];
			diag_tmp = diagonal[rank];

			if (rank - span >= 0)
			{
				lower_tmp = -lower[rank] / diagonal[rank - span];
				diag_tmp += lower_tmp * upper[rank - span];
				result_tmp += lower_tmp * equal[rank - span];
				lower_tmp *= lower[rank - span];
			}

			if (rank + span < dim)
			{
				upper_tmp = -upper[rank] / diagonal[rank + span];
				diag_tmp += upper_tmp * lower[rank + span];
				result_tmp += upper_tmp * equal[rank + span];
				upper_tmp *= upper[rank + span];
			}
		}
		__syncthreads();

		if (rank < dim)
		{
			lower[rank] = lower_tmp;
			upper[rank] = upper_tmp;
			equal[rank] = result_tmp;
			diagonal[rank] = diag_tmp;
		}
		__syncthreads();
	}

	if (rank < dim)
	{
		equal_glb[rank] = equal[rank] / diagonal[rank];
	}
}


void parallel_cyclic_reduction(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	/* Get device */
	HANDLE_CUDA_ERROR(cudaSetDevice(0));

	/* Allocate and copy memory */
	float *dev_equal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_equal, dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_equal, equal, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_lower;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_lower, dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_lower, lower, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_diagonal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_diagonal, dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_diagonal, diagonal, dim * sizeof(float), cudaMemcpyHostToDevice));
	
	float *dev_upper;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_upper, dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_upper, upper, dim * sizeof(float), cudaMemcpyHostToDevice));

	/* Run kernel */
	if (dim > max_threads_per_block)
	{
		printf("Thead count (%i) exceeds maximum", dim);
		return;
	}
	parallel_cyclic_reduction_device<<<1, dim>>>(dev_lower, dev_diagonal, dev_upper, dev_equal, dim);

	/* Copy data back */
	HANDLE_CUDA_ERROR(cudaDeviceSynchronize());
	HANDLE_CUDA_ERROR(cudaMemcpy(equal, dev_equal, dim * sizeof(float), cudaMemcpyDeviceToHost));

	/* Clean up */
	HANDLE_CUDA_ERROR(cudaFree(dev_equal));
	HANDLE_CUDA_ERROR(cudaFree(dev_lower));
	HANDLE_CUDA_ERROR(cudaFree(dev_diagonal));
	HANDLE_CUDA_ERROR(cudaFree(dev_upper));

	/* Flush profiling info */
	HANDLE_CUDA_ERROR(cudaDeviceReset());
}


__global__ void hybrid_cyclic_reduction_device(float *lower_glb, float *diagonal_glb, float *upper_glb, float *equal_glb, const int dim)
{
	__shared__ float lower[512];
	__shared__ float diagonal[512];
	__shared__ float upper[512];
	__shared__ float equal[512];

	lower[threadIdx.x] = lower_glb[threadIdx.x];
	diagonal[threadIdx.x] = diagonal_glb[threadIdx.x];
	upper[threadIdx.x] = upper_glb[threadIdx.x];
	equal[threadIdx.x] = equal_glb[threadIdx.x];
	__syncthreads();


	/* Cyclic forward reduction */
	int step = cyclic_reduction_forward_reduction(lower, diagonal, upper, equal, dim, 1, 128);

	/* Parallel cyclic reduction to solve system */
	float lower_tmp;
	float upper_tmp;
	float result_tmp;
	float diag_tmp;
	const int rank = (threadIdx.x * step) + step - 1;
	for (int span = step; span < dim; span <<= 1)
	{
		if (rank < dim)
		{
			result_tmp = equal[rank];
			diag_tmp = diagonal[rank];

			if (rank - span >= 0)
			{
				lower_tmp = -lower[rank] / diagonal[rank - span];
				diag_tmp += lower_tmp * upper[rank - span];
				result_tmp += lower_tmp * equal[rank - span];
				lower_tmp *= lower[rank - span];
			}

			if (rank + span < dim)
			{
				upper_tmp = -upper[rank] / diagonal[rank + span];
				diag_tmp += upper_tmp * lower[rank + span];
				result_tmp += upper_tmp * equal[rank + span];
				upper_tmp *= upper[rank + span];
			}
		}
		__syncthreads();

		if (rank < dim)
		{
			lower[rank] = lower_tmp;
			upper[rank] = upper_tmp;
			equal[rank] = result_tmp;
			diagonal[rank] = diag_tmp;
		}
		__syncthreads();
	}

	if (rank < dim)
	{
		equal[rank] /= diagonal[rank];
	}
	__syncthreads();


	/* Cyclic backward substitution */
	cyclic_reduction_back_substitution(lower, diagonal, upper, equal, dim, step >> 1, 0);
	equal_glb[threadIdx.x] = equal[threadIdx.x];
}


void hybrid_cyclic_reduction(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	const int log_dim = static_cast<int>(ceil(log(static_cast<float>(dim)) / log(2.0f)));

	/* Get device */
	HANDLE_CUDA_ERROR(cudaSetDevice(0));

	/* Allocate and copy memory */
	float *dev_equal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_equal, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_equal, equal, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_lower;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_lower, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_lower, lower, dim * sizeof(float), cudaMemcpyHostToDevice));

	float *dev_diagonal;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_diagonal, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_diagonal, diagonal, dim * sizeof(float), cudaMemcpyHostToDevice));
	
	float *dev_upper;
	HANDLE_CUDA_ERROR(cudaMalloc((void **)&dev_upper, dim * log_dim * sizeof(float)));
	HANDLE_CUDA_ERROR(cudaMemcpy(dev_upper, upper, dim * sizeof(float), cudaMemcpyHostToDevice));

	/* Run kernel */
	if (dim > max_threads_per_block)
	{
		printf("Thead count (%i) exceeds maximum", dim);
		return;
	}
	hybrid_cyclic_reduction_device<<<1, dim>>>(dev_lower, dev_diagonal, dev_upper, dev_equal, dim);

	/* Copy data back */
	HANDLE_CUDA_ERROR(cudaDeviceSynchronize());
	HANDLE_CUDA_ERROR(cudaMemcpy(equal, dev_equal, dim * sizeof(float), cudaMemcpyDeviceToHost));

	/* Clean up */
	HANDLE_CUDA_ERROR(cudaFree(dev_equal));
	HANDLE_CUDA_ERROR(cudaFree(dev_lower));
	HANDLE_CUDA_ERROR(cudaFree(dev_diagonal));
	HANDLE_CUDA_ERROR(cudaFree(dev_upper));

	/* Flush profiling info */
	HANDLE_CUDA_ERROR(cudaDeviceReset());
}
