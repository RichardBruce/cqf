
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <windows.h>



void serialFir(float *filtered, const float *values, const float *coeffs, unsigned int size, unsigned int taps)
{
	const time_t start_time = time(0);
	fprintf(stdout, "Serial start time: %d\n", start_time);

	// For each value
	for (unsigned int i = 0; i < size - taps; ++i)
	{
		// MAC accross tap coeffs
		float res = 0.f;
		for (unsigned int j = 0; j < taps; ++j)
		{
			const float v = values[i + j];
			const float c = coeffs[j];
			res += v * c;
		}
		filtered[i] = res;
	}

	fprintf(stdout, "Serial end time: %d\n", time(0));
	fprintf(stdout, "Serial run time: %d\n", time(0) - start_time);
}


// Parallel fir kernel
__global__ void firKernel(float *filtered, const float *values, const float *coeffs, const int size, const int taps)
{
	// Copy the coeffs
	extern __shared__ float s_coeffs[];
	const int coeff_blocks = (int)((taps / (float)blockDim.x) + 0.5f);
	for (unsigned int i = 0; i < (coeff_blocks * blockDim.x); i += blockDim.x)
	{
		if (i < taps)
		{
			s_coeffs[i] = coeffs[i];
		}
	}

	// Wait for all threads
	__syncthreads();

    int i = blockIdx.x * blockDim.x + threadIdx.x;
	//int i = threadIdx.x;

	// MAC accross tap coeffs
	float res = 0.f;
	for (unsigned int j = 0; j < taps; ++j)
	{
		res += values[i + j] * s_coeffs[j];
	}
	//filtered[(blockIdx.x * blockDim.x) + i] = res;
	filtered[i] = res;
}


// Helper function for using CUDA to add vectors in parallel.
cudaError_t firWithCuda(float *filtered, const float *values, const float *coeffs, unsigned int size, unsigned int taps)
{
    float *dev_values = 0;
    float *dev_coeffs = 0;
    float *dev_filtered = 0;
    cudaError_t cudaStatus;
	
    const int threadCount = 64;
    const int blockCount = ((size - taps) / threadCount) + ((size - taps) % threadCount == 0?0:1); 


	const time_t start_time = time(0);
	fprintf(stdout, "Parallel start time: %d\n", start_time);

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed: %s", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_values, size * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_coeffs, taps * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_filtered, (size - taps) * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_values, values, size * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_coeffs, coeffs, taps * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    firKernel<<<blockCount, threadCount, taps * sizeof(float)>>>(dev_filtered, dev_values, dev_coeffs, size, taps);

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error: %s!\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(filtered, dev_filtered, (size - taps) * sizeof(float), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	fprintf(stdout, "Parallel end time: %d\n", time(0));
	fprintf(stdout, "Parallel run time: %d\n", time(0) - start_time);

Error:
    cudaFree(dev_values);
    cudaFree(dev_coeffs);
    cudaFree(dev_filtered);
    
    return cudaStatus;
}


/*int main()
{
	const float EPSILON = 1.0e-3f;

    const int size = 1048576;//131072;//4194304;
	const int taps = 1024;//8192;//8192;
    float *values;
    float *coeffs;
    float *serial_filtered;
	float *parallel_filtered;

	values = (float *)malloc(size * sizeof(float));
	coeffs = (float *)malloc(taps * sizeof(float));
	serial_filtered = (float *)malloc((size - taps) * sizeof(float));
	parallel_filtered = (float *)malloc((size - taps) * sizeof(float));

	// Initialise data
	for (unsigned int i = 0; i < size; ++i)
	{
		values[i] = (float)sin(2. * M_PI * (i / (double)size));
	}

	for (unsigned int i = 0; i < taps; ++i)
	{
		coeffs[i] = (float)cos(2. * M_PI * (i / (double)taps));
	}

	// Run serial kernel
	serialFir(serial_filtered, values, coeffs, size, taps);

    // Run parallel fir
    cudaError_t cudaStatus = firWithCuda(parallel_filtered, values, coeffs, size, taps);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "firWithCuda failed!");
        return 1;
    }

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

	// Check output
	for (unsigned int i = 0; i < size - taps; ++i)
	{
		if (fabs(serial_filtered[i] - parallel_filtered[i]) > EPSILON)
		{
			fprintf(stdout, "Mismatch at %i, expected: %f, got %f\n", i, serial_filtered[i], parallel_filtered[i]);
		}
	}

	// Clean up
	free(values);
	free(coeffs);
	free(serial_filtered);
	free(parallel_filtered);

    return 0;
}*/
