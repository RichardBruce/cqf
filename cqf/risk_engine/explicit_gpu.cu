#include <stdio.h>

#include "cuda_profiler_api.h"

/* Scratch space positions */
const int max_ns    = 1024;

const int lower_delta_pos       = 0;
const int mid_delta_pos         = lower_delta_pos + max_ns;
const int upper_delta_pos       = mid_delta_pos + max_ns;

const int lower_gamma_pos       = upper_delta_pos + max_ns;
const int mid_gamma_pos         = lower_gamma_pos + max_ns;
const int upper_gamma_pos       = mid_gamma_pos + max_ns;

const int matrix_equal_pos      = upper_gamma_pos + max_ns;

const int scratch_space_size    = matrix_equal_pos + max_ns;


void print_cuda_error(cudaError_t err, char *at)
{
    if (err)
    {
        printf("Error from CUDA at : %s\n", at);
        printf("Message: %s\n", cudaGetErrorString(err));
    }
}


__device__ float call_payoff(const float s, const float k)
{
    return fmaxf(0.0f, s - k);
}


__device__ void get_coeffs(const float *const grid, float *const scratch, const int ns, const int i)
{
    /* Difference vs. the grid below */
    float d0;
    float d1;
    if (i == 0)
    {
        d0 = grid[1] - grid[0];
        d1 = grid[2] - grid[1];
    }
    else if (i == (ns - 1))
    {
        d0 = grid[i - 1] - grid[i - 2];
        d1 = grid[i] - grid[i - 1];
    }
    else
    {
        d0 = grid[i] - grid[i - 1];
        d1 = grid[i + 1] - grid[i];
    }
    const float d1_p_d2 = d0 + d1;

    /* Delta coeffs */
    /* Middle */
    if ((i != 0) & (i != (ns - 1)))
    {
        scratch[lower_delta_pos + i] = -d1 / (d0 * d1_p_d2);
        scratch[mid_delta_pos   + i] = (d1 - d0) / (d0 * d1);
        scratch[upper_delta_pos + i] = d0 / (d1 * d1_p_d2);
    }
    /* Lower boundary */
    else if (i == 0)
    {
        scratch[lower_delta_pos + i] = (-2.0f * d0 - d1) / (d0 * d1_p_d2);
        scratch[mid_delta_pos   + i] = d1_p_d2 / (d0 * d1);
        scratch[upper_delta_pos + i] = -d0 / (d1 * d1_p_d2);
    }
    /* Upper boundary */
    else if (i == (ns - 1))
    {
        scratch[lower_delta_pos + i] = d1 / (d0 * d1_p_d2);
        scratch[mid_delta_pos   + i] = (-d0 - d1) / (d0 * d1);
        scratch[upper_delta_pos + i] = (d0 + 2.0f * d1) / (d1 * d1_p_d2);
    }

    /* Gamma coeffs */
    /* Middle */
    if ((i != 0) & (i != (ns - 1)))
    {
        scratch[lower_gamma_pos + i]  =  2.0f / (d0 * d1_p_d2);
        scratch[mid_gamma_pos   + i]  = -2.0f / (d0 * d1);
        scratch[upper_gamma_pos + i]  =  2.0f / (d1 * d1_p_d2);
    }
    __syncthreads();
}


/* Populate the matrix */
__device__ void explicit_step(float *const scratch, float *const matrix_equal, const float *const tp1, const float *const grid, 
    const float half_sigma_sq, const float r, const float t_inc, const int ns, const int i)
{
    /* Boundary conditions */
    /* s = 0.0 */
    if (i == 0)
    {
        const float b = -r * t_inc;
        matrix_equal[0] = (1.0f + b) * tp1[0];
    }
    /* s = s_max*/
    else if (i == (ns - 1))
    {
        const float r_s = r * grid[ns - 1];

        const float a = -r_s * t_inc;
        const float b = -(r - r_s) * t_inc;
        matrix_equal[ns - 1]  = a * tp1[ns - 2];
        matrix_equal[ns - 1] += (1.0f + b) * tp1[ns - 1];
    }
    else if (i < ns)
    {
        const float g = half_sigma_sq * grid[i] * grid[i];
        const float r_s = r * grid[i];
            
        const float a = ((scratch[lower_delta_pos + i] * r_s) + (scratch[lower_gamma_pos + i] * g))     * t_inc;
        const float b = ((scratch[mid_delta_pos + i]   * r_s) + (scratch[mid_gamma_pos + i]   * g) - r) * t_inc;
        const float c = ((scratch[upper_delta_pos + i] * r_s) + (scratch[upper_gamma_pos + i] * g))     * t_inc;
        matrix_equal[i]  = a * tp1[i - 1];
        matrix_equal[i] += (1.0f + b) * tp1[i];
        matrix_equal[i] += c * tp1[i + 1];
    }

    __syncthreads();
}


__global__ void explicit_scheme(const float *const grid, float *const scratch, const float half_sigma_sq, const float r, 
    const float t_inc, const float k, const int ns, const int nt)
{
    const int i = threadIdx.x;
    if (ns & 0x1f)
    {
        /* Only multiple of 32 space steps are supported */
        return;
    }
    
    /* Move grid to shared memory, needed for off by 1 access and reused */
    __shared__ float shared_equal[max_ns];
    shared_equal[i] = grid[i];
    __syncthreads();

    /* Build grid based coeffs, completely parrallel */
    __shared__ float shared_tp1[max_ns];
    shared_tp1[i] = call_payoff(shared_equal[i], k);
    get_coeffs(shared_equal, scratch, ns, i);
    
    /* Solve back in time */
    for (unsigned int j = 0; j < nt >> 1; ++j)
    {
        explicit_step(scratch, shared_equal, shared_tp1, grid, half_sigma_sq, r, t_inc, ns, i);
        shared_equal[i] = fmaxf(shared_equal[i], call_payoff(shared_equal[i], k));

        explicit_step(scratch, shared_tp1, shared_equal, grid, half_sigma_sq, r, t_inc, ns, i);
        shared_tp1[i] = fmaxf(shared_tp1[i], call_payoff(shared_tp1[i], k));
    }

    scratch[matrix_equal_pos + i] = shared_tp1[i];
}


void american_call_test()
{
    /* Pricing set up */
    printf("Pricing American Call\n");
    const unsigned int ns = 1024;   /* Want multiples of warp size (32) */

    const float s = 100.0f;
    const float r = 0.05f;
    const float sigma = 0.2f;
    const float half_sigma_sq = 0.5f * sigma * sigma;

    const float k = 100.0f;
    const float t = 1.0f;
    const float t_inc   = 0.9f / (static_cast<float>(ns * ns) * sigma * sigma);
    const int nt        = static_cast<int>(t / t_inc) + 1;

    /* Build regular grid based at 0 */
    float *grid = new float [ns];
    const float s_inc = (s * 3.0f) / ns;
    for (unsigned int i = 0; i < ns; ++i)
    {
        grid[i] = i * s_inc;
    }

    print_cuda_error(cudaSetDevice(0), "Set device");

    /* Prepare device memory */
    float *dev_grid;
    print_cuda_error(cudaMalloc((void **)&dev_grid, ns * sizeof(float)), "Malloc grid");
    print_cuda_error(cudaMemcpy(dev_grid, grid, ns * sizeof(float), cudaMemcpyHostToDevice), "Copy grid to device");

    float *dev_scratch;
    print_cuda_error(cudaMalloc((void **)&dev_scratch, scratch_space_size * sizeof(float)), "Malloc scratch");
    print_cuda_error(cudaMemset(dev_scratch, 0, scratch_space_size * sizeof(float)), "Clear scratch");

    /* Call kernels */
    cudaProfilerStart();
    explicit_scheme<<<1, ns>>>(dev_grid, dev_scratch, half_sigma_sq, r, t_inc, k, ns, nt);
    cudaProfilerStop();
    print_cuda_error(cudaGetLastError(), "Kernel execution");

    float *res = new float [ns];
    print_cuda_error(cudaMemcpy(res, &dev_scratch[matrix_equal_pos], ns * sizeof(float), cudaMemcpyDeviceToHost), "Copy grid to host");
    for (unsigned int i = 0; i < ns; ++i)
    {
        printf("%.2f: %.2f\n", grid[i], res[i]);
    }

    /* Clean up */
    print_cuda_error(cudaFree(dev_grid), "Free grid");
    print_cuda_error(cudaFree(dev_scratch), "Free scratch");

    print_cuda_error(cudaDeviceReset(), "Device reset");

    delete [] grid;
    delete [] res;
}


int main()
{
    american_call_test();

    return 0;
}
