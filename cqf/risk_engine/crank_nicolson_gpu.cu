#include <algorithm>

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

/* Shared space partition */
const int matrix_lower_pos      = 0;
const int matrix_mid_pos        = matrix_lower_pos + max_ns;
const int matrix_upper_pos      = matrix_mid_pos + max_ns;



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


__device__ float asian_call_payoff(const float s, const float a, const float k, const int asianings)
{
    const float val = (a / (asianings - 1)) + (s / asianings);
    return fmaxf(0.0f, val - k);
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
__device__ void populate_matrix(float *const scratch, float *const matrix, float *const matrix_equal, const float *const tp1, const float *const grid, 
    const float half_sigma_sq, const float r, const float t_inc, const int ns, const int i)
{
    /* Boundary conditions */
    /* s = 0.0 */
    if (i == 0)
    {
        const float b = -r * 0.5f * t_inc;
        matrix[matrix_mid_pos  ] = 1.0f - b;
        matrix[matrix_upper_pos] = 0.0f;
        
        matrix_equal[0] = (1.0f + b) * tp1[0];
    }
    /* s = s_max*/
    else if (i == (ns - 1))
    {
        const float r_s = r * grid[ns - 1];

        const float a = -r_s * 0.5f * t_inc;
        const float b = -(r - r_s) * 0.5f * t_inc;
        
        matrix[matrix_lower_pos + ns - 1] = -a;
        matrix[matrix_mid_pos   + ns - 1] = 1.0f - b;
        
        matrix_equal[ns - 1]  = a * tp1[ns - 2];
        matrix_equal[ns - 1] += (1.0f + b) * tp1[ns - 1];
    }
    else if (i < ns)
    {
        const float g = half_sigma_sq * grid[i] * grid[i];
        const float r_s = r * grid[i];
            
        const float a = ((scratch[lower_delta_pos + i] * r_s) + (scratch[lower_gamma_pos + i] * g))     * 0.5f * t_inc;
        const float b = ((scratch[mid_delta_pos + i]   * r_s) + (scratch[mid_gamma_pos + i]   * g) - r) * 0.5f * t_inc;
        const float c = ((scratch[upper_delta_pos + i] * r_s) + (scratch[upper_gamma_pos + i] * g))     * 0.5f * t_inc;
            
        matrix[matrix_lower_pos + i] = -a;
        matrix[matrix_mid_pos   + i] = 1.0f - b;
        matrix[matrix_upper_pos + i] = -c;
            
        matrix_equal[i]  = a * tp1[i - 1];
        matrix_equal[i] += (1.0f + b) * tp1[i];
        matrix_equal[i] += c * tp1[i + 1];
    }

    __syncthreads();
}


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
                lower[addr]     = alpha * lower[addr - step];
            }

            if (addr + step < dim)
            {
                const float gamma = -upper[addr] / diagonal[addr + step];
                equal[addr]    += (gamma * equal[addr + step]);
                diagonal[addr] += (gamma * lower[addr + step]);
                upper[addr]     = gamma * upper[addr + step];
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


__device__ void cyclic_reduction_device(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
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
}


__device__ void parallel_cyclic_reduction(float *l, float *d, float *u, float *h, const int ns, const int i)
{
    for (int step = 1; step < ns; step <<= 1)
    {
        float h_tmp = h[i];
        float d_tmp = d[i];

        float l_tmp;
        if (i - step >= 0)
        {
            l_tmp = -l[i] / d[i - step ];

            d_tmp += l_tmp * u[i - step];
            h_tmp += l_tmp * h[i - step];
            l_tmp *= l[i - step];
        }

        float u_tmp;
        if (i + step < ns)
        {
            u_tmp = -u[i] / d[i + step];

            d_tmp += u_tmp * l[i + step];
            h_tmp += u_tmp * h[i + step];
            u_tmp *= u[i + step];
        }
        __syncthreads();

        l[i] = l_tmp;
        u[i] = u_tmp;
        h[i] = h_tmp;
        d[i] = d_tmp;
        __syncthreads();
    }

    h[i] /= d[i];
    __syncthreads();
}


__device__ void solve_tridiagonal(float *const matrix, float *const matrix_equal, const int ns, const int i)
{
    //cyclic_reduction_device(&matrix[matrix_lower_pos], &matrix[matrix_mid_pos], &matrix[matrix_upper_pos], matrix_equal, ns);
    parallel_cyclic_reduction(&matrix[matrix_lower_pos], &matrix[matrix_mid_pos], &matrix[matrix_upper_pos], matrix_equal, ns, i);
}


__global__ void crank_nicolson(const float *const grid, float *const scratch, const float half_sigma_sq, const float r, 
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
    __shared__ float shared_matrix[3 * max_ns];
    for (unsigned int j = 0; j < nt >> 1; ++j)
    {
        populate_matrix(scratch, shared_matrix, shared_equal, shared_tp1, grid, half_sigma_sq, r, t_inc, ns, i);
        solve_tridiagonal(shared_matrix, shared_equal, ns, i);
        shared_equal[i] = fmaxf(shared_equal[i], call_payoff(shared_equal[i], k));
        __syncthreads();

        populate_matrix(scratch, shared_matrix, shared_tp1, shared_equal, grid, half_sigma_sq, r, t_inc, ns, i);
        solve_tridiagonal(shared_matrix, shared_tp1, ns, i);
        shared_tp1[i] = fmaxf(shared_tp1[i], call_payoff(shared_tp1[i], k));
        __syncthreads();
    }

    scratch[matrix_equal_pos + i] = shared_tp1[i];
}


__global__  void transpose(float *const trans, const int x, const int y)
{
    /* Position of this thread in a tranpose block */
    const int blk_size = 32;
    const int x_offset = threadIdx.x >> 5;
    const int y_offset = threadIdx.x & (blk_size - 1);

    const int x_blks = x & ~(blk_size - 1);
    const int y_blks = y & ~(blk_size - 1);

    __shared__ float shared_trans[blk_size][blk_size + 1][2]; /* 32x32 block of transposed data */

    /* Transpose off diagonal blocks */
    for (int i = blk_size; i < x_blks; i += blk_size)
    {
        for (int j = 0; j < i; j += blk_size)
        {
            /* Streaming load and within block transposed save into shared */
            shared_trans[y_offset][x_offset][0] = trans[((i + x_offset) * y) + j + y_offset];
            shared_trans[y_offset][x_offset][1] = trans[((j + x_offset) * y) + i + y_offset];

            __syncthreads();

            /* Streaming save */
            trans[((i + x_offset) * y) + j + y_offset] = shared_trans[x_offset][y_offset][1];
            trans[((j + x_offset) * y) + i + y_offset] = shared_trans[x_offset][y_offset][0];
        }
    }

    /* Transpose diagonal blocks */
    for (int i = 0; i < min(x_blks, y_blks); i += blk_size)
    {
        /* Streaming load and within block transposed save into shared */
        const int blk_addr = ((i + x_offset) * y) + i + y_offset;
        shared_trans[y_offset][x_offset][0] = trans[blk_addr];

        __syncthreads();

        /* Streaming save */
        trans[blk_addr] = shared_trans[x_offset][y_offset][0];
    }
}


void transpose_test()
{
    const int x = 1024;
    const int y = 1024;
    float *trans = new float [x * y];
    for (int i = 0; i < x * y; ++i)
    {
        trans[i] = i;
    }

    print_cuda_error(cudaSetDevice(0), "Set device");

    /* Prepare device memory */
    float *dev_trans;
    print_cuda_error(cudaMalloc((void **)&dev_trans, x * y * sizeof(float)), "Malloc matrix");
    print_cuda_error(cudaMemcpy(dev_trans, trans, x * y * sizeof(float), cudaMemcpyHostToDevice), "Copy martix to device");

    /* Call kernel */
    cudaProfilerStart();
    transpose<<<1, 1024>>>(dev_trans, x, y);
    cudaProfilerStop();
    print_cuda_error(cudaGetLastError(), "Kernel execution");

    print_cuda_error(cudaMemcpy(trans, dev_trans, x * y * sizeof(float), cudaMemcpyDeviceToHost), "Copy matrix to host");
    // for (int i = 0; i < x; ++i)
    // {
    //     for (int j = 0; j < y; ++j)
    //     {
    //         printf("%.2f ", trans[(i * y) + j]);
    //     }
    //     printf("\n");
    // }

    /* Clean up */
    print_cuda_error(cudaFree(dev_trans), "Free matrix");
    
    print_cuda_error(cudaDeviceReset(), "Device reset");

    delete [] trans;
}


__device__ __host__ float interpolate(const float *const x, const float *const y, const float a, const int ns, const int a_idx, const int s_idx)
{
    /* Linear interpolation */
    const int y_idx = (a_idx * ns) + s_idx;
    return y[y_idx - ns] + ((y[y_idx] - y[y_idx - ns]) * ((a - x[a_idx - 1]) / (x[a_idx] - x[a_idx - 1])));
}


// __global__ void permute_asian_values(float *const v_grid1, const float *const a_grid1, const float *const a_grid0, const float *const s_grid, 
//     const float *const v_grid0, const float a_fac, const float s_fac, const float a_inc_inv, const int ns, const int na0, const int na1)
// {
//     const int y_lower = blockIdx.y * blockDim.x;
//     const int y_upper = (blockIdx.y * blockDim.x) + blockDim.x;

//     const int i = (blockIdx.x * blockDim.x) + threadIdx.x;

//     int cache_idx = threadIdx.x;
//     __shared__ float cache[64 * 64];
//     for (int j = y_lower; j < y_upper; ++j)
//     {
//         cache[cache_idx] = v_grid0[(j * ns) + i];
//         cache_idx += blockDim.x;
//     }

//     /* I would fix at s_value */
//     const float s_value = s_grid[i];

//     int a1_idx = 0;
//     const float prefix_asian_value = a_grid1[a1_idx];

//     /* After fixing I was at postfix_asian_value */
//     const float postfix_asian_value = (prefix_asian_value * a_fac) + (s_value * s_fac);
    
//     /* Find and interpolate around the postfix_asian_value */
//     int a0_idx = max(y_lower + 1, min(y_upper, static_cast<int>(postfix_asian_value * a_inc_inv)));
//     for (int j = y_lower; j < y_upper; ++j)
//     {
//         while (a0_idx == j)
//         {
//             v_grid1[(a1_idx * ns) + i] = interpolate(a_grid0, cache, postfix_asian_value, blockDim.x, a0_idx - y_lower, threadIdx.x);

//             const float prefix_asian_value = a_grid1[++a1_idx];
//             const float postfix_asian_value = (prefix_asian_value * a_fac) + (s_value * s_fac);
//             a0_idx = max(y_lower + 1, min(y_upper, static_cast<int>(postfix_asian_value * a_inc_inv)));
//         }
//     }
// }


__global__ void permute_asian_values(float *const v_grid1, const float *const a_grid1, const float *const a_grid0, const float *const s_grid, 
    const float *const v_grid0, const float a_fac, const float s_fac, const float a_inc_inv, const int ns, const int na0, const int na1)
{
    const int i = threadIdx.x;

    /* I would fix at s_value */
    const float s_value = s_grid[i];

    int a1_idx = 0;
    const float prefix_asian_value = a_grid1[a1_idx];

    /* After fixing I was at postfix_asian_value */
    const float postfix_asian_value = (prefix_asian_value * a_fac) + (s_value * s_fac);
    
    /* Find and interpolate around the postfix_asian_value */
    int a0_idx = max(1, min(na0, static_cast<int>(postfix_asian_value * a_inc_inv)));
    for (int j = 0; j < na0; ++j)
    {
        while (a0_idx == j)
        {
            v_grid1[(a1_idx * ns) + i] = interpolate(a_grid0, v_grid0, postfix_asian_value, ns, a0_idx, i);

            const float prefix_asian_value = a_grid1[++a1_idx];
            const float postfix_asian_value = (prefix_asian_value * a_fac) + (s_value * s_fac);
            a0_idx = max(1, min(na0, static_cast<int>(postfix_asian_value * a_inc_inv)));
        }
    }
}


__host__ int update_asian_grid(const float *const grid_ping, float *const grid_pong, const int na, const int t_idx)
{
    const int t_idx_m1          = t_idx - 1;
    const float flt_t_idx       = static_cast<float>(t_idx);
    const float t_idx_inv       = 1.0f / flt_t_idx;
    const float t_idx_m1_inv    = 1.0f / (flt_t_idx - 1.0f);
    const float a_fac           = (flt_t_idx - 1.0f) * t_idx_inv;

    /* Asianing update */
    /* Rebuild uniform asian grid */
    /* This may be slightly too big, but not too much to worry about */
    /* This loop will get all, but the last t_idx - 1 or less points */
    const unsigned int whole_iters = (na - 1) / t_idx;
    for (unsigned int i = 0; i < whole_iters; i++)
    {
        const unsigned int l_idx = i * t_idx;
        const unsigned int h_idx = l_idx + t_idx;
        const float step = (grid_ping[h_idx] - grid_ping[l_idx]) * t_idx_m1_inv;
        for (unsigned int j = 0; j < t_idx_m1; j++)
        {
            grid_pong[(i * t_idx_m1) + j] = grid_ping[l_idx] + (step * static_cast<float>(j));
        }
    }
    
    /* This loop will get the remaining points */
    const unsigned int part_iters = na - (whole_iters * t_idx);
    const unsigned int l_idx = na - part_iters - 1;
    const float step = (grid_ping[na - 1] - grid_ping[l_idx]) * t_idx_m1_inv;
    for (unsigned int i = 0; i < part_iters; i++)
    {
        grid_pong[(whole_iters * t_idx_m1) + i] = grid_ping[l_idx] + (step * static_cast<float>(i));
    }
    grid_pong[(whole_iters * t_idx_m1) + part_iters - 1] = grid_ping[na - 1];

    return (na * a_fac) + 1;
}


__global__ void crank_nicolson_asian(const float *const grid, float *const glb_values, float *const glb_scratch, const float half_sigma_sq, 
    const float r, const float t_inc, const float k, const int ns, const int nt, const int asianings, const bool create_payoff)
{
    const int i = threadIdx.x;
    const int asian_idx = blockIdx.x;
    float *const scratch = &glb_scratch[asian_idx * scratch_space_size];
    float *const values = &glb_values[asian_idx * ns];
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
    if (create_payoff)
    {
        shared_tp1[i] = asian_call_payoff(shared_equal[i], grid[ns + asian_idx], k, asianings);
    }
    else
    {
        shared_tp1[i] = values[i];
    }
    get_coeffs(shared_equal, scratch, ns, i);
    
    /* Solve back in time */
    __shared__ float shared_matrix[3 * max_ns];
    for (unsigned int j = 0; j < nt >> 1; ++j)
    {
        populate_matrix(scratch, shared_matrix, shared_equal, shared_tp1, grid, half_sigma_sq, r, t_inc, ns, i);
        solve_tridiagonal(shared_matrix, shared_equal, ns, i);
        shared_equal[i] = fmaxf(shared_equal[i], call_payoff(shared_equal[i], k));
        __syncthreads();

        populate_matrix(scratch, shared_matrix, shared_tp1, shared_equal, grid, half_sigma_sq, r, t_inc, ns, i);
        solve_tridiagonal(shared_matrix, shared_tp1, ns, i);
        shared_tp1[i] = fmaxf(shared_tp1[i], call_payoff(shared_tp1[i], k));
        __syncthreads();
    }

    values[i] = shared_tp1[i];
}


void asian_call_test()
{
    /* Pricing set up */
    printf("Pricing Asian Call\n");
    const unsigned int ns = 1024;   /* Want multiples of warp size (32) */

    const int asianings = 2;
    const float t_inc = 0.01;
    const float t[asianings] = { 0.9f, 1.0f };
    const float k = 100.0f;

    const float s = 100.0f;
    const float r = 0.05f;
    const float sigma = 0.2f;
    const float half_sigma_sq = 0.5f * sigma * sigma;

    /* Build regular grid based at 0 */
    int na = 1024;
    const int grid_size = ns + na;
    float *grid_ping = new float [grid_size];

    const float s_inc = (s * 3.0f) / ns;
    for (int i = 0; i < ns; ++i)
    {
        grid_ping[i] = i * s_inc;
    }

    const float a_inc = (s * 3.0f) / na;
    for (int i = 0; i < na; ++i)
    {
        grid_ping[ns + i] = i * a_inc;
    }

    float *grid_pong = new float [grid_size];
    memcpy(grid_pong, grid_ping, ns * sizeof(float));


    /* Set cuda device */
    print_cuda_error(cudaSetDevice(0), "Set device");

    /* Prepare device memory */
    float *dev_grid;
    print_cuda_error(cudaMalloc((void **)&dev_grid, 2 * grid_size * sizeof(float)), "Malloc grid");
    print_cuda_error(cudaMemcpy(dev_grid, grid_ping, grid_size * sizeof(float), cudaMemcpyHostToDevice), "Copy grid to device");
    float *dev_grid_ping = &dev_grid[0];
    float *dev_grid_pong = &dev_grid[grid_size];

    float *values = new float [ns * na];

    float *dev_values;
    print_cuda_error(cudaMalloc((void **)&dev_values, 2 * ns * na * sizeof(float)), "Malloc values");
    print_cuda_error(cudaMemset(dev_values, 0, ns * na * sizeof(float)), "Clear values");
    float *dev_values_ping = &dev_values[0];
    float *dev_values_pong = &dev_values[ns * na];

    float *dev_scratch;
    print_cuda_error(cudaMalloc((void **)&dev_scratch, na * scratch_space_size * sizeof(float)), "Malloc scratch");
    print_cuda_error(cudaMemset(dev_scratch, 0, na * scratch_space_size * sizeof(float)), "Clear scratch");

    /* Call kernel and update asian grid in parallel */
    cudaProfilerStart();

    /* Work back through the asianings */
    for (int i = asianings; i > 1; --i)
    {
        /* Time step the asian slices */
        const int nt = (t[asianings - 1] - t[asianings - 2]) / t_inc;
        crank_nicolson_asian<<<na, ns>>>(dev_grid_ping, dev_values_ping, dev_scratch, half_sigma_sq, r, t_inc, k, ns, nt, i, (i == asianings));
        
        /* Update the asian grid in parallel */
        const float flt_t_idx   = static_cast<float>(i);
        const float t_idx_inv   = 1.0f / flt_t_idx;
        const float a_fac       = (flt_t_idx - 1.0f) * t_idx_inv;
        const int na_lst        = na;
        na = update_asian_grid(&grid_ping[ns], &grid_pong[ns], na, i);
        std::swap(grid_ping, grid_pong);
        std::swap(dev_grid_ping, dev_grid_pong);
        std::swap(dev_values_ping, dev_values_pong);
        print_cuda_error(cudaMemcpy(dev_grid_ping, grid_ping, grid_size * sizeof(float), cudaMemcpyHostToDevice), "Copy grid to device");

        /* Wait for the device to complete */
        print_cuda_error(cudaDeviceSynchronize(), "Synchronise device");
        // print_cuda_error(cudaMemcpy(values, dev_values_pong, ns * na_lst * sizeof(float), cudaMemcpyDeviceToHost), "Copy values to host");
        // printf("      ");
        // for (int i = ns; i < ns + na_lst; ++i)
        // {
        //     printf("%.2f ", grid_pong[i]);
        // }
        // printf("\n");

        // for (unsigned int i = 0; i < ns; ++i)
        // {
        //     printf("%.2f: ", grid_ping[i]);
        //     for (int j = 0; j < na_lst; ++j)
        //     {
        //         printf("%.2f ", values[(j * ns)+ i]);
        //     }
        //     printf("\n");
        // }
        // printf("\n\n");

        /* Update the grid values for asianing */
        dim3 grid(ns >> 6, na >> 6);
        const int a_inc_inv = na_lst / (s * 3.0f);
        permute_asian_values<<<1, ns>>>(dev_values_ping, &dev_grid_ping[ns], &dev_grid_pong[ns], dev_grid_pong, dev_values_pong, 
            a_fac, t_idx_inv, a_inc_inv, ns, na_lst, na);

        /* Wait for the device to complete */
        print_cuda_error(cudaDeviceSynchronize(), "Synchronise device");
        // print_cuda_error(cudaMemcpy(values, dev_values_ping, ns * na * sizeof(float), cudaMemcpyDeviceToHost), "Copy values to host");
        // printf("      ");
        // for (int i = ns; i < ns + na; ++i)
        // {
        //     printf("%.2f ", grid_ping[i]);
        // }
        // printf("\n");

        // for (unsigned int i = 0; i < ns; ++i)
        // {
        //     printf("%.2f: ", grid_ping[i]);
        //     for (int j = 0; j < na; ++j)
        //     {
        //         printf("%.2f ", values[(j * ns)+ i]);
        //     }
        //     printf("\n");
        // }
        // printf("\n\n");
    }

    /* Final period */
    const int nt = t[0] / t_inc;
    crank_nicolson_asian<<<1, ns>>>(dev_grid_ping, dev_values_ping, dev_scratch, half_sigma_sq, r, t_inc, k, ns, nt, 0, false);

    cudaProfilerStop();
    print_cuda_error(cudaGetLastError(), "Kernel execution");


    float *res = new float [ns];
    print_cuda_error(cudaMemcpy(res, &dev_values_ping[0], ns * sizeof(float), cudaMemcpyDeviceToHost), "Copy grid to host");
    for (unsigned int i = 0; i < ns; ++i)
    {
        printf("%.2f: %.2f\n", grid_ping[i], res[i]);
    }

    /* Clean up */
    print_cuda_error(cudaFree(dev_grid), "Free grid");
    print_cuda_error(cudaFree(dev_scratch), "Free scratch");

    print_cuda_error(cudaDeviceReset(), "Device reset");

    delete [] grid_ping;
    delete [] grid_pong;
    delete [] res;
    delete [] values;
}


void american_call_test()
{
    /* Pricing set up */
    printf("Pricing American Call\n");
    const unsigned int ns = 1024;   /* Want multiples of warp size (32) */
    const unsigned int nt = 100;

    const float k = 100.0f;
    const float t = 1.0f;
    const float t_inc = t / nt;

    const float s = 100.0f;
    const float r = 0.05f;
    const float sigma = 0.2f;
    const float half_sigma_sq = 0.5f * sigma * sigma;

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
    crank_nicolson<<<1, ns>>>(dev_grid, dev_scratch, half_sigma_sq, r, t_inc, k, ns, nt);
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
    //asian_call_test();
    //transpose_test();

    return 0;
}


/* Return the first index not less than a */
__device__ __host__ int search(const float *const x, const float a, int i, const int s)
{
    if (x[i] < a)
    {
        while ((i < (s - 1)) && (x[i] < a))
        {
            ++i;
        }
    }
    /* Values are correlates and down, begin linear search downwards */
    else
    {
        while ((i > 1) && (x[i - 1] >= a))
        {
            --i;
        }
    }
    
    return i;
}
