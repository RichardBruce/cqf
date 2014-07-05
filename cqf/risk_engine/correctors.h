#ifndef __CORRECTORS_H__
#define __CORRECTORS_H__


#include "matrix.h"

#include "grid.h"


/* Explicit predictor step */
template <class T>
class explicit_predictor
{
    public :
        /* CTOR */
        explicit_predictor(const grid<T> &g, const T theta, const T t_inc,
            const T sigma_s1, const T sigma_s2, const T r_s1, const T r_s2, const T corr)
            : _g(g), _theta(theta), _t_inc(t_inc),
              _sigma_s1(sigma_s1), _sigma_s2(sigma_s2), 
              _half_sigma_s1_sq(0.5 * sigma_s1 * sigma_s1), 
              _half_sigma_s2_sq(0.5 * sigma_s2 * sigma_s2), 
              _r_s1(r_s1), _r_s2(r_s2), _corr(corr)
        { }

        /* DTOR */
        ~explicit_predictor() {  };
        
        /* Implicit correction algorithm of the form F(A) = B + _theta * _t_inc * (F(A) - F(C)). 
           F is the 1D BS PDE. A, B and C may be any data sets to correct from */
        void operator()(const T *const src_offset0, T *const dst_offset)
        {
            T coeffs[6];
            const T s2_min = _g.s_grid(0)[0];
            const T s2_max = _g.s_grid(0)[_g.size(0) - 1];

//const T t_inc           = p.get_time_increment();
            
            /* Fill in the row */
            for (int i = 1; i < _g.size(1) - 1; i++)
            {
                const int s1_idx        = i * _g.size(0);
                const int s1_idx_p1     = s1_idx + _g.size(0);
                const int s1_idx_m1     = s1_idx - _g.size(0);
                const T s1              = _g.s_grid(1)[i];
    
                /* Fill in the column */
                for (int j = 1; j < _g.size(0) - 1; j++)
                {
                    const T s2 = _g.s_grid(0)[j];
                    
                    /* S1 coefficients */
                    const T g1 = _half_sigma_s1_sq * s1 * s1;
                    const T r_s1 = _r_s1 * s1;
                
                    _g.get_delta_coeffs_middle(&coeffs[0], 1, i);
                    _g.get_gamma_coeffs_middle(&coeffs[3], 1, i);
                    const T s1_a = ((coeffs[0] * r_s1) + (coeffs[3] * g1)) * _theta * _t_inc;
                    const T s1_b = ((coeffs[1] * r_s1) + (coeffs[4] * g1) - _r_s1) * _theta * _t_inc;
                    const T s1_c = ((coeffs[2] * r_s1) + (coeffs[5] * g1)) * _theta * _t_inc;
                    
                    const T s1_s = (s1_a * src_offset0[s1_idx_m1 + j]) + 
                                   (s1_b * src_offset0[s1_idx    + j]) + 
                                   (s1_c * src_offset0[s1_idx_p1 + j]);
        
                    /* S2 coefficients */
                    const T g2 = _half_sigma_s2_sq * s2 * s2;
                    const T r_s2 = _r_s2 * s2;
                
                    _g.get_delta_coeffs_middle(&coeffs[0], 0, j);
                    _g.get_gamma_coeffs_middle(&coeffs[3], 0, j);
                    const T s2_a = ((coeffs[0] * r_s2) + (coeffs[3] * g2)) * _theta * _t_inc;
                    const T s2_b = ((coeffs[1] * r_s2) + (coeffs[4] * g2) - _r_s2) * _theta * _t_inc;
                    const T s2_c = ((coeffs[2] * r_s2) + (coeffs[5] * g2)) * _theta * _t_inc;
        
                    const T s2_s = (s2_a * src_offset0[s1_idx + j - 1]) + 
                                   (s2_b * src_offset0[s1_idx + j    ]) + 
                                   (s2_c * src_offset0[s1_idx + j + 1]);
                    
                    /* S1 greeks */
//                    const T delta_s1   = (src_offset0[s1_idx_p1 + j] - src_offset0[s1_idx_m1 + j]) * _s1_inc2_inv;
                    dst_offset[s1_idx + j] = src_offset0[s1_idx + j] + (s1_s + s2_s);
                }
                /* Boundary conditions */
                /* s2 = s2_min */
                {        
                    const std::pair<bool, T> &b = _g.get_value_coeffs_lower();
                    if (b.first)
                    {
                        dst_offset[s1_idx] = b.second;
                    }
                    else
                    {
                        /* S1 coefficients */
                        const T g1 = _half_sigma_s1_sq * s1 * s1;
                        const T r_s1 = _r_s1 * s1;
                
                        _g.get_delta_coeffs_middle(&coeffs[0], 1, i);
                        _g.get_gamma_coeffs_middle(&coeffs[3], 1, i);
                        const T s1_a = ((coeffs[0] * r_s1) + (coeffs[3] * g1)) * _theta * _t_inc;
                        const T s1_b = ((coeffs[1] * r_s1) + (coeffs[4] * g1) - _r_s1) * _theta * _t_inc;
                        const T s1_c = ((coeffs[2] * r_s1) + (coeffs[5] * g1)) * _theta * _t_inc;
                    
                        const T s1_s = (s1_a * src_offset0[s1_idx_m1]) + 
                                       (s1_b * src_offset0[s1_idx   ]) + 
                                       (s1_c * src_offset0[s1_idx_p1]);
        
                        /* S2 coefficients */
                        const T g2 = _half_sigma_s2_sq * s2_min * s2_min;
                        const T r_s2 = _r_s2 * s2_min;
                
                        _g.get_delta_coeffs_lower(&coeffs[0], 0, 0);
                        const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                        const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                        const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;
        
                        const T s2_s = (s2_a * src_offset0[s1_idx    ]) + 
                                       (s2_b * src_offset0[s1_idx + 1]) + 
                                       (s2_c * src_offset0[s1_idx + 2]);
                    
                        /* Update */
                        dst_offset[s1_idx] = src_offset0[s1_idx] + (s1_s + s2_s);
                    }
                }

                /* s2 = s2_max*/
                {
                    const std::pair<bool, T> &b = _g.get_value_coeffs_upper();
                    if (b.first)
                    {
                        dst_offset[s1_idx + _g.size(0) - 1] = b.second;
                    }
                    else
                    {
                        /* S1 coefficients */
                        const T g1 = _half_sigma_s1_sq * s1 * s1;
                        const T r_s1 = _r_s1 * s1;
                
                        _g.get_delta_coeffs_middle(&coeffs[0], 1, i);
                        _g.get_gamma_coeffs_middle(&coeffs[3], 1, i);
                        const T s1_a = ((coeffs[0] * r_s1) + (coeffs[3] * g1)) * _theta * _t_inc;
                        const T s1_b = ((coeffs[1] * r_s1) + (coeffs[4] * g1) - _r_s1) * _theta * _t_inc;
                        const T s1_c = ((coeffs[2] * r_s1) + (coeffs[5] * g1)) * _theta * _t_inc;
                    
                        const T s1_s = (s1_a * src_offset0[s1_idx_m1 + _g.size(0) - 1]) + 
                                       (s1_b * src_offset0[s1_idx    + _g.size(0) - 1]) + 
                                       (s1_c * src_offset0[s1_idx_p1 + _g.size(0) - 1]);
        
                        /* S2 coefficients */
                        const T g2 = _half_sigma_s2_sq * s2_max * s2_max;
                        const T r_s2 = _r_s2 * s2_max;
                
                        _g.get_delta_coeffs_upper(&coeffs[0], 0, 0);
                        const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                        const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                        const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;

                        const T s2_s = (s2_a * src_offset0[s1_idx_m1 + _g.size(0) - 3]) + 
                                       (s2_b * src_offset0[s1_idx    + _g.size(0) - 2]) + 
                                       (s2_c * src_offset0[s1_idx_p1 + _g.size(0) - 1]);
                    
                        /* Update */
                        dst_offset[s1_idx + _g.size(0) - 1] = src_offset0[s1_idx + _g.size(0) - 1] + (s1_s + s2_s);
                    }
                }
            }
    
            /* Upper and lower boundary condition in s1 */
            const int s1_idx_upper  = (_g.size(0) - 1) * _g.size(0);
            const int s1_idx_m1     = s1_idx_upper - _g.size(0);
            const int s1_idx_m2     = s1_idx_m1 - _g.size(0);

            const int s1_idx_lower  = 0;
            const int s1_idx_p1     = _g.size(0);
            const int s1_idx_p2     = (_g.size(0) << 1);
            
            const T s1_min = _g.s_grid(1)[0];
            const T s1_max = _g.s_grid(1)[_g.size(1) - 1];
            for (int i = 1; i < _g.size(0) - 1; i++)
            {
                const T s2 = _g.s_grid(0)[i];

                /* s1 = s1_max */
                const auto &ub = _g.get_value_coeffs_upper();
                if (ub.first)
                {
                    dst_offset[s1_idx_upper + i] = ub.second;
                }
                else
                {
                    /* S1 coefficients */
                    const T g1 = _half_sigma_s1_sq * s1_max * s1_max;
                    const T r_s1 = _r_s1 * s1_max;
                
                    _g.get_delta_coeffs_upper(&coeffs[0], 1, i);
                    const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                    const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                    const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
                    
                    const T s1_s = (s1_a * src_offset0[s1_idx_m2    + i]) + 
                                   (s1_b * src_offset0[s1_idx_m1    + i]) + 
                                   (s1_c * src_offset0[s1_idx_upper + i]);
        
                    /* S2 coefficients */
                    const T g2 = _half_sigma_s2_sq * s2 * s2;
                    const T r_s2 = _r_s2 * s2;
                
                    _g.get_delta_coeffs_middle(&coeffs[0], 0, i);
                    _g.get_gamma_coeffs_middle(&coeffs[3], 0, i);
                    const T s2_a = ((coeffs[0] * r_s2) + (coeffs[3] * g2)) * _theta * _t_inc;
                    const T s2_b = ((coeffs[1] * r_s2) + (coeffs[4] * g2) - _r_s2) * _theta * _t_inc;
                    const T s2_c = ((coeffs[2] * r_s2) + (coeffs[5] * g2)) * _theta * _t_inc;
        
                    const T s2_s = (s2_a * src_offset0[s1_idx_upper + i - 1]) + 
                                   (s2_b * src_offset0[s1_idx_upper + i    ]) + 
                                   (s2_c * src_offset0[s1_idx_upper + i + 1]);
                    
                    /* Update */
                    dst_offset[s1_idx_upper + i] = src_offset0[s1_idx_upper + i] + (s1_s + s2_s);
                }
                
                
                /* s1 = s1_min */
                const auto &lb = _g.get_value_coeffs_lower();
                if (lb.first)
                {
                    dst_offset[s1_idx_lower + i] = lb.second;
                }
                else
                {
                    /* S1 coefficients */
                    const T g1 = _half_sigma_s1_sq * s1_min * s1_min;
                    const T r_s1 = _r_s1 * s1_min;
                
                    _g.get_delta_coeffs_lower(&coeffs[0], 1, i);
                    const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                    const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                    const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
                    
                    const T s1_s = (s1_a * src_offset0[s1_idx_lower + i]) + 
                                   (s1_b * src_offset0[s1_idx_p1    + i]) + 
                                   (s1_c * src_offset0[s1_idx_p2    + i]);
                    
        
                    /* S2 coefficients */
                    const T g2 = _half_sigma_s2_sq * s2 * s2;
                    const T r_s2 = _r_s2 * s2;
                
                    _g.get_delta_coeffs_middle(&coeffs[0], 0, i);
                    _g.get_gamma_coeffs_middle(&coeffs[3], 0, i);
                    const T s2_a = ((coeffs[0] * r_s2) + (coeffs[3] * g2)) * _theta * _t_inc;
                    const T s2_b = ((coeffs[1] * r_s2) + (coeffs[4] * g2) - _r_s2) * _theta * _t_inc;
                    const T s2_c = ((coeffs[2] * r_s2) + (coeffs[5] * g2)) * _theta * _t_inc;
        
                    const T s2_s = (s2_a * src_offset0[s1_idx_lower + i - 1]) + 
                                   (s2_b * src_offset0[s1_idx_lower + i    ]) + 
                                   (s2_c * src_offset0[s1_idx_lower + i + 1]);
                    
                    /* Update */
                    dst_offset[s1_idx_lower + i] = src_offset0[s1_idx_lower + i] + (s1_s + s2_s);
                }
            }
        
            /* The corners */
            /* Apply the s1 = s1_min , s2 = s2_min boundary condition */
            {
                /* S1 coefficients */
                const T g1 = _half_sigma_s1_sq * s1_min * s1_min;
                const T r_s1 = _r_s1 * s1_min;
              
                _g.get_delta_coeffs_lower(&coeffs[0], 1);
                const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
              
                const T s1_s = (s1_a * src_offset0[0]) + 
                               (s1_b * src_offset0[1]) + 
                               (s1_c * src_offset0[2]);

                /* S2 coefficients */
                const T g2 = _half_sigma_s2_sq * s2_min * s2_min;
                const T r_s2 = _r_s2 * s2_min;
             
                _g.get_delta_coeffs_lower(&coeffs[0], 0);
                const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;
        
                const T s2_s = (s2_a * src_offset0[0]) + 
                               (s2_b * src_offset0[1]) + 
                               (s2_c * src_offset0[2]);
           
                /* Update */
                dst_offset[0] = src_offset0[0] + (s1_s + s2_s);
            }
        
            /* Apply the s1 = s1_min, s2 = s2_max boundary condition */
            {
                /* S1 coefficients */
                const T g1 = _half_sigma_s1_sq * s1_min * s1_min;
                const T r_s1 = _r_s1 * s1_min;
              
                _g.get_delta_coeffs_lower(&coeffs[0], 1);
                const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
              
                const T s1_s = (s1_a * src_offset0[(_g.size(0)  * 1) - 1]) + 
                               (s1_b * src_offset0[(_g.size(0)  * 2) - 1]) + 
                               (s1_c * src_offset0[(_g.size(0)  * 3) - 1]);

                /* S2 coefficients */
                const T g2 = _half_sigma_s2_sq * s2_max * s2_max;
                const T r_s2 = _r_s2 * s2_max;
            
                _g.get_delta_coeffs_upper(&coeffs[0], 0);
                const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;

                const T s2_s = (s2_a * src_offset0[_g.size(0) - 3]) + 
                               (s2_b * src_offset0[_g.size(0) - 2]) + 
                               (s2_c * src_offset0[_g.size(0) - 1]);
                    
                /* Update */
                dst_offset[_g.size(0) - 1] = src_offset0[_g.size(0) - 1] + (s1_s + s2_s);
            }
                
            /* Apply the s1 = s1_max, s2 = s2_min boundary condition */
            {
                /* S1 coefficients */
                const T g1 = _half_sigma_s1_sq * s1_max * s1_max;
                const T r_s1 = _r_s1 * s1_max;
               
                _g.get_delta_coeffs_upper(&coeffs[0], 1);
                const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
               
                const T s1_s = (s1_a * src_offset0[s1_idx_m2   ]) + 
                               (s1_b * src_offset0[s1_idx_m1   ]) + 
                               (s1_c * src_offset0[s1_idx_upper]);

                /* S2 coefficients */
                const T g2 = _half_sigma_s2_sq * s2_min * s2_min;
                const T r_s2 = _r_s2 * s2_min;
             
                _g.get_delta_coeffs_lower(&coeffs[0], 0);
                const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;
        
                const T s2_s = (s2_a * src_offset0[s1_idx_upper]) + 
                               (s2_b * src_offset0[s1_idx_upper + 1]) + 
                               (s2_c * src_offset0[s1_idx_upper + 2]);           
                /* Update */
                dst_offset[s1_idx_upper] = src_offset0[s1_idx_upper] + (s1_s + s2_s);
            }

            /* Apply the s1 = s1_max, s2 = s2_max boundary condition */
            {
                const int nsn       = _g.size(0) * _g.size(1);
                const int nsn_m1    = nsn - 1;
                
                /* S1 coefficients */
                const T g1 = _half_sigma_s1_sq * s1_max * s1_max;
                const T r_s1 = _r_s1 * s1_max;
               
                _g.get_delta_coeffs_upper(&coeffs[0], 1);
                const T s1_a = ((coeffs[0] * r_s1)) * _theta * _t_inc;
                const T s1_b = ((coeffs[1] * r_s1) - _r_s1) * _theta * _t_inc;
                const T s1_c = ((coeffs[2] * r_s1)) * _theta * _t_inc;
               
                const T s1_s = (s1_a * src_offset0[nsn_m1 - (_g.size(0) << 1)]) + 
                               (s1_b * src_offset0[nsn_m1 - _g.size(0)       ]) + 
                               (s1_c * src_offset0[nsn_m1                    ]);

                /* S2 coefficients */
                const T g2 = _half_sigma_s2_sq * s2_max * s2_max;
                const T r_s2 = _r_s2 * s2_max;
              
                _g.get_delta_coeffs_upper(&coeffs[0], 0);
                const T s2_a = ((coeffs[0] * r_s2)) * _theta * _t_inc;
                const T s2_b = ((coeffs[1] * r_s2) - _r_s2) * _theta * _t_inc;
                const T s2_c = ((coeffs[2] * r_s2)) * _theta * _t_inc;

                const T s2_s = (s2_a * src_offset0[nsn - 3]) + 
                               (s2_b * src_offset0[nsn - 2]) + 
                               (s2_c * src_offset0[nsn - 1]);
              
                /* Update */
                dst_offset[nsn - 1] = src_offset0[nsn - 1] + (s1_s + s2_s);
            }
        }

    private :
        const grid<T>   &   _g;
        const T             _theta;
        const T             _t_inc;
        const T             _sigma_s1;
        const T             _sigma_s2;
        const T             _half_sigma_s1_sq;
        const T             _half_sigma_s2_sq;
        const T             _r_s1;
        const T             _r_s2;
        const T             _corr; 
};


template <class T>
class implicit_corrector
{
    public :
        /* CTOR */
        implicit_corrector(const grid<T> &g, const T r, const T theta, const T t_inc, 
            const T sigma, const int size, const int step, const int dim)
            : _g(g), matrix_data_k(new T [_size * _size]), matrix_data_kp1(new T [_size * _size]), 
              corrected(new T [_size]), r(r), _theta(theta), _t_inc(t_inc), 
              half_sigma_sq(0.5 * sigma * sigma), _size(_size), _step(step), _dim(dim)
        { }

        /* DTOR */
        ~implicit_corrector()
        {
            delete [] corrected;
            delete [] matrix_data_k;
            delete [] matrix_data_kp1;
        }
        
        /* Implicit correction algorithm of the form F(A) = B + _theta * _t_inc * (F(A) - F(C)). 
           F is the 1D BS PDE. A, B and C may be any data sets to correct from */
        void operator()(const T *const src_offset0, const T *const src_offset1, T *const dst_offset)
        {
            T coeffs[6];
            
            /* TODO -- Perhaps not needed now the matrix class is fixed */
            memset(matrix_data_k, 0, (_size * _size * sizeof(T)));

            /* Populate the matrix */
            for (int i = 1; i < _size - 1; i++)
            {
                const T s = _g.s_grid(_dim)[i];

                const T g = half_sigma_sq * s * s;
                const T r_s = r * s;

                _g.get_delta_coeffs_middle(&coeffs[0], 2 - _dim, i);
                _g.get_gamma_coeffs_middle(&coeffs[3], 2 - _dim, i);
                const T s_a = ((coeffs[0] * r_s) + (coeffs[3] * g)) * _theta * _t_inc;
                const T s_b = ((coeffs[1] * r_s) + (coeffs[4] * g) - r) * _theta * _t_inc;
                const T s_c = ((coeffs[2] * r_s) + (coeffs[5] * g)) * _theta * _t_inc;
              
                matrix_data_k[(i * _size) + i - 1] = -s_a;
                matrix_data_k[(i * _size) + i    ] = 1.0 - s_b;
                matrix_data_k[(i * _size) + i + 1] = -s_c;

                matrix_data_kp1[i]  = src_offset1[_step * i]; 
                matrix_data_kp1[i] -= (s_a * src_offset0[(_step * (i - 1))]);
                matrix_data_kp1[i] -= (s_b * src_offset0[(_step * i)]);
                matrix_data_kp1[i] -= (s_c * src_offset0[(_step * (i + 1))]);
            }

            /* Boundary conditions */
            /* s = s_min */
            {        
                const auto &b = _g.get_value_coeffs_lower();
                if (b.first)
                {
                    matrix_data_k[0] = 1.0;
                    matrix_data_k[1] = 0.0;
                    matrix_data_k[2] = 0.0;
        
                    matrix_data_kp1[0] = b.second;
                }
                else
                {
                    _g.get_delta_coeffs_lower(&coeffs[0], _dim);
    
                    const T s = _g.s_grid(_dim)[0];
                    const T r_s = r * s;
            
                    const T s_a = ((coeffs[0] * r_s)) * _theta * _t_inc;
                    const T s_b = ((coeffs[1] * r_s) - r) * _theta * _t_inc;
                    const T s_c = ((coeffs[2] * r_s)) * _theta * _t_inc;
        
                    matrix_data_k[0] = (1.0 - s_a);
                    matrix_data_k[1] = -s_b;
                    matrix_data_k[2] = -s_c;
        
                    matrix_data_kp1[0]  = src_offset1[0]; 
                    matrix_data_kp1[0] -= (1.0 + s_a) * src_offset0[0];
                    matrix_data_kp1[0] -= s_b         * src_offset0[1];
                    matrix_data_kp1[0] -= s_c         * src_offset0[2];
                }
            }

            /* s = s_max*/
            {
                const auto &b = _g.get_value_coeffs_upper();
                const int matrix_offset = ((_size - 1) * _size) + (_size - 1);
                if (b.first)
                {
                    matrix_data_k[matrix_offset - 2] = 0.0;
                    matrix_data_k[matrix_offset - 1] = 0.0;
                    matrix_data_k[matrix_offset    ] = 1.0;
        
                    matrix_data_kp1[_size - 1] = b.second;
                }
                else
                {
                    _g.get_delta_coeffs_upper(&coeffs[0], _dim);
            
                    const T s = _g.s_grid(_dim)[_size - 1];
                    const T r_s = r * s;
            
                    const T s_a = ((coeffs[0] * r_s)) * _theta * _t_inc;
                    const T s_b = ((coeffs[1] * r_s) - r) * _theta * _t_inc;
                    const T s_c = ((coeffs[2] * r_s)) * _theta * _t_inc;
        
                    matrix_data_k[matrix_offset - 2] = -s_a;
                    matrix_data_k[matrix_offset - 1] = -s_b;
                    matrix_data_k[matrix_offset    ] = (1.0 - s_c);
        
                    matrix_data_kp1[_size - 1]  = src_offset1[(_step * (_size - 1))]; 
                    matrix_data_kp1[_size - 1] -= s_a         * src_offset0[(_step * (_size - 3))]; 
                    matrix_data_kp1[_size - 1] -= s_b         * src_offset0[(_step * (_size - 2))]; 
                    matrix_data_kp1[_size - 1] -= (1.0 + s_c) * src_offset0[(_step * (_size - 1))]; 
                }
            }
        
            /* Solve matrix */
            matrix<T> matrix_kp1(&matrix_data_kp1[0], _size, 1, true, true);
            matrix<T> matrix_k(&matrix_data_k[0], _size, _size, true, true);
            matrix_k.gauss_solve(matrix_kp1, &corrected[0], Tridiagonal);
            for (int i = 0; i < _size; i++)
            {
                dst_offset[(i * _step)] = corrected[i];
            }
        }

    private :
        const grid<T>   &   _g;
        T *const            matrix_data_k;
        T *const            matrix_data_kp1;
        T *const            corrected;
        const T             r;
        const T             _theta;
        const T             _t_inc;
        const T             half_sigma_sq;
        const int           _size;
        const int           _step;
        const int           _dim;
};


template<class T>
class explicit_corrector
{
    public :
        /* CTOR */
        explicit_corrector(const grid<T> &g, const T theta, const T t_inc,
            const T r_s1, const T r_s2, const T sigma_s1, const T sigma_s2, const T corr_s1_s2)
            : _g(g), _theta(theta), _t_inc(t_inc),
              _r_s1(r_s1), _r_s2(r_s2), _sigma_s1(sigma_s1), _sigma_s2(sigma_s2), 
              _corr_s1_s2(corr_s1_s2), _half_sigma_s1_sq(0.5 * sigma_s1 * sigma_s1), 
              _half_sigma_s2_sq(0.5 * sigma_s2 * sigma_s2)
        { }

        /* DTOR */
        ~explicit_corrector() { }
        
        /* Cross correction algorithm of the form F(A) = A + _theta * _t_inc * (F(B) - F(C)). 
           F is a cross term in a multi-asset BS PDE. A, B and C may be any data sets 
           to correct from */
        void operator()(const T *const src_offset0, const T *const src_offset1, T *const dst_offset, const int step)
        {
            T coeffs[6];
            
            const int s1_idx        = step * _g.size(1);
            const int s1_idx_p1     = s1_idx + _g.size(1);
            const int s1_idx_m1     = s1_idx - _g.size(1);
            const T s1              = _g.s_grid(1)[step];/* Todo - wrong dir */
            
            /* Fill in the column */
            for (int i = 1; i < _g.size(0) - 1; i++)
            {
                const T s2 = _g.s_grid(0)[i];/* Todo - wrong dir */
             
                /* S1 coefficients */
                const T g1 = _half_sigma_s1_sq * s1 * s1;
                const T r_s1 = _r_s1 * s1;
             
                _g.get_delta_coeffs_middle(&coeffs[0], 1, i);
                _g.get_gamma_coeffs_middle(&coeffs[3], 1, i);
                const T s1_a = ((coeffs[0] * r_s1) + (coeffs[3] * g1)) * _theta * _t_inc;
                const T s1_b = ((coeffs[1] * r_s1) + (coeffs[4] * g1) - _r_s1) * _theta * _t_inc;
                const T s1_c = ((coeffs[2] * r_s1) + (coeffs[5] * g1)) * _theta * _t_inc;
             
                const T s1_s = (s1_a * (src_offset0[s1_idx_p1 + i] - src_offset1[s1_idx_p1 + i])) + 
                               (s1_b * (src_offset0[s1_idx    + i] - src_offset1[s1_idx    + i])) + 
                               (s1_c * (src_offset0[s1_idx_m1 + i] - src_offset1[s1_idx_m1 + i]));
             
        
                /* S2 coefficients */
                const T g2 = _half_sigma_s2_sq * s2 * s2;
                const T r_s2 = _r_s2 * s2;
             
                _g.get_delta_coeffs_middle(&coeffs[0], 0, i);
                _g.get_gamma_coeffs_middle(&coeffs[3], 0, i);
                const T s2_a = ((coeffs[0] * r_s2) + (coeffs[3] * g2)) * _theta * _t_inc;
                const T s2_b = ((coeffs[1] * r_s2) + (coeffs[4] * g2) - _r_s2) * _theta * _t_inc;
                const T s2_c = ((coeffs[2] * r_s2) + (coeffs[5] * g2)) * _theta * _t_inc;
        
                const T s2_s = (s2_a * (src_offset0[s1_idx + i + 1] - src_offset1[s1_idx + i + 1])) + 
                               (s2_b * (src_offset0[s1_idx + i    ] - src_offset1[s1_idx + i    ])) + 
                               (s2_c * (src_offset0[s1_idx + i - 1] - src_offset1[s1_idx + i - 1]));
                                   

                /* Update */
                dst_offset[s1_idx + i] = (s1_s + s2_s) + dst_offset[s1_idx + i];
            }
    
            /* Apply the s2 = 0 boundary condition */
            {
                /* S1 greeks */
                const T delta_s1    = ((src_offset0[s1_idx_p1] - src_offset0[s1_idx_m1]) -
                                       (src_offset1[s1_idx_p1] - src_offset1[s1_idx_m1])) * _s1_inc2_inv;
                const T rv_s1       = (src_offset0[s1_idx] - src_offset1[s1_idx]) * _r_s1;

                /* s2 greeks */
                const T s2          = 0.0;
                const T delta_s2    = ((src_offset0[s1_idx + 1] - src_offset0[s1_idx]) -
                                       (src_offset1[s1_idx + 1] - src_offset1[s1_idx])) * s2_inc_inv;
                const T rv_s2       = (src_offset0[s1_idx] - src_offset1[s1_idx]) * _r_s2;
                
                /* Cross partial */
                const T delta_m1    = ((src_offset0[s1_idx_m1 + 1] - src_offset0[s1_idx_m1]) -
                                       (src_offset1[s1_idx_m1 + 1] - src_offset1[s1_idx_m1])) * s2_inc_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + 1] - src_offset0[s1_idx_p1]) -
                                       (src_offset1[s1_idx_p1 + 1] - src_offset1[s1_idx_p1])) * s2_inc_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T corre       = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                      (_r_s2 * s2 * delta_s2) - rv_s2 +
                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * _corr_s1_s2 * x_gamma);
                dst_offset[s1_idx]   = (corre * _t_inc * _theta) + dst_offset[s1_idx];
            }
        
            /* Apply the s2 = _s2_max boundary condition */
            {
                /* S1 greeks */
                const T delta_s1    = ((src_offset0[s1_idx_p1 + _g.size(0)] - src_offset0[s1_idx_m1 + _g.size(0)]) -
                                       (src_offset1[s1_idx_p1 + _g.size(0)] - src_offset1[s1_idx_m1 + _g.size(0)])) * _s1_inc2_inv;
                const T rv_s1       = (src_offset0[s1_idx + _g.size(0)] - src_offset1[s1_idx + _g.size(0)]) * _r_s1;

                /* s2 greeks */
                const T s2          = _g.s_grid(0)[_g.size(0)];
                const T delta_s2    = ((src_offset0[s1_idx + _g.size(0)] - src_offset0[s1_idx + _g.size(0) - 1]) -
                                       (src_offset1[s1_idx + _g.size(0)] - src_offset1[s1_idx + _g.size(0) - 1])) * s2_inc_inv;
                const T rv_s2       = (src_offset0[s1_idx + _g.size(0)] - src_offset1[s1_idx + _g.size(0)]) * _r_s2;
                
                /* Cross partial */
                const T delta_m1    = ((src_offset0[s1_idx_m1 + _g.size(0)] - src_offset0[s1_idx_m1 + _g.size(0) - 1]) -
                                       (src_offset1[s1_idx_m1 + _g.size(0)] - src_offset1[s1_idx_m1 + _g.size(0) - 1])) * s2_inc_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + _g.size(0)] - src_offset0[s1_idx_p1 + _g.size(0) - 1]) -
                                       (src_offset1[s1_idx_p1 + _g.size(0)] - src_offset1[s1_idx_p1 + _g.size(0) - 1])) * s2_inc_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T corre       = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                      (_r_s2 * s2 * delta_s2) - rv_s2 +
                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * _corr_s1_s2 * x_gamma);
                dst_offset[s1_idx_p1 - 1]   = (corre * _t_inc * _theta) + dst_offset[s1_idx_p1 - 1];
            }
        }

    private :
        const grid<T>   &   _g;
        const T             _theta;
        const T             _t_inc;
        const T             _r_s1;
        const T             _r_s2;
        const T             _sigma_s1;
        const T             _sigma_s2;
        const T             _corr_s1_s2;
        const T             _half_sigma_s1_sq;
        const T             _half_sigma_s2_sq;
};


template<class T>
class cross_term_corrector
{
    public :
        /* CTOR */
        cross_term_corrector(const grid<T> &g, const T theta, const T t_inc,
            const T sigma_s1, const T sigma_s2, const T corr_s1_s2)
            : _g(g), _theta(theta), _t_inc(t_inc), 
              _sigma_s1_s2(sigma_s1 * sigma_s1), _corr_s1_s2(corr_s1_s2)
        { }

        /* DTOR */
        ~cross_term_corrector() { }
        
        /* Cross correction algorithm of the form F(A) = A + _theta * _t_inc * (F(B) - F(C)). 
           F is a cross term in a multi-asset BS PDE. A, B and C may be any data sets 
           to correct from */
        void operator()(const T *const src_offset0, const T *const src_offset1, T *const dst_offset, const int step)
        {
            const int s1_idx    = step * _g.size(1);
            const int s1_idx_p1 = s1_idx + _g.size(1);
            const int s1_idx_m1 = s1_idx - _g.size(1);
            const T s1          = _g.s_grid(1)[step];
                
            /* Fill in the column */
            for (int i = 1; i < _g.size(0) - 1; i++)
            {
                /* Cross partial */
                const T s2          = _g.s_grid(0)[i];
                const T delta_m1    = ((src_offset0[s1_idx_m1 + i + 1] - src_offset0[s1_idx_m1 + i - 1]) -
                                       (src_offset1[s1_idx_m1 + i + 1] - src_offset1[s1_idx_m1 + i - 1])) * _s2_inc2_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + i + 1] - src_offset0[s1_idx_p1 + i - 1]) -
                                       (src_offset1[s1_idx_p1 + i + 1] - src_offset1[s1_idx_p1 + i - 1])) * _s2_inc2_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T cross           = _sigma_s1_s2 * s1 * s2 * _corr_s1_s2 * x_gamma;
                dst_offset[s1_idx + i]  = (cross * _theta * _t_inc) + dst_offset[s1_idx + i];
            }
    
            /* Apply the s2 = 0 boundary condition */
            //dst_offset[s1_idx] = (2.0 * dst_offset[s1_idx + 1]) - dst_offset[s1_idx + 2];
//            {
//                /* Cross partial */
//                const T s2          = 0.0;
//                const T delta_m1    = ((src_offset0[s1_idx_m1 + i + 1] - src_offset0[s1_idx_m1 + i - 1]) -
//                                       (src_offset1[s1_idx_m1 + i + 1] - src_offset1[s1_idx_m1 + i - 1])) * s2_inc_inv;
//                const T delta_p1    = ((src_offset0[s1_idx_p1 + i + 1] - src_offset0[s1_idx_p1 + i - 1]) -
//                                       (src_offset1[s1_idx_p1 + i + 1] - src_offset1[s1_idx_p1 + i - 1])) * s2_inc_inv;
//                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;
//
//                /* Update */
//                const T cross           = _sigma_s1_s2 * s1 * s2 * _corr_s1_s2 * x_gamma;
//                dst_offset[s1_idx + i]  = (cross * _theta * _t_inc) + dst_offset[s1_idx + i];
//            }
        
            /* Apply the s2 = _s2_max boundary condition */
            //dst_offset[s1_idx_p1 - 1] = (2.0 * dst_offset[s1_idx_p1 - 2]) - dst_offset[s1_idx_p1 - 3];
        }

    private :
        const grid<T>   &   _g;
        const T             _theta;
        const T             _t_inc;
        const T             _sigma_s1_s2;
        const T             _corr_s1_s2;
};

#endif
