#ifndef __CORRECTORS_H__
#define __CORRECTORS_H__


#include "matrix.h"


/* Explicit predictor step */
template <class T>
class explicit_predictor
{
    public :
        /* CTOR */
        explicit_predictor(const T theta, const T t_inc, const T s1_inc, const T s2_inc,
            const T s1_max, const T s2_max,
            const T sigma_s1, const T sigma_s2, const T r_s1, const T r_s2, const T corr,
            const int nas1, const int nas2)
            : _theta(theta), _t_inc(t_inc), _s1_inc(s1_inc), _s2_inc(s2_inc),
              _s1_inc2_inv(0.5 / s1_inc), _s2_inc2_inv(0.5 / s2_inc),
              _s1_inc_sq_inv(1.0 / (s1_inc * s1_inc)),
              _s2_inc_sq_inv(1.0 / (s2_inc * s2_inc)),
              _sigma_s1(sigma_s1), _sigma_s2(sigma_s2), 
              _half_sigma_s1_sq(0.5 * sigma_s1 * sigma_s1), 
              _half_sigma_s2_sq(0.5 * sigma_s2 * sigma_s2), 
              _r_s1(r_s1), _r_s2(r_s2), _corr(corr),
              _s1_max(s1_max), _s2_max(s2_max),
              _nas1(nas1), _nas2(nas2), _nsn(nas1 * nas2)
        { }

        /* DTOR */
        ~explicit_predictor() {  };
        
        /* Implicit correction algorithm of the form F(A) = B + _theta * _t_inc * (F(A) - F(C)). 
           F is the 1D BS PDE. A, B and C may be any data sets to correct from */
        void operator()(const T *const src_offset0, T *const dst_offset)
        {
            /* Fill in the row */
            for (int i = 1; i < this->_nas1 - 1; i++)
            {
                const int s1_idx        = i * this->_nas2;
                const int s1_idx_p1     = s1_idx + this->_nas2;
                const int s1_idx_m1     = s1_idx - this->_nas2;
                const T s1              = static_cast<T>(i) * _s1_inc;
    
                /* Fill in the column */
                for (int j = 1; j < this->_nas2 - 1; j++)
                {
                    /* S1 greeks */
                    const T delta_s1   = (src_offset0[s1_idx_p1 + j] - src_offset0[s1_idx_m1 + j]) * _s1_inc2_inv;
                    const T gamma_s1   = (src_offset0[s1_idx_p1 + j] - (2.0 * src_offset0[s1_idx + j]) + src_offset0[s1_idx_m1 + j]) * _s1_inc_sq_inv;
                    const T rv_s1      = _r_s1 * src_offset0[s1_idx + j];
    
                    /* s2 greeks */
                    const T s2         = static_cast<T>(j) * _s2_inc;
                    const T delta_s2   = (src_offset0[s1_idx + j + 1] - src_offset0[s1_idx + j - 1]) * _s2_inc2_inv;
                    const T gamma_s2   = (src_offset0[s1_idx + j + 1] - (2.0 * src_offset0[s1_idx + j]) + src_offset0[s1_idx + j - 1]) * _s2_inc_sq_inv;
                    const T rv_s2      = _r_s2 * src_offset0[s1_idx + j];
                    
                    /* Cross partial */
                    const T delta_m1   = (src_offset0[s1_idx_m1 + j + 1] - src_offset0[s1_idx_m1 + j - 1]) * _s2_inc2_inv;
                    const T delta_p1   = (src_offset0[s1_idx_p1 + j + 1] - src_offset0[s1_idx_p1 + j - 1]) * _s2_inc2_inv;
                    const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;
    
                    /* Update */
                    const T _theta     = (_half_sigma_s1_sq * s1 * s1 * gamma_s1) + (_r_s1 * s1 * delta_s1) - rv_s1 +
                                         (_half_sigma_s2_sq * s2 * s2 * gamma_s2) + (_r_s2 * s2 * delta_s2) - rv_s2 +
                                         (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                    dst_offset[s1_idx + j]   = (_theta * _t_inc) + src_offset0[s1_idx + j];
                }
        
                /* Apply the s2 = 0 boundary condition */
                //dst_offset[s1_idx] = (2.0 * dst_offset[s1_idx + 1]) - dst_offset[s1_idx + 2];
                {
                    /* S1 greeks */
                    const T delta_s1   = (src_offset0[s1_idx_p1] - src_offset0[s1_idx_m1]) * _s1_inc2_inv;
                    const T gamma_s1   = (src_offset0[s1_idx_p1] - (2.0 * src_offset0[s1_idx]) + src_offset0[s1_idx_m1]) * _s1_inc_sq_inv;
                    const T rv_s1      = _r_s1 * src_offset0[s1_idx];
    
                    /* s2 greeks */
                    const T s2         = 0.0;
                    const T delta_s2   = (src_offset0[s1_idx + 2] - src_offset0[s1_idx]) * _s2_inc2_inv;
                    const T rv_s2      = _r_s2 * src_offset0[s1_idx];
                
                    /* Cross partial */
//                  const T delta_m1   = (src_offset0[s1_idx_m1 + 2] - src_offset0[s1_idx_m1]) * _s2_inc2_inv;
//                  const T delta_p1   = (src_offset0[s1_idx_p1 + 2] - src_offset0[s1_idx_p1]) * _s2_inc2_inv;
//                  const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;
    
                    /* Update */
                    const T _theta      = (_half_sigma_s1_sq * s1 * s1 * gamma_s1) + (_r_s1 * s1 * delta_s1) - rv_s1 +
                                         (_r_s2 * s2 * delta_s2) - rv_s2;// +
//                                              (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                    dst_offset[s1_idx]   = (_theta * _t_inc) + src_offset0[s1_idx];
                }
        
                /* Apply the s2 = _s2_max boundary condition */
                //dst_offset[s1_idx_p1 - 1] = (2.0 * dst_offset[s1_idx_p1 - 2]) - dst_offset[s1_idx_p1 - 3];
                {
                    /* S1 greeks */
                    const T delta_s1   = (src_offset0[s1_idx_p1 + this->_nas2 - 1] - src_offset0[s1_idx_m1 + this->_nas2 - 1]) * _s1_inc2_inv;
                    const T gamma_s1   = (src_offset0[s1_idx_p1 + this->_nas2 - 1] - (2.0 * src_offset0[s1_idx + this->_nas2 - 1]) + src_offset0[s1_idx_m1 + this->_nas2 - 1]) * _s1_inc_sq_inv;
                    const T rv_s1      = _r_s1 * src_offset0[s1_idx + this->_nas2 - 1];
    
                    /* s2 greeks */
                    const T s2         = _s2_max;
                    const T delta_s2   = (src_offset0[s1_idx + this->_nas2 - 1] - src_offset0[s1_idx + this->_nas2 - 3]) * _s2_inc2_inv;
                    const T rv_s2      = _r_s2 * src_offset0[s1_idx + this->_nas2 - 1];
                    
                    /* Cross partial */
//                  const T delta_m1   = (src_offset0[s1_idx_m1 + this->_nas2 - 1] - src_offset0[s1_idx_m1 + this->_nas2 - 3]) * _s2_inc2_inv;
//                  const T delta_p1   = (src_offset0[s1_idx_p1 + this->_nas2 - 1] - src_offset0[s1_idx_p1 + this->_nas2 - 3]) * _s2_inc2_inv;
//                  const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;

                    /* Update */
                    const T _theta      = (_half_sigma_s1_sq * s1 * s1 * gamma_s1) + (_r_s1 * s1 * delta_s1) - rv_s1 +
                                         (_r_s2 * s2 * delta_s2) - rv_s2;// +
//                                          (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                    dst_offset[s1_idx + this->_nas2 - 1]   = (_theta * _t_inc) + src_offset0[s1_idx + this->_nas2 - 1];
                }
            }
    
            /* Upper and lower boundary condition in s1 */
            const int s1_idx_upper  = (this->_nas1 - 1) * this->_nas2;
            const int s1_idx_m1     = s1_idx_upper - this->_nas2;
            const int s1_idx_m2     = s1_idx_m1 - this->_nas2;

            const int s1_idx_lower  = 0;
//            const int s1_idx_p1     = this->_nas2;
            const int s1_idx_p2     = (this->_nas2 << 1);
            for (int i = 1; i < this->_nas2 - 1; i++)
            {
                //dst_offset[s1_idx_upper + i] = (2.0 * dst_offset[s1_idx_m1 + i]) - dst_offset[s1_idx_m2 + i];
                {
                    /* S1 greeks */
                    const T s1         = _s1_max;
                    const T delta_s1   = (src_offset0[s1_idx_upper + i] - src_offset0[s1_idx_m2 + i]) * _s1_inc2_inv;
                    const T rv_s1      = _r_s1 * src_offset0[s1_idx_upper + i];
    
                    /* s2 greeks */
                    const T s2         = static_cast<T>(i) * _s2_inc;
                    const T delta_s2   = (src_offset0[s1_idx_upper + i + 1] - src_offset0[s1_idx_upper + i - 1]) * _s2_inc2_inv;
                    const T gamma_s2   = (src_offset0[s1_idx_upper + i + 1] - (2.0 * src_offset0[s1_idx_upper + i]) + src_offset0[s1_idx_upper + i - 1]) * _s2_inc_sq_inv;
                    const T rv_s2      = _r_s2 * src_offset0[s1_idx_upper + i];
                    
                    /* Cross partial */
//                    const T delta_m1   = (src_offset0[s1_idx_m2 + i + 1] - src_offset0[s1_idx_m2 + i - 1]) * _s2_inc2_inv;
//                    const T delta_p1   = (src_offset0[s1_idx_upper + i + 1] - src_offset0[s1_idx_upper + i - 1]) * _s2_inc2_inv;
//                    const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;

                    /* Update */
                    const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                         (_half_sigma_s2_sq * s2 * s2 * gamma_s2) + (_r_s2 * s2 * delta_s2) - rv_s2;// +
//                                          (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                    dst_offset[s1_idx_upper + i]   = (_theta * _t_inc) + src_offset0[s1_idx_upper + i];
                }
            
                //dst_offset[s1_idx_lower + i] = (2.0 * dst_offset[s1_idx_p1 + i]) - dst_offset[s1_idx_p2 + i];
                {
                    /* S1 greeks */
                    const T s1         = 0.0;
                    const T delta_s1   = (src_offset0[s1_idx_p2 + i] - src_offset0[s1_idx_lower + i]) * _s1_inc2_inv;
                    const T rv_s1      = _r_s1 * src_offset0[s1_idx_lower + i];

                    /* s2 greeks */
                    const T s2         = static_cast<T>(i) * _s2_inc;
                    const T delta_s2   = (src_offset0[s1_idx_lower + i + 1] - src_offset0[s1_idx_lower + i - 1]) * _s2_inc2_inv;
                    const T gamma_s2   = (src_offset0[s1_idx_lower + i + 1] - (2.0 * src_offset0[s1_idx_lower + i]) + src_offset0[s1_idx_lower + i - 1]) * _s2_inc_sq_inv;
                    const T rv_s2      = _r_s2 * src_offset0[s1_idx_lower + i];
                
                    /* Cross partial */
//                    const T delta_m1   = (src_offset0[s1_idx_lower + i + 1] - src_offset0[s1_idx_lower + i - 1]) * _s2_inc2_inv;
//                    const T delta_p1   = (src_offset0[s1_idx_p2 + i + 1] - src_offset0[s1_idx_p2 + i - 1]) * _s2_inc2_inv;
//                    const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;

                    /* Update */
                    const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                         (_half_sigma_s2_sq * s2 * s2 * gamma_s2) + (_r_s2 * s2 * delta_s2) - rv_s2;// +
//                                          (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                    dst_offset[s1_idx_lower + i]   = (_theta * _t_inc) + src_offset0[s1_idx_lower + i];
                }
            }
        
            /* The corners */
            /* Apply the s1 = 0 , s2 = 0 boundary condition */
//            fd_grid[t_offset] = (2.0 * dst_offset[1]) - dst_offset[2];
            {
                /* S1 greeks */
                const T s1         = 0.0;
                const T delta_s1   = (src_offset0[(this->_nas2 << 1)] - src_offset0[0]) * _s1_inc2_inv;
                const T rv_s1      = _r_s1 * src_offset0[0];

                /* s2 greeks */
                const T s2         = 0.0;
                const T delta_s2   = (src_offset0[2] - src_offset0[0]) * _s2_inc2_inv;
                const T rv_s2      = _r_s2 * src_offset0[0];
           
                /* Update */
                const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                     (_r_s2 * s2 * delta_s2) - rv_s2;
                dst_offset[0] = (_theta * _t_inc) + src_offset0[0];
            }
        
            /* Apply the s1 = 0, s2 = _s2_max boundary condition */
            //dst_offset[this->_nas2 - 1] = (2.0 * dst_offset[this->_nas2 - 2]) - dst_offset[this->_nas2 - 3];
            {
                /* S1 greeks */
                const T s1         = 0.0;
                const T delta_s1   = (src_offset0[(this->_nas2  * 3) - 1] - src_offset0[this->_nas2 - 1]) * _s1_inc2_inv;
                const T rv_s1      = _r_s1 * src_offset0[this->_nas2 - 1];

                /* s2 greeks */
                const T s2         = _s2_max;
                const T delta_s2   = (src_offset0[this->_nas2 - 1] - src_offset0[this->_nas2 - 3]) * _s2_inc2_inv;
                const T rv_s2      = _r_s2 * src_offset0[this->_nas2 - 1];
           
                /* Update */
                const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                     (_r_s2 * s2 * delta_s2) - rv_s2;
                dst_offset[this->_nas2 - 1] = (_theta * _t_inc) + src_offset0[this->_nas2 - 1];
            }
                
            /* Apply the s1 = _s1_max, s2 = 0 boundary condition */
            //dst_offset[s1_idx_upper] = (2.0 * dst_offset[s1_idx_upper + 1]) - dst_offset[s1_idx_upper + 2];
            {
                /* S1 greeks */
                const T s1         = _s1_max;
                const T delta_s1   = (src_offset0[s1_idx_upper] - src_offset0[s1_idx_m2]) * _s1_inc2_inv;
                const T rv_s1      = _r_s1 * src_offset0[s1_idx_upper];

                /* s2 greeks */
                const T s2         = 0.0;
                const T delta_s2   = (src_offset0[s1_idx_upper + 2] - src_offset0[s1_idx_upper]) * _s2_inc2_inv;
                const T rv_s2      = _r_s2 * src_offset0[s1_idx_upper];
           
                /* Update */
                const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                     (_r_s2 * s2 * delta_s2) - rv_s2;
                dst_offset[s1_idx_upper] = (_theta * _t_inc) + src_offset0[s1_idx_upper];
            }

            /* Apply the s1 = _s1_max, s2 = _s2_max boundary condition */
            //dst_offset[_nsn - 1] = (2.0 * dst_offset[_nsn - 2]) - dst_offset[_nsn - 3];
            {
                /* S1 greeks */
                const T s1         = _s1_max;
                const T delta_s1   = (src_offset0[_nsn - 1] - src_offset0[_nsn - 1 - (this->_nas2 << 1)]) * _s1_inc2_inv;
                const T rv_s1      = _r_s1 * src_offset0[_nsn - 1];

                /* s2 greeks */
                const T s2         = _s2_max;
                const T delta_s2   = (src_offset0[_nsn - 1] - src_offset0[_nsn - 3]) * _s2_inc2_inv;
                const T rv_s2      = _r_s2 * src_offset0[_nsn - 1];
            
                /* Cross partial */
//                const T delta_m1   = (src_offset0[_nsn - 1 - (this->_nas2 << 1)] - src_offset0[_nsn - 3 - (this->_nas2 << 1)]) * _s2_inc2_inv;
//                const T delta_p1   = (src_offset0[_nsn - 1              ] - src_offset0[_nsn - 3              ]) * _s2_inc2_inv;
//                const T x_gamma    = (delta_p1 - delta_m1) * _s1_inc2_inv;
           
                /* Update */
                const T _theta      = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                     (_r_s2 * s2 * delta_s2) - rv_s2;// +
//                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * this->_corr * x_gamma);
                dst_offset[_nsn - 1] = (_theta * _t_inc) + src_offset0[_nsn - 1];
            }
        }

    private :
        const T     _theta;
        const T     _t_inc;
        const T     _s1_inc;
        const T     _s2_inc;
        const T     _s1_inc2_inv;
        const T     _s2_inc2_inv;
        const T     _s1_inc_sq_inv;
        const T     _s2_inc_sq_inv;
        const T     _sigma_s1;
        const T     _sigma_s2;
        const T     _half_sigma_s1_sq;
        const T     _half_sigma_s2_sq;
        const T     _r_s1;
        const T     _r_s2;
        const T     _corr;
        const T     _s1_max;
        const T     _s2_max;
        const int   _nas1; 
        const int   _nas2; 
        const int   _nsn;
};


template <class T>
class implicit_corrector
{
    public :
        /* CTOR */
        implicit_corrector(const T r, const T theta, const T s_inc, const T t_inc, 
            const T sigma, const int step, const int size)
            : matrix_data_k(new T [size * size]), matrix_data_kp1(new T [size * size]), 
              corrected(new T [size]), r(r), _theta(theta), s_inc(s_inc), 
              _t_inc(t_inc), half_sigma_sq(0.5 * sigma * sigma), t_inc_s_inc_sq(t_inc / (s_inc * s_inc)), 
              t_inc_s_inc(t_inc / s_inc), step(step), size(size)
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
            memset(matrix_data_k, 0, (size * size * sizeof(T)));

            /* Populate the matrix */
            for (int i = 1; i < size - 1; i++)
            {
                const T s      = static_cast<T>(i) * s_inc;
                const T s_sq   = s * s;

                const T g = half_sigma_sq * s_sq;
                const T r_s = r * s;

                const T a = (t_inc_s_inc_sq * g) - (0.5 * t_inc_s_inc * r_s);
                const T b = (_t_inc * -r * 0.5)   - (2.0 * t_inc_s_inc_sq * g);
                const T c = (t_inc_s_inc_sq * g) + (0.5 * t_inc_s_inc * r_s);

                matrix_data_k[(i * size) + i - 1] = -a * _theta;
                matrix_data_k[(i * size) + i    ] = 1.0 - b * _theta;
                matrix_data_k[(i * size) + i + 1] = -c * _theta;

                matrix_data_kp1[i]  = src_offset1[step * i]; 
                matrix_data_kp1[i] -= (a * _theta * src_offset0[(step * (i - 1))]);
                matrix_data_kp1[i] -= (b * _theta * src_offset0[(step * i)]);
                matrix_data_kp1[i] -= (c * _theta * src_offset0[(step * (i + 1))]);
            }

            /* Boundary conditions */
            /* s = 0.0 */
            {
                const T r_s = 0.0;
                const T b   = (_t_inc * -r * 0.5) - (t_inc_s_inc * r_s);
                const T c   = (t_inc_s_inc * r_s);
                matrix_data_kp1[0]  = src_offset1[0];
                matrix_data_kp1[0] -= (b * _theta * src_offset0[0]);
                matrix_data_kp1[0] -= (c * _theta * src_offset0[1]);
                
                matrix_data_k[0] = 1.0 - b * _theta;
                matrix_data_k[1] = -c * _theta;
            }
            
            /* s = s_max*/
            {
                const T s   = static_cast<T>(size) * s_inc;
                const T r_s = s * r;
                const T a   = -(t_inc_s_inc * r_s);
                const T b   =  (_t_inc * -r * 0.5) + (t_inc_s_inc * r_s);
                matrix_data_kp1[size - 1]  = src_offset1[(size - 1) * step];
                matrix_data_kp1[size - 1] -= (a * _theta * src_offset0[step * (size - 2)]);
                matrix_data_kp1[size - 1] -= (b * _theta * src_offset0[step * (size - 1)]);

                matrix_data_k[(size * size) - 2] = -a * _theta;
                matrix_data_k[(size * size) - 1] = 1.0 - b * _theta;
            }
        
            /* Solve matrix */
            matrix<T> matrix_kp1(&matrix_data_kp1[0], size, 1, true, true);
            matrix<T> matrix_k(&matrix_data_k[0], size, size, true, true);
            matrix_k.gauss_solve(matrix_kp1, &corrected[0], Tridiagonal);
            for (int i = 0; i < size; i++)
            {
                dst_offset[(i * step)] = corrected[i];
            }
        }

    private :
        T *const    matrix_data_k;
        T *const    matrix_data_kp1;
        T *const    corrected;
        const T     r;
        const T     _theta;
        const T     s_inc;
        const T     _t_inc;
        const T     half_sigma_sq;
        const T     t_inc_s_inc_sq;
        const T     t_inc_s_inc;
        const int   step;
        const int   size;
};


template<class T>
class explicit_corrector
{
    public :
        /* CTOR */
        explicit_corrector(const T theta, const T t_inc, const T s1_inc, const T s2_inc,
            const T r_s1, const T r_s2, const T sigma_s1, const T sigma_s2, const T corr_s1_s2, 
            const int size)
            : _theta(theta), _t_inc(t_inc), _s1_inc(s1_inc), _s2_inc(s2_inc), 
              _r_s1(r_s1), _r_s2(r_s2), _sigma_s1(sigma_s1), _sigma_s2(sigma_s2), 
              corr_s1_s2(corr_s1_s2), s2_inc_inv(1.0 / s2_inc), _s1_inc2_inv(0.5 / s1_inc), _s2_inc2_inv(0.5 / s2_inc), 
              _s1_inc_sq_inv(1.0 / (s1_inc * s1_inc)), _s2_inc_sq_inv(1.0 / (s2_inc * s2_inc)), 
              _half_sigma_s1_sq(0.5 * sigma_s1 * sigma_s1), _half_sigma_s2_sq(0.5 * sigma_s2 * sigma_s2), 
              size(size)
        { }

        /* DTOR */
        ~explicit_corrector() { }
        
        /* Cross correction algorithm of the form F(A) = A + _theta * _t_inc * (F(B) - F(C)). 
           F is a cross term in a multi-asset BS PDE. A, B and C may be any data sets 
           to correct from */
        void operator()(const T *const src_offset0, const T *const src_offset1, T *const dst_offset, const int step)
        {
            const int s1_idx        = step * size;
            const int s1_idx_p1     = s1_idx + size;
            const int s1_idx_m1     = s1_idx - size;
            const T s1              = static_cast<T>(step) * _s1_inc;
            
            /* Fill in the column */
            for (int i = 1; i < size - 1; i++)
            {
                /* S1 greeks */
                const T delta_s1    = ((src_offset0[s1_idx_p1 + i] - src_offset0[s1_idx_m1 + i]) -
                                       (src_offset1[s1_idx_p1 + i] - src_offset1[s1_idx_m1 + i])) * _s1_inc2_inv;
                const T gamma_s1    = ((src_offset0[s1_idx_p1 + i] - (2.0 * src_offset0[s1_idx + i]) + src_offset0[s1_idx_m1 + i]) -
                                       (src_offset1[s1_idx_p1 + i] - (2.0 * src_offset1[s1_idx + i]) + src_offset1[s1_idx_m1 + i])) * _s1_inc_sq_inv;
                const T rv_s1       = (src_offset0[s1_idx + i] - src_offset1[s1_idx + i]) * _r_s1;

                /* s2 greeks */
                const T s2          = static_cast<T>(i) * _s2_inc;
                const T delta_s2    = ((src_offset0[s1_idx + i + 1] - src_offset0[s1_idx + i - 1]) -
                                       (src_offset1[s1_idx + i + 1] - src_offset1[s1_idx + i - 1])) * _s2_inc2_inv;
                const T gamma_s2    = ((src_offset0[s1_idx + i + 1] - (2.0 * src_offset0[s1_idx + i]) + src_offset0[s1_idx + i - 1]) -
                                       (src_offset1[s1_idx + i + 1] - (2.0 * src_offset1[s1_idx + i]) + src_offset1[s1_idx + i - 1])) * _s2_inc_sq_inv;
                const T rv_s2       = (src_offset0[s1_idx + i] - src_offset1[s1_idx + i]) * _r_s2;
                
                /* Cross partial */
                const T delta_m1    = ((src_offset0[s1_idx_m1 + i + 1] - src_offset0[s1_idx_m1 + i - 1]) -
                                       (src_offset1[s1_idx_m1 + i + 1] - src_offset1[s1_idx_m1 + i - 1])) * _s2_inc2_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + i + 1] - src_offset0[s1_idx_p1 + i - 1]) -
                                       (src_offset1[s1_idx_p1 + i + 1] - src_offset1[s1_idx_p1 + i - 1])) * _s2_inc2_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T corre       = (_half_sigma_s1_sq * s1 * s1 * gamma_s1) + (_r_s1 * s1 * delta_s1) - rv_s1 +
                                      (_half_sigma_s2_sq * s2 * s2 * gamma_s2) + (_r_s2 * s2 * delta_s2) - rv_s2 +
                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * corr_s1_s2 * x_gamma);
                dst_offset[s1_idx + i]   = (corre * _t_inc * _theta) + dst_offset[s1_idx + i];
            }
    
            /* Apply the s2 = 0 boundary condition */
            //dst_offset[s1_idx] = (2.0 * dst_offset[s1_idx + 1]) - dst_offset[s1_idx + 2];
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
                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * corr_s1_s2 * x_gamma);
                dst_offset[s1_idx]   = (corre * _t_inc * _theta) + dst_offset[s1_idx];
            }
        
            /* Apply the s2 = _s2_max boundary condition */
            //dst_offset[s1_idx_p1 - 1] = (2.0 * dst_offset[s1_idx_p1 - 2]) - dst_offset[s1_idx_p1 - 3];
            {
                /* S1 greeks */
                const T delta_s1    = ((src_offset0[s1_idx_p1 + size] - src_offset0[s1_idx_m1 + size]) -
                                       (src_offset1[s1_idx_p1 + size] - src_offset1[s1_idx_m1 + size])) * _s1_inc2_inv;
                const T rv_s1       = (src_offset0[s1_idx + size] - src_offset1[s1_idx + size]) * _r_s1;

                /* s2 greeks */
                const T s2          = static_cast<T>(size) * _s2_inc;;
                const T delta_s2    = ((src_offset0[s1_idx + size] - src_offset0[s1_idx + size - 1]) -
                                       (src_offset1[s1_idx + size] - src_offset1[s1_idx + size - 1])) * s2_inc_inv;
                const T rv_s2       = (src_offset0[s1_idx + size] - src_offset1[s1_idx + size]) * _r_s2;
                
                /* Cross partial */
                const T delta_m1    = ((src_offset0[s1_idx_m1 + size] - src_offset0[s1_idx_m1 + size - 1]) -
                                       (src_offset1[s1_idx_m1 + size] - src_offset1[s1_idx_m1 + size - 1])) * s2_inc_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + size] - src_offset0[s1_idx_p1 + size - 1]) -
                                       (src_offset1[s1_idx_p1 + size] - src_offset1[s1_idx_p1 + size - 1])) * s2_inc_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T corre       = (_r_s1 * s1 * delta_s1) - rv_s1 +
                                      (_r_s2 * s2 * delta_s2) - rv_s2 +
                                      (_sigma_s1 * s1 * _sigma_s2 * s2 * corr_s1_s2 * x_gamma);
                dst_offset[s1_idx_p1 - 1]   = (corre * _t_inc * _theta) + dst_offset[s1_idx_p1 - 1];
            }
        }

    private :
        const T     _theta;
        const T     _t_inc;
        const T     _s1_inc;
        const T     _s2_inc;
        const T     _r_s1;
        const T     _r_s2;
        const T     _sigma_s1;
        const T     _sigma_s2;
        const T     corr_s1_s2;
        const T     s2_inc_inv;
        const T     _s1_inc2_inv;
        const T     _s2_inc2_inv;
        const T     _s1_inc_sq_inv;
        const T     _s2_inc_sq_inv;
        const T     _half_sigma_s1_sq;
        const T     _half_sigma_s2_sq;
        const int   size; 
};


template<class T>
class cross_term_corrector
{
    public :
        /* CTOR */
        cross_term_corrector(const T theta, const T t_inc, const T s1_inc, const T s2_inc,
            const T sigma_s1, const T sigma_s2, const T corr_s1_s2, const int size)
            : _theta(theta), _t_inc(t_inc), _s1_inc(s1_inc), _s2_inc(s2_inc), 
              _sigma_s1(sigma_s1), _sigma_s2(sigma_s2), corr_s1_s2(corr_s1_s2), 
              _s1_inc2_inv(0.5 / s1_inc), _s2_inc2_inv(0.5 / s2_inc), size(size)
        { }

        /* DTOR */
        ~cross_term_corrector() { }
        
        /* Cross correction algorithm of the form F(A) = A + _theta * _t_inc * (F(B) - F(C)). 
           F is a cross term in a multi-asset BS PDE. A, B and C may be any data sets 
           to correct from */
        void operator()(const T *const src_offset0, const T *const src_offset1, T *const dst_offset, const int step)
        {
            const int s1_idx    = step * size;
            const int s1_idx_p1 = s1_idx + size;
            const int s1_idx_m1 = s1_idx - size;
            const T s1          = static_cast<T>(step) * _s1_inc;
                
            /* Fill in the column */
            for (int i = 1; i < size - 1; i++)
            {
                /* Cross partial */
                const T s2          = static_cast<T>(i) * _s2_inc;
                const T delta_m1    = ((src_offset0[s1_idx_m1 + i + 1] - src_offset0[s1_idx_m1 + i - 1]) -
                                       (src_offset1[s1_idx_m1 + i + 1] - src_offset1[s1_idx_m1 + i - 1])) * _s2_inc2_inv;
                const T delta_p1    = ((src_offset0[s1_idx_p1 + i + 1] - src_offset0[s1_idx_p1 + i - 1]) -
                                       (src_offset1[s1_idx_p1 + i + 1] - src_offset1[s1_idx_p1 + i - 1])) * _s2_inc2_inv;
                const T x_gamma     = (delta_p1 - delta_m1) * _s1_inc2_inv;

                /* Update */
                const T cross           = _sigma_s1 * s1 * _sigma_s2 * s2 * corr_s1_s2 * x_gamma;
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
//                const T cross           = _sigma_s1 * s1 * _sigma_s2 * s2 * corr_s1_s2 * x_gamma;
//                dst_offset[s1_idx + i]  = (cross * _theta * _t_inc) + dst_offset[s1_idx + i];
//            }
        
            /* Apply the s2 = _s2_max boundary condition */
            //dst_offset[s1_idx_p1 - 1] = (2.0 * dst_offset[s1_idx_p1 - 2]) - dst_offset[s1_idx_p1 - 3];
        }

    private :
        const T     _theta;
        const T     _t_inc;
        const T     _s1_inc;
        const T     _s2_inc;
        const T     _sigma_s1;
        const T     _sigma_s2;
        const T     corr_s1_s2;
        const T     _s1_inc2_inv;
        const T     _s2_inc2_inv;
        const int   size;       
};

#endif
