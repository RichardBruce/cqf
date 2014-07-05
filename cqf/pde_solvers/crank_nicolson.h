#ifndef __CRANK_NICOLSON_H__
#define __CRANK_NICOLSON_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>


#include "matrix.h"
#include "pde_solver.h"
#include "implicit.h"


template<class T>
class crank_nicolson_pde_solver : public pde_solver<T>
{
    public :
        crank_nicolson_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<risks<T>*> &risks,
            const boundary_condition<T> &bc,
            const std::string &log_file, 
            const T max_dt, bool log_final_only)
            : pde_solver<T>(cashflows, economics, risks, bc, log_file, max_dt, log_final_only)
            {
                assert(this->_grids.size() == 1);
            };
         
        ~crank_nicolson_pde_solver() { };
        
        T iteration(T *const t_offset, const T *const tp1_offset, const T t_inc, bool disc_cashflow) const;

    private :
        T *const            matrix_data_k;
        static const bool   use_sor         = false;
        static const int    sor_max_it      = 5000;
        static const T      sor_w           = 1.85;
        static const T      sor_tol         = 1e-6;
};

    
template<class T>    
T crank_nicolson_pde_solver<T>::iteration(T *const t_offset, const T *const tp1_offset, const T t_inc, bool disc_cashflow) const
{    
    memset(matrix_data_k, 0, (nas * nas * sizeof(T)));
    if (disc_cashflow)
    {
        implicit_step<double>(tp1_offset, t_offset, 
            matrix_data_k, grid, t_inc * 0.5, r, half_sigma_sq, sor_w, sor_tol, 
            sor_max_it, 2, use_sor);
    }
    else
    {
        /* Populate the matrix */
        for (int j = 1; j < nas - 1; j++)
        {
            const T s = grid.points()[j];
            const T g = half_sigma_sq * s * s;
            const T r_s = r * s;
        
            T coeffs[6];
            grid.get_delta_coeffs(&coeffs[0], j);
            grid.get_gamma_coeffs(&coeffs[3], j);

            const T a = ((coeffs[0] * r_s) + (coeffs[3] * g)) * 0.5 * t_inc;
            const T b = ((coeffs[1] * r_s) + (coeffs[4] * g) - r) * 0.5 * t_inc;
            const T c = ((coeffs[2] * r_s) + (coeffs[5] * g)) * 0.5 * t_inc;
        
            const int matrix_offset = (j * nas) + j;
            matrix_data_k[matrix_offset - 1] = -a;
            matrix_data_k[matrix_offset    ] = 1.0 - b;
            matrix_data_k[matrix_offset + 1] = -c;
        
            matrix_data_kp1[j]  = a * tp1_offset[j - 1];
            matrix_data_kp1[j] += (1.0 + b) * tp1_offset[j];
            matrix_data_kp1[j] += c * tp1_offset[j + 1];
        }

        /* Boundary conditions */
        /* s = 0.0 */
        {
            const T b = -r * 0.5 * t_inc;
            matrix_data_k[0] = 1.0 - b;
            matrix_data_k[1] = 0.0;
            
            matrix_data_kp1[0] = (1.0 + b) * tp1_offset[0];
        }

        /* s = s_max*/
        {
            const T s   = grid.points()[nas - 1];
            const T r_s = r * s;

            const T a = -r_s * 0.5 * t_inc;
            const T b = -(r - r_s) * 0.5 * t_inc;
            
            const int matrix_offset = (nas - 1) * (nas + 1);
            matrix_data_k[matrix_offset - 1] = -a;
            matrix_data_k[matrix_offset    ] = 1.0 - b;
            
            matrix_data_kp1[nas - 1]  = a * tp1_offset[nas - 2];
            matrix_data_kp1[nas - 1] += (1.0 + b) * tp1_offset[nas - 1];
        }

        matrix<T> matrix_kp1(&matrix_data_kp1[0], nas, 1, false, true);
        matrix<T> matrix_k(&matrix_data_k[0], nas, nas, false, true);
        matrix_k.gauss_solve(matrix_kp1, t_offset[0], Tridiagonal);
    }
}

#endif
