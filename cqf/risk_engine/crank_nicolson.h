#ifndef __CRANK_NICOLSON_H__
#define __CRANK_NICOLSON_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <assert.h>

#include "matrix.h"

#include "pde_solver.h"
#include "implicit.h"

#include "scenario_point.h"


template<class T>
class crank_nicolson_pde_solver : public pde_solver<T>
{
    public :
        crank_nicolson_pde_solver(const bool use_sor = false, const int sor_max_it = 5000, 
            const T sor_w = 1.85, const T sor_tol = 1e-6)
                : _use_sor(use_sor), _sor_max_it(sor_max_it), _sor_w(sor_w), _sor_tol(sor_tol) { };
       
        void solve(scenario_point<T> &p);

    private :
        std::vector<T>  _matrix_data_k;
        std::vector<T>  _matrix_data_kp1;
        const bool      _use_sor;
        const int       _sor_max_it;
        const T         _sor_w;
        const T         _sor_tol;
};

    
template<class T>    
void crank_nicolson_pde_solver<T>::solve(scenario_point<T> &p)
{    
    grid<T> &grid = p.get_grid();
    assert((grid.dimensions() == 1) || (grid.dimensions() == 2 && grid.outer_grid_disconnected()));
    
    /* Resize matrix data */
    const int nas = grid.size(grid.dimensions() - 1);
    if (_matrix_data_kp1.size() < static_cast<unsigned int>(nas))
    {
        _matrix_data_k.resize(nas * nas);
        _matrix_data_kp1.resize(nas);
    }
    
    /* This is needed in case 2 grids of different sizes are processed */
    /* When this happens some off diagonal elements arent set in the below loop */
    /* These elements probably shouldnt be pick up by the solve, to be investigated */
    /* Should now be fixed, lets test */
//    memset(&_matrix_data_k[0], 0, nas * nas * sizeof(T));
    
    /* Get common variables and precomputes */
    const T r               = p.get_interest_rate();
    const T t_inc           = p.get_time_increment();
    const T sigma           = p.get_volatility();
    const T half_sigma_sq   = 0.5 * sigma * sigma;
    const int dim_m1        = grid.dimensions() - 1;            /* Always the regular grid */
    
    /* Process discontinious cashflows with implicit */
    if (p.is_disc_cashflow())
    {
        implicit_step<double>(_matrix_data_k, grid, t_inc * 0.5, r, half_sigma_sq, 
            _sor_w, _sor_tol, dim_m1, _sor_max_it, 2, _use_sor);
    }
    else
    {
        T coeffs[6];
        const int nr_of_grids = (dim_m1 == 0) ? 1 : grid.size(0);
        
        /* Process disconnected grids one at a time */
        for (int i = 0; i < nr_of_grids; i++)
        {
            /* Populate the matrix */
            const int dim_offset = i * grid.size(dim_m1);
            for (int j = 1; j < nas - 1; j++)
            {
                const T s = grid.s_grid(dim_m1)[j];
                const T g = half_sigma_sq * s * s;
                const T r_s = r * s;
        
                grid.get_delta_coeffs_middle(&coeffs[0], dim_m1, j);
                grid.get_gamma_coeffs_middle(&coeffs[3], dim_m1, j);
                const T a = ((coeffs[0] * r_s) + (coeffs[3] * g)) * 0.5 * t_inc;
                const T b = ((coeffs[1] * r_s) + (coeffs[4] * g) - r) * 0.5 * t_inc;
                const T c = ((coeffs[2] * r_s) + (coeffs[5] * g)) * 0.5 * t_inc;
        
                const int matrix_offset = (j * nas) + j;
                _matrix_data_k[matrix_offset - 1] = -a;
                _matrix_data_k[matrix_offset    ] = 1.0 - b;
                _matrix_data_k[matrix_offset + 1] = -c;
            
                _matrix_data_kp1[j]  = a * grid.read_v_grid()[dim_offset + j - 1];
                _matrix_data_kp1[j] += (1.0 + b) * grid.read_v_grid()[dim_offset + j];
                _matrix_data_kp1[j] += c * grid.read_v_grid()[dim_offset + j + 1];
            }
            /* Boundary conditions */
            /* s = s_min */
            {        
                const std::pair<bool, T> &b = grid.get_value_coeffs_lower();
                const int matrix_offset = 0;
                if (b.first)
                {
                    _matrix_data_k[matrix_offset    ] = 1.0;
                    _matrix_data_k[matrix_offset + 1] = 0.0;
                    _matrix_data_k[matrix_offset + 2] = 0.0;
        
                    _matrix_data_kp1[0] = b.second;
                }
                else
                {
                    grid.get_delta_coeffs_lower(&coeffs[0], dim_m1);
    
                    const T s = grid.s_grid(dim_m1)[0];
                    const T r_s = r * s;
            
                    const T a = ((coeffs[0] * r_s)) * 0.5 * t_inc;
                    const T b = ((coeffs[1] * r_s) - r) * 0.5 * t_inc;
                    const T c = ((coeffs[2] * r_s)) * 0.5 * t_inc;
        
                    _matrix_data_k[matrix_offset    ] = (1.0 - a);
                    _matrix_data_k[matrix_offset + 1] = -b;
                    _matrix_data_k[matrix_offset + 2] = -c;
        
                    _matrix_data_kp1[0]  = (1.0 + a) * grid.read_v_grid()[dim_offset    ];
                    _matrix_data_kp1[0] += b         * grid.read_v_grid()[dim_offset + 1];
                    _matrix_data_kp1[0] += c         * grid.read_v_grid()[dim_offset + 2];
                }
            }

            /* s = s_max*/
            {
                const std::pair<bool, T> &b = grid.get_value_coeffs_upper();
                const int matrix_offset = ((nas - 1) * nas) + (nas - 1);
                if (b.first)
                {
                    _matrix_data_k[matrix_offset - 2] = 0.0;
                    _matrix_data_k[matrix_offset - 1] = 0.0;
                    _matrix_data_k[matrix_offset    ] = 1.0;
        
                    _matrix_data_kp1[nas - 1] = b.second;
                }
                else
                {
                    grid.get_delta_coeffs_upper(&coeffs[0], dim_m1);
            
                    const T s = grid.s_grid(dim_m1)[nas - 1];
                    const T r_s = r * s;
            
                    const T a = ((coeffs[0] * r_s)) * 0.5 * t_inc;
                    const T b = ((coeffs[1] * r_s) - r) * 0.5 * t_inc;
                    const T c = ((coeffs[2] * r_s)) * 0.5 * t_inc;
        
                    _matrix_data_k[matrix_offset - 2] = -a;
                    _matrix_data_k[matrix_offset - 1] = -b;
                    _matrix_data_k[matrix_offset    ] = (1.0 - c);
        
                    _matrix_data_kp1[nas - 1]  = a         * grid.read_v_grid()[dim_offset + nas - 3];
                    _matrix_data_kp1[nas - 1] += b         * grid.read_v_grid()[dim_offset + nas - 2];
                    _matrix_data_kp1[nas - 1] += (1.0 + c) * grid.read_v_grid()[dim_offset + nas - 1];
                }
            }
        
            matrix<T> matrix_kp1(&_matrix_data_kp1[0], nas, 1, false, true);
            matrix<T> matrix_k(&_matrix_data_k[0], nas, nas, false, true);
            matrix_k.gauss_solve(matrix_kp1, &grid.write_v_grid()[dim_offset], TridiagonalPlus);
        }
    }
}

#endif
