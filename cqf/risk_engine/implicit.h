#ifndef __IMPLICIT_H__
#define __IMPLICIT_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <assert.h>

#include "matrix.h"

#include "pde_solver.h"

#include "grid_mapper.h"


template<class T>
void implicit_step(std::vector<T> &mat_dat, grid<T> &grid, 
    const T t_inc, const T r, const T half_sigma_sq, const T sor_w, const T sor_tol, 
    const int dim_m1, const int sor_max_it, const int nts, const bool use_sor)
{
    assert((grid.dimensions() == 1) || (grid.dimensions() == 2 && grid.outer_grid_disconnected()));

    /* Run multiple steps into the same result */
    T coeffs[6];
    const int nas = grid.size(dim_m1);
    const int nr_of_grids = (dim_m1 == 0) ? 1 : grid.size(0);
    for (int i = 0; i < nts; i++)
    {
        for (int j = 0; j < nr_of_grids; j++)
        {
            /* Populate the matrix */
            for (int k = 1; k < nas - 1; k++)
            {
                const T s = grid.s_grid(dim_m1)[k];
                const T g = half_sigma_sq * s * s;
                const T r_s = r * s;
             
                grid.get_delta_coeffs_middle(&coeffs[0], dim_m1, k);
                grid.get_gamma_coeffs_middle(&coeffs[3], dim_m1, k);
    
                const T a = (coeffs[0] * r_s) + (coeffs[3] * g);
                const T b = (coeffs[1] * r_s) + (coeffs[4] * g) - r;
                const T c = (coeffs[2] * r_s) + (coeffs[5] * g);
                
                const int matrix_offset = (k * nas) + k;
                mat_dat[matrix_offset - 1] = -a * t_inc;
                mat_dat[matrix_offset    ] = 1.0 - b * t_inc;
                mat_dat[matrix_offset + 1] = -c * t_inc;
            }
            /* Boundary conditions */
            /* s = s_min */
            {
                const std::pair<bool, T> &b = grid.get_value_coeffs_lower();
                const int matrix_offset = 0;
                if (b.first)
                {
                    mat_dat[matrix_offset    ] = 1.0;
                    mat_dat[matrix_offset + 1] = 0.0;
                    mat_dat[matrix_offset + 2] = 0.0;
                }
                else
                {            
                    const T s = grid.s_grid(dim_m1)[0];
                    const T r_s = r * s;
            
                    grid.get_delta_coeffs_lower(&coeffs[0], dim_m1);
                    const T a = ((coeffs[0] * r_s)) * t_inc;
                    const T b = ((coeffs[1] * r_s) - r) * t_inc;
                    const T c = ((coeffs[2] * r_s)) * t_inc;
            
                    mat_dat[matrix_offset    ] = (1.0 - a);
                    mat_dat[matrix_offset + 1] = -b;
                    mat_dat[matrix_offset + 2] = -c;
                }
            }
            /* s = s_max*/
            {
                const std::pair<bool, T> &b = grid.get_value_coeffs_upper();
                const int matrix_offset = ((nas - 1) * nas) + (nas - 1);
                if (b.first)
                {
                    mat_dat[matrix_offset - 2] = 0.0;
                    mat_dat[matrix_offset - 1] = 0.0;
                    mat_dat[matrix_offset    ] = 1.0;
                }
                else
                {            
                    const T s = grid.s_grid(dim_m1)[nas - 1];
                        const T r_s = r * s;
            
                    grid.get_delta_coeffs_upper(&coeffs[0], dim_m1);
                    const T a = ((coeffs[0] * r_s)) * t_inc;
                    const T b = ((coeffs[1] * r_s) - r) * t_inc;
                    const T c = ((coeffs[2] * r_s)) * t_inc;
            
                    mat_dat[matrix_offset - 2] = -a;
                    mat_dat[matrix_offset - 1] = -b;
                    mat_dat[matrix_offset    ] = (1.0 - c);
                }
            }
            
            /* Solve the matrix */
            const int dim_offset = j * grid.size(0);
            matrix<T> mat(mat_dat.data(), nas, nas, false, true);
            if (use_sor)
            {
                mat.sor_solve(&grid.read_v_grid()[dim_offset], &grid.write_v_grid()[dim_offset], sor_w, sor_tol, sor_max_it, Tridiagonal);
            }
            else
            {
                mat.gauss_solve(&grid.read_v_grid()[dim_offset], &grid.write_v_grid()[dim_offset], Tridiagonal);
            }
        }
            
        /* Repeated iterations into dst */
        if (i != (nts - 1))
        {
            grid.flip();
        }
    }
}



//template<class T>
//class implicit_pde_solver : public pde_solver<T>
//{
//    public :
//        implicit_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
//            const std::vector<equity_economics<T>*> &economics,
//            const std::vector<grid<T>*> &grids,
//            const boundary_condition<T> &bc,
//            const std::string &log_file, 
//            const T max_dt, bool log_final_only)
//            : pde_solver<T>(cashflows, economics, grids, bc, log_file, max_dt, log_final_only)
//            {
//                assert(this->_grids.size() == 1);
//            };
//         
//        ~implicit_pde_solver() { };
//        
//        T solve() const;
//};
//
//    
//template<class T>    
//T implicit_pde_solver<T>::solve() const
//{
//    /* Matrix solver params */
//    bool use_sor    = false;
//    int sor_max_it  = 5000;
//    T sor_w         = 1.85;
//    T sor_tol       = 1e-6;
//   
//    /* Economics */
//    const T r               = this->_economics[0]->get_interest_rate();
//    const T sigma           = this->_economics[0]->get_volatility();
//    const T half_sigma_sq   = 0.5 * sigma * sigma;
//    
//    /* Allocate grid for V */
//    const grid<T> &grid = (*(this->_grids[0]));
//    const int nas       = grid.size();
//    T *fd_grid          = new T[2 * nas];
//    
//    /* Set the end condition (ie payoff) */
//    int tp1_offset = nas;
//    int t_offset   = 0;
//    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
//    {
//        if (this->_cashflows[i]->exists(this->_t_steps->back()))
//        {
//            this->_cashflows[i]->add(&fd_grid[tp1_offset], this->_grids, this->_t_steps->back());
//        }
//    }
//    
//    /* Logging */
//    std::ofstream file;
//    if (!this->_log_file.empty())
//    {
//        file.open(this->_log_file.c_str());
//        assert(file.is_open());
//        dump_row_to_gnuplot(file, &fd_grid[tp1_offset], grid.points(), this->_t_steps->back());
//    }
//        
//    T *matrix_data_k = new T [nas * nas];
//    for (int i = this->_t_steps->size() - 1; i > 0; i--)
//    {
//        memset(matrix_data_k, 0, (nas * nas * sizeof(T)));
//        const T t       = (*this->_t_steps)[i - 1];
//        const T t_inc   = (*this->_t_steps)[i] - t;
//        implicit_step<double>(&fd_grid[tp1_offset], &fd_grid[t_offset], 
//                matrix_data_k, grid, t_inc, r, half_sigma_sq, sor_w, sor_tol, 
//                sor_max_it, 1, use_sor);
//        
//        /* Process cashflows */
//        for (unsigned int i = 0; i < this->_cashflows.size(); i++)
//        {
//            if (this->_cashflows[i]->exists(t))
//            {
//                this->_cashflows[i]->add(&fd_grid[t_offset], this->_grids, t);
//            }
//        }
//    
//        /* Logging */
//        if (!this->_log_file.empty())
//        {
//            dump_row_to_gnuplot(file, &fd_grid[t_offset], grid.points(), t);
//        }
//        
//        std::swap(t_offset, tp1_offset);
//    }
//
//    /* Logging */
//    if (!this->_log_file.empty())
//    {
//        file.close();
//    }
//    
//    /* Get result */
//    const T result = grid.interpolate_result(&fd_grid[tp1_offset], this->_economics[0]->get_spot());
//    
//    /* Clean up */
//    delete [] matrix_data_k;
//    delete [] fd_grid;
//    
//    return result;
//}

#endif
