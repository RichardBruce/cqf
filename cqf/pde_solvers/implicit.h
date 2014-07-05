#ifndef __IMPLICIT_H__
#define __IMPLICIT_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>

#include "matrix.h"
#include "pde_solver.h"


template<class T>
void implicit_step(T *const src, T *const dst, T *const mat_dat, const grid<T> &grid, 
    const T t_inc, const T r, const T half_sigma_sq, const T sor_w, const T sor_tol, 
    const int sor_max_it, const int nts, const bool use_sor)
{
    /* Run multiple steps into the same result */
    const int nas = grid.points().size();
    for (int i = 0; i < nts; i++)
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

            const T a = (coeffs[0] * r_s) + (coeffs[3] * g);
            const T b = (coeffs[1] * r_s) + (coeffs[4] * g) - r;
            const T c = (coeffs[2] * r_s) + (coeffs[5] * g);
            
            const int matrix_offset = (j * nas) + j;
            mat_dat[matrix_offset - 1] = -a * t_inc;
            mat_dat[matrix_offset    ] = 1.0 - b * t_inc;
            mat_dat[matrix_offset + 1] = -c * t_inc;
        }
        
        /* Boundary conditions */
        /* s = 0.0 */
        {
            mat_dat[0] = 1.0 + r * t_inc;
            mat_dat[1] = 0.0;
        }
       
        /* s = s_max*/
        {
            const T s = grid.points()[nas - 1];
            const int matrix_offset = (nas - 1) * (nas + 1);
            mat_dat[matrix_offset - 1] = s * r * t_inc;
            mat_dat[matrix_offset    ] = 1.0 + (r -s * r) * t_inc;
        }
        
        /* Solve the matrix */
        matrix<T> mat(mat_dat, nas, nas, false, true);
        if (use_sor)
        {
            mat.sor_solve(&src[0], &dst[0], sor_w, sor_tol, sor_max_it, Tridiagonal);
        }
        else
        {
            mat.gauss_solve(&src[0], &dst[0], Tridiagonal);
        }
        
        /* Repeated iterations into dst */
        src = dst;
    }
}


template<class T>
void implicit_step(T *src, T *dst, const T t_inc, const T r, const T sigma_sq,
    const T sor_w, const T sor_tol, const int sor_max_it, 
    const int nts, const int nas, const bool use_sor)
{
    /* Run multiple steps into the same result */
    for (int i = 0; i < nts; i++)
    {
        /* Populate the matrix */
        T *matrix_data = new T [nas * nas];
        for (int j = 1; j < nas - 1; j++)
        {
            const int matrix_offset = (j * nas) + j;
            const T flt_j = static_cast<T>(j);
            matrix_data[matrix_offset - 1] = 0.5 * flt_j * (r - sigma_sq * flt_j) * t_inc;
            matrix_data[matrix_offset    ] = 1.0 + (r + sigma_sq * flt_j * flt_j) * t_inc;
            matrix_data[matrix_offset + 1] = 0.5 * flt_j * (-r - sigma_sq * flt_j) * t_inc;
        }
        
        /* Boundary conditions */
        /* s = 0.0 */
        {
            matrix_data[0] = 1.0 + r * t_inc;
            matrix_data[1] = 0.0;
        }
       
        /* s = s_max*/
        {
            const T flt_j = static_cast<T>(nas - 1);
            const int matrix_offset = (nas - 1) * (nas + 1);
            matrix_data[matrix_offset - 1] = flt_j * r * t_inc;
            matrix_data[matrix_offset    ] = 1.0 + (r -flt_j * r) * t_inc;
        }
        
        /* Solve the matrix */
        matrix<T> mat(matrix_data, nas, nas, true);
        if (use_sor)
        {
            mat.sor_solve(&src[0], &dst[0], sor_w, sor_tol, sor_max_it, Tridiagonal);
        }
        else
        {
            mat.gauss_solve(&src[0], &dst[0], Tridiagonal);
        }
        
        /* Repeated iterations into dst */
        src = dst;
    }
}


template<class T>
class implicit_pde_solver : public pde_solver<T>
{
    public :
        implicit_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<grid<T>*> &grids,
            const boundary_condition<T> &bc,
            const std::string &log_file, 
            const T max_dt, bool log_final_only)
            : pde_solver<T>(cashflows, economics, grids, bc, log_file, max_dt, log_final_only)
            {
                assert(this->_grids.size() == 1);
            };
         
        ~implicit_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T implicit_pde_solver<T>::solve() const
{
    /* Matrix solver params */
    bool use_sor    = false;
    int sor_max_it  = 5000;
    T sor_w         = 1.85;
    T sor_tol       = 1e-6;
   
    /* Economics */
    const T r               = this->_economics[0]->get_interest_rate();
    const T sigma           = this->_economics[0]->get_volatility();
    const T half_sigma_sq   = 0.5 * sigma * sigma;
    
    /* Allocate grid for V */
    const grid<T> &grid = (*(this->_grids[0]));
    const int nas       = grid.size();
    T *fd_grid          = new T[2 * nas];
    
    /* Set the end condition (ie payoff) */
    int tp1_offset = nas;
    int t_offset   = 0;
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        if (this->_cashflows[i]->exists(this->_t_steps->back()))
        {
            this->_cashflows[i]->add(&fd_grid[tp1_offset], this->_grids, this->_t_steps->back());
        }
    }
    
    /* Logging */
    std::ofstream file;
    if (!this->_log_file.empty())
    {
        file.open(this->_log_file.c_str());
        assert(file.is_open());
        dump_row_to_gnuplot(file, &fd_grid[tp1_offset], grid.points(), this->_t_steps->back());
    }
        
    T *matrix_data_k = new T [nas * nas];
    for (int i = this->_t_steps->size() - 1; i > 0; i--)
    {
        memset(matrix_data_k, 0, (nas * nas * sizeof(T)));
        const T t       = (*this->_t_steps)[i - 1];
        const T t_inc   = (*this->_t_steps)[i] - t;
        implicit_step<double>(&fd_grid[tp1_offset], &fd_grid[t_offset], 
                matrix_data_k, grid, t_inc, r, half_sigma_sq, sor_w, sor_tol, 
                sor_max_it, 1, use_sor);
        
        /* Process cashflows */
        for (unsigned int i = 0; i < this->_cashflows.size(); i++)
        {
            if (this->_cashflows[i]->exists(t))
            {
                this->_cashflows[i]->add(&fd_grid[t_offset], this->_grids, t);
            }
        }
    
        /* Logging */
        if (!this->_log_file.empty())
        {
            dump_row_to_gnuplot(file, &fd_grid[t_offset], grid.points(), t);
        }
        
        std::swap(t_offset, tp1_offset);
    }

    /* Logging */
    if (!this->_log_file.empty())
    {
        file.close();
    }
    
    /* Get result */
    const T result = grid.interpolate_result(&fd_grid[tp1_offset], this->_economics[0]->get_spot());
    
    /* Clean up */
    delete [] matrix_data_k;
    delete [] fd_grid;
    
    return result;
}

#endif
