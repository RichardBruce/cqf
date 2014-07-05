#ifndef __DOUGLAS_H__
#define __DOUGLAS_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>


#include "matrix.h"
#include "dumping.h"
#include "pde_solver.h"


template<class T>
class douglas_pde_solver : public pde_solver<T>
{
    public :
        douglas_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<grid<T>*> &grids,
            const boundary_condition<T> &bc,
            const std::string &log_file, 
            const T max_dt, bool log_final_only)
            : pde_solver<T>(cashflows, economics, grids, bc, log_file, max_dt, log_final_only)
            {
                assert(this->_grids.size() == 1);
            };
         
        ~douglas_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T douglas_pde_solver<T>::solve() const
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
    int tp1_offset      = nas;
    int t_offset        = 0;
    bool use_rannacher  = false;
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        if (this->_cashflows[i]->exists(this->_t_steps->back()))
        {
            use_rannacher |= this->_cashflows[i]->add(&fd_grid[tp1_offset], this->_grids, this->_t_steps->back());
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
    
    /* Work backwards in time */
    T *matrix_data_k    = new T [nas * nas];
    T *matrix_data_kp1  = new T [nas];
    for (int i = this->_t_steps->size() - 1; i > 0; i--)
    {
        memset(matrix_data_k, 0, (nas * nas * sizeof(T)));
        const T t       = (*this->_t_steps)[i - 1];
        const T t_inc   = (*this->_t_steps)[i] - t;
        const T prop    = 0.0;//(s_inc * s_inc) / (12.0 * t_inc);
        const T theta   = (0.5 + prop) * t_inc;
        const T m_theta = (0.5 - prop) * t_inc;
        assert(theta > 0.0);
        assert(theta < 1.0);
        if (use_rannacher)
        {
            implicit_step<double>(&fd_grid[tp1_offset], &fd_grid[t_offset], 
                matrix_data_k, grid, t_inc * 0.5, r, half_sigma_sq, sor_w, sor_tol, 
                sor_max_it, 2, use_sor);
            use_rannacher = false;
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

                const T a = (coeffs[0] * r_s) + (coeffs[3] * g);
                const T b = (coeffs[1] * r_s) + (coeffs[4] * g) - r;
                const T c = (coeffs[2] * r_s) + (coeffs[5] * g);
            
                const int matrix_offset = (j * nas) + j;
                matrix_data_k[matrix_offset - 1] = -a * theta;
                matrix_data_k[matrix_offset    ] = 1.0 - b * theta;
                matrix_data_k[matrix_offset + 1] = -c * theta;
            
                matrix_data_kp1[j]  = a * m_theta * fd_grid[tp1_offset + j - 1];
                matrix_data_kp1[j] += (1.0 + b * m_theta) * fd_grid[tp1_offset + j];
                matrix_data_kp1[j] += c * m_theta * fd_grid[tp1_offset + j + 1];
            }
        
            /* Boundary conditions */
            /* s = 0.0 */
            {
                const T b = -r;
                matrix_data_k[0] = 1.0 - b * theta;
                matrix_data_k[1] = 0.0;
                
                matrix_data_kp1[0] = (1.0 + b * m_theta) * fd_grid[tp1_offset];
            }
       
            /* s = s_max*/
            {
                const T s   = grid.points()[nas - 1];
                const T r_s = r * s;

                const T a = -r_s;
                const T b = -(r - r_s);
                
                const int matrix_offset = (nas - 1) * (nas + 1);
                matrix_data_k[matrix_offset - 1] = -a * theta;
                matrix_data_k[matrix_offset    ] = 1.0 - b * theta;
                
                matrix_data_kp1[nas - 1]  = a * m_theta * fd_grid[tp1_offset + nas - 2];
                matrix_data_kp1[nas - 1] += (1.0 + b * m_theta) * fd_grid[tp1_offset + nas - 1];
            }
        
            matrix<T> matrix_kp1(&matrix_data_kp1[0], nas, 1, false, true);
            matrix<T> matrix_k(&matrix_data_k[0], nas, nas, false, true);
            matrix_k.gauss_solve(matrix_kp1, &fd_grid[t_offset], Tridiagonal);
        }
        
        /* Process cashflows */
        for (unsigned int i = 0; i < this->_cashflows.size(); i++)
        {
            if (this->_cashflows[i]->exists(t))
            {
                use_rannacher |= this->_cashflows[i]->add(&fd_grid[t_offset], this->_grids, t);
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
    delete [] matrix_data_kp1;
    delete [] fd_grid;
    
    return result;
}

#endif
