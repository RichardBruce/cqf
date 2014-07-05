#ifndef __ADE_2D_H__
#define __ADE_2D_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>

#include <boost/lexical_cast.hpp>

#include "pde_solver.h"
#include "correctors.h"
#include "matrix.h"
#include "dumping.h"

#include "cashflow.h"


template<class T>
class ade_2d_pde_solver : public pde_solver_2d<T>
{
    public :
        ade_pde_solver_2d(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const T corr, const int nts, 
            const int nas1, const int nas2, const T t, bool log_final_only)
            : pde_solver_2d<T>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only),
            _u_grid(new T[nas1 * nas2]), _v_grid(new T[nas1 * nas2]) { };
         
        ~ade_pde_solver_2d()
        {
            delete [] _u_grid;
            delete [] _v_grid;
        };
        
        T solve() const;

    private :
        T *_u_grid;
        T *_v_grid;
};

    
template<class T>    
T ade_pde_solver<T>::solve() const
{
    /* Economics */
    const T r_s1        = this->_economics[0]->get_interest_rate();
    const T sigma_s1    = this->_economics[0]->get_volatility();
    
    const T r_s2        = this->_economics[1]->get_interest_rate();
    const T sigma_s2    = this->_economics[1]->get_volatility();
    
    /* Grid set up */
    const T t_inc   = this->_t / static_cast<T>(this->_nts);
    const int nsn   = this->_nas1 * this->_nas2 * 2;
    T *fd_grid      = new T[nsn];
    
    const T s1_max      = 300.0;
    const T s1_inc      = s1_max / static_cast<T>(this->_nas1);
    const T s1_inc_inv  = 1.0 / s1_inc;
    
    const T s2_max      = 300.0;
    const T s2_inc      = s2_max / static_cast<T>(this->_nas2);
    const T s2_inc_inv  = 1.0 / s2_inc;

        
    /* Set the end condition (ie payoff) */
    int tp1_offset    = nsn;
    int t_offset      = 0; 
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[tp1_offset], this->_t);
    }

    /* Work backwards in time */
    const T lambda_s1 = t_inc / (s1_inc * s1_inc);
    const T lambda_s2 = t_inc / (s2_inc * s2_inc);
    for (int i = 0; i < this->_nts; i++)
    {
            /* Lower Boundary Condition*/
            this->_u_grid[0] = exp(-r_s1 * t_inc) * fd_grid[tp1_offset];
            this->_v_grid[0] = this->_u_grid[0];
        
            /* Upper Boundary Condition */
            {
                const T s       = static_cast<T>(this->_nas1) * s1_inc;
                const T delta   = (fd_grid[tp1_offset + this->_nas1 - 1] - fd_grid[tp1_offset + this->_nas1 - 2]) * s1_inc_inv;
                const T rv      = r_s1 * fd_grid[tp1_offset + this->_nas1 - 1];
                const T theta   = (r_s1 * s * delta) - rv;            

                this->_u_grid[this->_nas1 - 1] = (theta * t_inc) + fd_grid[tp1_offset + this->_nas1 - 1];
                this->_v_grid[this->_nas1 - 1] = this->_u_grid[this->_nas1 - 1];
            }
        

        /* Update u (upward sweep) */
        for (int j = 1; j < this->_nas1 - 1; j++)
        {
            const int s1_idx    = j * this->_nas2;
            const int s1_idx_m1 = s1_idx - this->_nas2;
            const int s1_idx_p1 = s1_idx + this->_nas2;
            
            const T s1 = static_cast<T>(j) * s1_inc;
            for (int k = 1; k < this->_nas2 - 1; k++)
            {
                const T alpha_s1    = lambda * (0.5 * sigma_s1 * sigma_s1 * s1 * s1);
                const T beta_s1     = lambda * r_s1 * s1;
                const T gamma_s1    = 1.0 / (1.0 + alpha_s1 + (r_s1 * t_inc));
            
                const T s2          = static_cast<T>(k) * s2_inc;
                const T alpha_s2    = lambda * (0.5 * sigma_s2 * sigma_s2 * s2 * s2);
                const T beta_s2     = lambda * r_s2 * s2;
                const T gamma_s2    = 1.0 / (1.0 + alpha_s2 + (r_s2 * t_inc));
            
                this->_u_grid[k] = ((fd_grid[tp1_offset + s1_idx + k   ] * (1.0 - alpha_s1)      +
                                     fd_grid[tp1_offset + s1_idx_p1 + k] * (alpha_s1 + beta_s1)  +
                                     this->_u_grid[s1_idx_m1 + k]        * (alpha_s1 - beta_s1)) * gamma_s1) +
                                     
                                   ((fd_grid[tp1_offset + s1_idx + k    ] * (1.0 - alpha_s2)      +
                                     fd_grid[tp1_offset + s1_idx + k + 1] * (alpha_s2 + beta_s2)  +
                                     this->_u_grid[s1_idx + k - 1]        * (alpha_s2 - beta_s2)) * gamma_s2);
            }
        }

        /* Update v (downward sweep) */
        for (int j = this->_nas1 - 2; j > 0; j--)
        {
            const int s1_idx    = j * this->_nas2;
            const int s1_idx_m1 = s1_idx - this->_nas2;
            const int s1_idx_p1 = s1_idx + this->_nas2;
            
            const T s1 = static_cast<T>(j) * s1_inc;
            for (int k = this->_nas2 - 2; k > 0; k--)
            {
                const T alpha_s1    = lambda * (0.5 * sigma_s1 * sigma_s1 * s1 * s1);
                const T beta_s1     = lambda * r_s1 * s1;
                const T gamma_s1    = 1.0 / (1.0 + alpha_s1 + (r_s1 * t_inc));
            
                const T s2          = static_cast<T>(k) * s2_inc;
                const T alpha_s2    = lambda * (0.5 * sigma_s2 * sigma_s2 * s2 * s2);
                const T beta_s2     = lambda * r_s2 * s2;
                const T gamma_s2    = 1.0 / (1.0 + alpha_s2 + (r_s2 * t_inc));

                this->_v_grid[k] = ((fd_grid[tp1_offset + s1_idx + k   ] * (1.0 - alpha_s1)      +
                                     fd_grid[tp1_offset + s1_idx_m1 + k] * (alpha_s1 - beta_s1)  +
                                     this->_v_grid[s1_idx_p1 + k]        * (alpha_s1 + beta_s1)) * gamma_s1) +
                                   
                                   ((fd_grid[tp1_offset + s1_idx + k    ] * (1.0 - alpha_s2)      +
                                     fd_grid[tp1_offset + s1_idx + k - 1] * (alpha_s2 - beta_s2)  +
                                     this->_v_grid[s1_idx + k + 1]        * (alpha_s2 + beta_s2)) * gamma_s2);
            }
        }
        
        /* Update fd grid (average) */
        for (int j = 0; j < this->_nas1; j++)
        {
            const int s1_idx = j * this->_nas2;
            for (int k = 0; k < this->_nas2; k++)
            {
                fd_grid[t_offset + s1_idx + k] = 0.5 * (this->_v_grid[s1_idx + k] + this->_u_grid[s1_idx + k]);
            }
        }
            
        /* Dump this time step */
        if (!this->_log_file.empty() && (!this->_log_final_only || (i == (this->_nts - 1))))
        {
            dump_2d_data_to_gnuplot(this->_log_file + boost::lexical_cast<std::string, int>(i), &fd_grid[t_offset], this->_nas1, this->_nas2, s1_inc, s2_inc);
        }

        std::swap(t_offset, tp1_offset);
    }
    
    /* Dump the grid */
    //dump_2d_data_to_gnuplot(std::cout, &fd_grid[tp1_offset], this->_nas1, this->_nas2, s1_inc, s2_inc);
    const T pos_s1 = this->_economics[0]->get_spot() / s1_max;
    const T pos_s2 = this->_economics[1]->get_spot() / s2_max;
    const T result = fd_grid[static_cast<int>((pos_s1 * static_cast<T>(this->_nas1 * this->_nas2)) + (pos_s2 * static_cast<T>(this->_nas2)))];
    
    /* Clean up */
    delete [] fd_grid;
    
    return result;
}

#endif
