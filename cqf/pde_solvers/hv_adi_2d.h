#ifndef __HV_ADI_2D_H__
#define __HV_ADI_2D_H__

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
class hv_adi_2d_pde_solver : public pde_solver_2d<T>
{
    public :
        hv_adi_2d_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const T corr, const int nts, 
            const int nas1, const int nas2, const T t, bool log_final_only)
            : pde_solver_2d<T>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only) { };
         
        ~hv_adi_2d_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T hv_adi_2d_pde_solver<T>::solve() const
{
    /* Economics */
    const T r_s1        = this->_economics[0]->get_interest_rate();
    const T sigma_s1    = this->_economics[0]->get_volatility();
    
    const T r_s2        = this->_economics[1]->get_interest_rate();
    const T sigma_s2    = this->_economics[1]->get_volatility();

    
    /* Grid set up */
    const T t_inc   = this->_t / static_cast<T>(this->_nts);
    const int nsn   = this->_nas1 * this->_nas2;
    T *fd_grid      = new T[nsn * 3];
    
    const T s1_max      = 300.0;
    const T s1_inc      = s1_max / static_cast<T>(this->_nas1);
    const T s1_inc_inv  = 1.0 / s1_inc;
    
    const T s2_max      = 300.0;
    const T s2_inc      = s2_max / static_cast<T>(this->_nas2);
    const T s2_inc_inv  = 1.0 / s2_inc;
    
//    const T theta       = 1.0 - (0.5 * std::sqrt(2.0));
    const T theta       = 0.5 + ((1.0 / 6.0) * std::sqrt(3.0));

    /* Set the end condition (ie payoff) */
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[nsn], this->_t);
    }

    
    /* Work backwards in time */
    int ytild_yzero     = 0; 
    int ytildj_un_unm1  = nsn;
    int yj              = (nsn << 1);
    
    explicit_predictor<T> pred(theta, t_inc, s1_inc, s2_inc, s1_max, s2_max, sigma_s1, sigma_s2, r_s1, r_s2, this->_corr, this->_nas1, this->_nas2);
    implicit_corrector<T> x_corr(r_s1, theta, s1_inc, t_inc, sigma_s1, this->_nas2, this->_nas1);
    implicit_corrector<T> y_corr(r_s2, theta, s2_inc, t_inc, sigma_s2, 1, this->_nas2);
    explicit_corrector<T> exp_corr(0.5, t_inc, s1_inc, s2_inc, r_s1, r_s2, sigma_s1, sigma_s2, this->_corr, this->_nas2);
    for (int i = this->_nts - 2; i >= 0; i--)
    {
        /* Explicit predictor step */
        pred(&fd_grid[ytildj_un_unm1], &fd_grid[ytild_yzero]);

        
        /* S1 direction corrector */
        for (int j = 0; j < this->_nas2; j++)
        {
            x_corr(&fd_grid[ytildj_un_unm1 + j], &fd_grid[ytild_yzero + j], &fd_grid[yj + j]);
        }

                
        /* S2 direction corrector */
        for (int j = 0; j < this->_nas1; j++)
        {
            y_corr(&fd_grid[ytildj_un_unm1 + (this->_nas2 * j)], &fd_grid[yj + (this->_nas2 * j)], &fd_grid[yj + (this->_nas2 * j)]);
        }
        
        
        /* Second predictor step */
        /* Fill in the row */
        for (int j = 1; j < this->_nas1 - 1; j++)
        {
            exp_corr(&fd_grid[yj], &fd_grid[ytildj_un_unm1], &fd_grid[ytild_yzero], j);
        }
    
        /* Upper and lower boundary condition in s1 */
        const int s1_idx_upper  = (this->_nas1 - 1) * this->_nas2;
        const int s1_idx_m1     = s1_idx_upper - this->_nas2;
        const int s1_idx_m2     = s1_idx_m1 - this->_nas2;
        
        const int s1_idx_lower  = 0;
        const int s1_idx_p1     = this->_nas2;
        const int s1_idx_p2     = (this->_nas2 << 1);
        for (int j = 1; j < this->_nas2 - 1; j++)
        {
            fd_grid[ytild_yzero + s1_idx_upper + j] = (2.0 * fd_grid[ytild_yzero + s1_idx_m1 + j]) - fd_grid[ytild_yzero + s1_idx_m2 + j];
            fd_grid[ytild_yzero + s1_idx_lower + j] = (2.0 * fd_grid[ytild_yzero + s1_idx_p1 + j]) - fd_grid[ytild_yzero + s1_idx_p2 + j];
        }
        
        /* The corners */
        /* Apply the s1 = 0 , s2 = 0 boundary condition */
        //fd_grid[ytild_yzero] = (2.0 * fd_grid[ytild_yzero + 1]) - fd_grid[ytild_yzero + 2];
        {
            /* S1 greeks */
            const T s1         = 0.0;
            const T delta_s1   = ((fd_grid[yj + this->_nas2] - fd_grid[yj]) -
                                       (fd_grid[ytildj_un_unm1 + this->_nas2] - fd_grid[ytildj_un_unm1])) * s1_inc_inv;
            const T rv_s1      = (fd_grid[yj] - fd_grid[ytildj_un_unm1]) * r_s1;

            /* s2 greeks */
            const T s2         = 0.0;
            const T delta_s2   = ((fd_grid[yj + 1] - fd_grid[yj]) -
                                       (fd_grid[ytildj_un_unm1 + 1] - fd_grid[ytildj_un_unm1])) * s2_inc_inv;
            const T rv_s2      = (fd_grid[yj] - fd_grid[ytildj_un_unm1]) * r_s2;
            
            /* Update */
            const T corre      = (r_s1 * s1 * delta_s1) - rv_s1 +
                                      (r_s2 * s2 * delta_s2) - rv_s2;
            fd_grid[ytild_yzero]       = (corre * t_inc * theta) + fd_grid[ytild_yzero];
        }
        
        /* Apply the s1 = s1_max, s2 = 0.0 boundary condition */
        //fd_grid[ytild_yzero + s1_idx_upper] = (2.0 * fd_grid[ytild_yzero + s1_idx_upper + 1]) - fd_grid[ytild_yzero + s1_idx_upper + 2];
        {
            /* S1 greeks */
            const T s1         = s1_max;
            const T delta_s1   = ((fd_grid[yj + s1_idx_upper] - fd_grid[yj + s1_idx_upper - this->_nas2]) -
                                       (fd_grid[ytildj_un_unm1 + s1_idx_upper] - fd_grid[ytildj_un_unm1 + s1_idx_upper - this->_nas2])) * s1_inc_inv;
            const T rv_s1      = (fd_grid[yj + s1_idx_upper] - fd_grid[ytildj_un_unm1 + s1_idx_upper]) * r_s1;

            /* s2 greeks */
            const T s2         = 0.0;
            const T delta_s2   = ((fd_grid[yj + s1_idx_upper + 1] - fd_grid[yj + s1_idx_upper]) -
                                       (fd_grid[ytildj_un_unm1 + s1_idx_upper + 1] - fd_grid[ytildj_un_unm1 + s1_idx_upper])) * s2_inc_inv;
            const T rv_s2      = (fd_grid[yj + s1_idx_upper] - fd_grid[ytildj_un_unm1 + s1_idx_upper]) * r_s2;
            
            /* Update */
            const T corre      = (r_s1 * s1 * delta_s1) - rv_s1 +
                                      (r_s2 * s2 * delta_s2) - rv_s2;
            fd_grid[ytild_yzero + s1_idx_upper]    = (corre * t_inc * theta) + fd_grid[ytild_yzero + s1_idx_upper];
        }
        
        /* Apply the s1 = 0.0, s2 = s2_max boundary condition */
        //fd_grid[ytild_yzero + this->_nas2 - 1] = (2.0 * fd_grid[ytild_yzero + this->_nas2 - 2]) - fd_grid[ytild_yzero + this->_nas2 - 3];
        {
            /* S1 greeks */
            const T s1         = 0.0;
            const T delta_s1   = ((fd_grid[yj + (this->_nas2 << 1) - 1] - fd_grid[yj + this->_nas2 - 1]) -
                                       (fd_grid[ytildj_un_unm1 + (this->_nas2 << 1) - 1] - fd_grid[ytildj_un_unm1 + this->_nas2 - 1])) * s1_inc_inv;
            const T rv_s1      = (fd_grid[yj + this->_nas2 - 1] - fd_grid[ytildj_un_unm1 + this->_nas2 - 1]) * r_s1;

            /* s2 greeks */
            const T s2         = s2_max;
            const T delta_s2   = ((fd_grid[yj + this->_nas2 - 1] - fd_grid[yj + this->_nas2 - 2]) -
                                       (fd_grid[ytildj_un_unm1 + this->_nas2 - 1] - fd_grid[ytildj_un_unm1 + this->_nas2 - 2])) * s2_inc_inv;
            const T rv_s2      = (fd_grid[yj + this->_nas2 - 1] - fd_grid[ytildj_un_unm1 + this->_nas2 - 1]) * r_s2;
            
            /* Update */
            const T corre      = (r_s1 * s1 * delta_s1) - rv_s1 +
                                      (r_s2 * s2 * delta_s2) - rv_s2;
            fd_grid[ytild_yzero + this->_nas2 - 1]    = (corre * t_inc * theta) + fd_grid[ytild_yzero + this->_nas2 - 1];
        }

        /* Apply the s1 = s1_max, s2 = s2_max boundary condition */
        //fd_grid[ytild_yzero + nsn - 1] = (2.0 * fd_grid[ytild_yzero + nsn - 2]) - fd_grid[ytild_yzero + nsn - 3];
        {
            /* S1 greeks */
            const T s1         = s1_max;
            const T delta_s1   = ((fd_grid[yj + nsn - 1] - fd_grid[yj + nsn - 1 - this->_nas2]) -
                                       (fd_grid[ytildj_un_unm1 + nsn - 1] - fd_grid[ytildj_un_unm1 + nsn - 1 - this->_nas2])) * s1_inc_inv;
            const T rv_s1      = (fd_grid[yj + nsn - 1] - fd_grid[ytildj_un_unm1 + nsn - 1]) * r_s1;

            /* s2 greeks */
            const T s2         = s2_max;
            const T delta_s2   = ((fd_grid[yj + nsn - 1] - fd_grid[yj + nsn - 2]) -
                                       (fd_grid[ytildj_un_unm1 + nsn - 1] - fd_grid[ytildj_un_unm1 + nsn - 2])) * s2_inc_inv;
            const T rv_s2      = (fd_grid[yj + nsn - 1] - fd_grid[ytildj_un_unm1 + nsn - 1]) * r_s2;
            
            /* Update */
            const T corre      = (r_s1 * s1 * delta_s1) - rv_s1 +
                                      (r_s2 * s2 * delta_s2) - rv_s2;
            fd_grid[ytild_yzero + nsn - 1]       = (corre * t_inc * theta) + fd_grid[ytild_yzero + nsn - 1];
        }


        /* S1 direction corrector */
        for (int j = 0; j < this->_nas2; j++)
        {
            x_corr(&fd_grid[yj + j], &fd_grid[ytild_yzero + j], &fd_grid[ytildj_un_unm1 + j]);
        }

                
        /* S2 direction corrector */
        for (int j = 0; j < this->_nas1; j++)
        {
            y_corr(&fd_grid[yj + (this->_nas2 * j)], &fd_grid[ytildj_un_unm1 + (this->_nas2 * j)], &fd_grid[ytildj_un_unm1 + (this->_nas2 * j)]);
        }
        
        
        //break;
        if (!this->_log_file.empty() && (!this->_log_final_only || (i == (this->_nts - 1))))
        {
            dump_2d_data_to_gnuplot(this->_log_file + boost::lexical_cast<std::string, int>(i), &fd_grid[ytildj_un_unm1], this->_nas1, this->_nas2, s1_inc, s2_inc);
        }
    }


    /* Dump the grid */
    //dump_2d_data_to_gnuplot(std::cout, &fd_grid[ytildj_un_unm1], this->_nas1, this->_nas2, s1_inc, s2_inc);
    const T pos_s1 = this->_economics[0]->get_spot() / s1_max;
    const T pos_s2 = this->_economics[1]->get_spot() / s2_max;
    const T result = fd_grid[static_cast<int>((pos_s1 * static_cast<T>(this->_nas1 * this->_nas2)) + (pos_s2 * static_cast<T>(this->_nas2)))];
    
    /* Clean up */
    delete [] fd_grid;
    
    return result;
}

#endif
