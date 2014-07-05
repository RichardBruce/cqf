#ifndef __DOUGLAS_ADI_2D_H__
#define __DOUGLAS_ADI_2D_H__

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
class douglas_adi_2d_pde_solver : public pde_solver_2d<T>
{
    public :
        douglas_adi_2d_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const T corr, const int nts, 
            const int nas1, const int nas2, const T t, bool log_final_only)
            : pde_solver_2d<T>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only) { };
         
        ~douglas_adi_2d_pde_solver() { };
    
        T solve() const;
        
    
};


template<class T>    
T douglas_adi_2d_pde_solver<T>::solve() const
{    
    /* Economics */
    const T r_s1        = this->_economics[0]->get_interest_rate();
    const T sigma_s1    = this->_economics[0]->get_volatility();
    
    const T r_s2        = this->_economics[1]->get_interest_rate();
    const T sigma_s2    = this->_economics[1]->get_volatility();

    
    /* Grid set up */
    const T t_inc   = this->_t / static_cast<T>(this->_nts);
    const int nsn   = this->_nas1 * this->_nas2;
    T *fd_grid      = new T[nsn * 2];
    
    const T s1_max  = 300.0;
    const T s1_inc  = s1_max / static_cast<T>(this->_nas1);
    
    const T s2_max  = 300.0;
    const T s2_inc  = s2_max / static_cast<T>(this->_nas2);

    const T theta   = 0.5;
        
    /* Set the end condition (ie payoff) */
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[nsn], this->_t);
    }

    
    /* Work backwards in time */
    int t_offset      = 0; 
    int tp1_offset    = nsn;
    
    explicit_predictor<T> pred(theta, t_inc, s1_inc, s2_inc, s1_max, s2_max, sigma_s1, sigma_s2, r_s1, r_s2, this->_corr, this->_nas1, this->_nas2);
    implicit_corrector<T> x_corr(r_s1, theta, s1_inc, t_inc, sigma_s1, this->_nas2, this->_nas1);
    implicit_corrector<T> y_corr(r_s2, theta, s2_inc, t_inc, sigma_s2, 1, this->_nas2);
    for (int i = 0; i < this->_nts; i++)
    {
        /* Explicit predictor step */
        pred(&fd_grid[tp1_offset], &fd_grid[t_offset]);

        //break;


        /* S1 direction corrector */
        for (int j = 0; j < this->_nas2; j++)
        {
            x_corr(&fd_grid[tp1_offset + j], &fd_grid[t_offset + j], &fd_grid[t_offset + j]);
        }
        
        //break;

                
        /* S2 direction corrector */
        for (int j = 0; j < this->_nas1; j++)
        {
            y_corr(&fd_grid[tp1_offset + (this->_nas2 * j)], &fd_grid[t_offset + (this->_nas2 * j)], &fd_grid[tp1_offset + (this->_nas2 * j)]);
        }
        
        
        //break;
        if (!this->_log_file.empty() && (!this->_log_final_only || (i == (this->_nts - 1))))
        {
            dump_2d_data_to_gnuplot(this->_log_file + boost::lexical_cast<std::string, int>(i), &fd_grid[tp1_offset], this->_nas1, this->_nas2, s1_inc, s2_inc);
        }
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
