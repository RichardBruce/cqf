#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>

#include "pde_solver.h"


template<class T>
class explicit_pde_solver : public pde_solver<T>
{
    public :
        explicit_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const int nts, const int nas1, 
            const T t, bool log_final_only)
            : pde_solver<T>(cashflows, economics, log_file, nts, nas1, t, log_final_only) { };
         
        ~explicit_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T explicit_pde_solver<T>::solve() const
{    
    /* Economics */
    const T r               = this->_economics[0]->get_interest_rate();
    const T sigma           = this->_economics[0]->get_volatility();
    const T half_sigma_sq   = 0.5 * sigma * sigma;
    
    
    /* Grid set up, not nts is redefined */
    T t_inc         = 0.9 / (static_cast<T>(this->_nas1 * this->_nas1) * sigma * sigma);
    const int nts   = (int)(this->_t / t_inc) + 1;
    t_inc           = this->_t / this->_nts;
    const int ngn   = this->_nas1 * this->_nts;
    T *fd_grid      = new T[ngn];
    
    const T s_max           = 300.0;
    const T s_inc           = s_max / static_cast<T>(this->_nas1);
    const T s_inc_inv       = 1.0 / s_inc;
    const T s_inc2_inv      = 0.5 * s_inc_inv;
    const T s_inc_sq_inv    = 1.0 / (s_inc * s_inc);
    
    
    /* Set the end condition (ie payoff) */
    const int payoff_offset = ngn - this->_nas1;
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[payoff_offset], this->_t);
    }
    

    /* Work backwards in time */
    for (int i = this->_nts - 2; i >= 0; i--)
    {
        const int t_offset      = i * this->_nas1; 
        const int tp1_offset    = t_offset + this->_nas1;
        
        /* Fill in the middle values */
        for (int j = 1; j < this->_nas1 - 1; j++)
        {
            const T s       = static_cast<T>(j) * s_inc;
            const T delta   = (fd_grid[tp1_offset + j + 1] - fd_grid[tp1_offset + j - 1]) * s_inc2_inv;
            const T gamma   = (fd_grid[tp1_offset + j + 1] - (2.0 * fd_grid[tp1_offset + j]) + fd_grid[tp1_offset + j - 1]) * s_inc_sq_inv;
            const T rv      = r * fd_grid[tp1_offset + j];
            const T theta   = (half_sigma_sq * s * s * gamma) + (r * s * delta) - rv;
            fd_grid[t_offset + j]   = (theta * t_inc) + fd_grid[tp1_offset + j];
        }
    
        /* Apply the s = 0 boundary condition */
        fd_grid[i * this->_nas1] = fd_grid[(i * this->_nas1) + this->_nas1] * std::exp(-r * t_inc);
        
        /* Apply the s = s_max boundary condition */
        fd_grid[tp1_offset - 1] = (2.0 * fd_grid[tp1_offset - 2]) - fd_grid[tp1_offset - 3];
    }
    
    dump_2d_data_to_gnuplot(std::cout, &fd_grid[0], this->_nts + 1, this->_nas1, t_inc, s1_inc);
    const T result = fd_grid[static_cast<int>((1.0/3.0) * static_cast<T>(this->_nas1))];
    
    /* Clean up */
    delete [] fd_grid;
    
    return result;
}
