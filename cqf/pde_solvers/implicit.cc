#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>

#include "matrix.h"
#include "pde_solver.h"


template<class T>
class implicit_pde_solver : public pde_solver<T>
{
    public :
        implicit_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const int nts, const int nas1, 
            const T t, bool log_final_only)
            : pde_solver<T>(cashflows, economics, log_file, nts, nas1, t, log_final_only) { };
         
        ~implicit_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T implicit_pde_solver<T>::solve() const
{
    /* Parse inputs */
    bool use_sor    = false;
    int sor_max_it  = 5000;
    T sor_w         = 1.85;
    T sor_tol       = 1e-6;
   
    
    /* Economics */
    const T r       = this->_economics[0]->get_interest_rate();
    const T sigma   = this->_economics[0]->get_volatility();
    
    
    /* Grid set up */
    const int ngn   = this->_nas1 * this->_nts;
    T *fd_grid = new T[ngn];
    
    const T s_max  = 300.0;
    const T s_inc  = s_max / static_cast<T>(this->_nas1);
    const T t_inc  = t     / static_cast<T>(this->_nts);
    
    /* Pre-computes */
    const T sigma_sq       = sigma * sigma;
    const T neg_r_t_inc    = -r * t_inc;

    
    /* Set the end condition (ie payoff) */
    const int payoff_offset = ngn - this->_nas1;
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[payoff_offset], this->_t);
    }
    
    for (int i = this->_nts - 2; i >= 0; i--)
    {
        /* Populate the matrix */
        T *matrix_data = new T [this->_nas1 * this->_nas1];
        for (int j = 1; j < this->_nas1 - 1; j++)
        {
            const int matrix_offset = (j * this->_nas1) + j;
            const T flt_j = static_cast<T>(j);
            matrix_data[matrix_offset - 1] = 0.5 * flt_j * (r - sigma_sq * flt_j) * t_inc;
            matrix_data[matrix_offset    ] = 1.0 + (r + sigma_sq * flt_j * flt_j) * t_inc;
            matrix_data[matrix_offset + 1] = 0.5 * flt_j * (-r - sigma_sq * flt_j) * t_inc;
        }
        
        matrix_data[0] = std::exp(neg_r_t_inc);
        matrix_data[(this->_nas1 * this->_nas1) - 1] = 2.0;
        matrix_data[(this->_nas1 * this->_nas1) - 2] = -1.0;
        
        /* Solve the matrix */
        matrix<T> mat(matrix_data, this->_nas1, this->_nas1, true);
        if (use_sor)
        {
            mat.sor_solve(&fd_grid[(i * this->_nas1) + this->_nas1], &fd_grid[i * this->_nas1], sor_w, sor_tol, sor_max_it, Tridiagonal);
        }
        else
        {
            mat.gauss_solve(&fd_grid[(i * this->_nas1) + this->_nas1], &fd_grid[i * this->_nas1], Tridiagonal);
        }
    }
    
    
    /* Dump the grid */
    for (int i = 0; i < this->_nas1; i++)
    {
        std::cout << (i * s_inc) << ", ";
        for (int j = 0; j < nts; j++)
        {
            std::cout << fd_grid[(j * this->_nas1) + i] << ", ";
        }
        std::cout << std::endl;
    }
    
    
    /* Clean up */
    delete [] fd_grid;
    
    return 0;
}
