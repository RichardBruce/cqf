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
            const std::string &log_file, const int nts, const int nas1, 
            const T t, bool log_final_only)
            : pde_solver<T>(cashflows, economics, log_file, nts, nas1, t, log_final_only) { };
         
        ~douglas_pde_solver() { };
        
        T solve() const;
};

    
template<class T>    
T douglas_pde_solver<T>::solve() const
{    
    /* Economics */
    const T r       = this->_economics[0]->get_interest_rate();
    const T sigma   = this->_economics[0]->get_volatility();
    
    /* Grid set up */
    const int ngn = this->_nas1 * this->_nts;
    
    const T s_max = 300.0;
    const T s_inc = s_max / static_cast<T>(this->_nas1);
    
    T *fd_grid = new T[ngn];
    
    /* Pre-computes */
    const T s_inc_sq           = s_inc * s_inc;
    const T half_sigma_sq      = 0.5 * sigma * sigma;
    const T t_inc              = this->_t / static_cast<T>(this->_nts);
    const T t_inc_s_inc_sq     = t_inc / s_inc_sq;
    const T t_inc_s_inc        = t_inc / s_inc;    
    
    /* Set the end condition (ie payoff) */
    const int payoff_offset = ngn - this->_nas1;
    for (unsigned int i = 0; i < this->_cashflows.size(); i++)
    {
        this->_cashflows[i]->add(&fd_grid[payoff_offset], this->_t);
    }
    
    /* Work backwards in time */
    const T theta      = 0.5 + (s_inc_sq / (12.0 * t_inc));
    const T m_theta    = 1.0 - theta;
    std::cout << "Theta is " << theta << std::endl;
    assert(theta > 0.0);
    assert(theta < 1.0);
    for (int i = this->_nts - 2; i >= 0; i--)
    {
        T *matrix_data_k    = new T [this->_nas1 * this->_nas1];
        T *matrix_data_kp1  = new T [this->_nas1];
        memset(matrix_data_k, 0, (this->_nas1 * this->_nas1 * sizeof(T)));
    
        /* Populate the matrix */
        for (int j = 1; j < this->_nas1 - 1; j++)
        {
            const T s      = static_cast<T>(j) * s_inc;
            const T s_sq   = s * s;
            
            const T g = half_sigma_sq * s_sq;
            const T r_s = r * s;
            
            const T a = (t_inc_s_inc_sq * g) - (0.5 * t_inc_s_inc * r_s);
            const T b = (t_inc * -r) - (2.0 * t_inc_s_inc_sq * g);
            const T c = (t_inc_s_inc_sq * g) + (0.5 * t_inc_s_inc * r_s);
            
            matrix_data_k[(j * this->_nas1) + j - 1] = -a * theta;
            matrix_data_k[(j * this->_nas1) + j    ] = 1.0 - b * theta;
            matrix_data_k[(j * this->_nas1) + j + 1] = -c * theta;
            
            matrix_data_kp1[j]  = a * m_theta * fd_grid[(i * this->_nas1) + this->_nas1 + j - 1];
            matrix_data_kp1[j] += (1.0 + b * m_theta) * fd_grid[(i * this->_nas1) + this->_nas1 + j];
            matrix_data_kp1[j] += c * m_theta * fd_grid[(i * this->_nas1) + this->_nas1 + j + 1];
        }
        
        matrix_data_k[0]                    = exp(-r * t_inc);
        matrix_data_k[(this->_nas1 * this->_nas1) - 1]      = 2.0;
        matrix_data_k[(this->_nas1 * this->_nas1) - 2]      = -1.0;
        matrix_data_kp1[0]                  = exp(-r * t_inc) * fd_grid[(i * this->_nas1) + this->_nas1];
        matrix_data_kp1[this->_nas1 - 1]            = 2.0 * fd_grid[(i * this->_nas1) + this->_nas1 + this->_nas1 - 1] - fd_grid[(i * this->_nas1 ) + this->_nas1 + this->_nas1 - 2];
        
        matrix<T> matrix_kp1(&matrix_data_kp1[0], this->_nas1, 1, true);
        matrix<T> matrix_k(&matrix_data_k[0], this->_nas1, this->_nas1, true);
        matrix_k.gauss_solve(matrix_kp1, &fd_grid[i * this->_nas1], Tridiagonal);
    }


    /* Get result */
    if (!this->_log_file.empty())
    {
        dump_2d_data_to_gnuplot(std::cout, &fd_grid[0], this->_nts + 1, this->_nas1, t_inc, s_inc);
    }
    const T pos_s  = this->_economics[0]->get_spot() / s_max;
    const T result = fd_grid[static_cast<int>(pos_s * static_cast<T>(this->_nas1))];
    
    /* Clean up */
    delete [] fd_grid;
    
    return result;
}
