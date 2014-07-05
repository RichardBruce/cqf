#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include <assert.h>

#include "pde_solver.h"


template<class T>
class explicit_pde_solver : public pde_solver<T>
{
    public :
        explicit_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<grid<T>*> &grids,
            const boundary_condition<T> &bc,
            const std::string &log_file, 
            const T max_dt, bool log_final_only)
            : pde_solver<T>(cashflows, economics, grids, bc, log_file, 
                std::min(max_dt, get_max_dt(economics[0]->get_volatility(), grids[0]->size())), 
                log_final_only)
            {
                assert(this->_grids.size() == 1);
            };

         
        ~explicit_pde_solver() { };
        
        T solve() const;
        
    private :
        T get_max_dt(const T sigma, const int nas)
        {
            return 0.9 / (static_cast<T>(nas * nas) * sigma * sigma);
        }
};

    
template<class T>    
T explicit_pde_solver<T>::solve() const
{
    /* Economics */
    const T r               = this->_economics[0]->get_interest_rate();
    const T sigma           = this->_economics[0]->get_volatility();
    const T half_sigma_sq   = 0.5 * sigma * sigma;
    
    /* Allocate grid for V, note t_incis redefined */
    const grid<T> &grid = (*(this->_grids[0]));
    const int nas       = grid.size();
    T *fd_grid          = new T[2 * nas];  
    
    /* Set the end condition (ie payoff) */
    int tp1_offset  = nas;
    int t_offset    = 0;
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
    

    /* Work backwards in time */
    for (int i = this->_t_steps->size() - 1; i > 0; i--)
    {
        /* Fill in the middle values */
        const T t       = (*this->_t_steps)[i - 1];
        const T t_inc   = (*this->_t_steps)[i] - t;
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
          
            fd_grid[t_offset + j]  = a * t_inc * fd_grid[tp1_offset + j - 1];
            fd_grid[t_offset + j] += (1.0 + b * t_inc) * fd_grid[tp1_offset + j];
            fd_grid[t_offset + j] += c * t_inc * fd_grid[tp1_offset + j + 1];
        }
    
        /* Apply the s = 0 boundary condition */
        fd_grid[t_offset] = fd_grid[tp1_offset] * std::exp(-r * t_inc);
        
        /* Apply the s = s_max boundary condition */
        fd_grid[t_offset + nas - 1] = (2.0 * fd_grid[t_offset + nas - 2]) - fd_grid[t_offset + nas - 3];
        
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
    delete [] fd_grid;
    
    return result;
}
