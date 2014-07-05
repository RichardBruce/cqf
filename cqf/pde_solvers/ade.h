#ifndef __ADE_H__
#define __ADE_H__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>

#include <boost/lexical_cast.hpp>
#include <boost/ref.hpp>
#include <boost/thread/thread.hpp>

#include "tbb/task.h"

#include "pde_solver.h"
#include "correctors.h"
#include "matrix.h"
#include "dumping.h"

#include "cashflow.h"


template<class T>
class ade_pde_solver : public pde_solver<T>
{
    public :
        ade_pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<grid<T>*> &grids,
            const boundary_condition<T> &bc, 
            const std::string &log_file, const T max_dt,
            bool log_final_only)
            : pde_solver<T>(cashflows, economics, grids, bc, log_file, max_dt, log_final_only),
            _u_grid(new T[grids[0]->size()]), _v_grid(new T[grids[0]->size()]) { };
         
        ~ade_pde_solver()
        {
            delete [] _u_grid;
            delete [] _v_grid;
        };
        
        T solve() const;

    private :
        T *_u_grid;
        T *_v_grid;
        
        class up_pass_task : public tbb::task
        {
            public :
                up_pass_task(const T *fd_grid, T *u_grid, const T s1_inc, const T s1_inc_inv, const T lambda, const T t_inc, const T sigma_s1, const T r_s1, const int nas) 
                    : _fd_grid(fd_grid), _u_grid(u_grid), _s1_inc(s1_inc), _s1_inc_inv(s1_inc_inv), _lambda(lambda), _t_inc(t_inc), _sigma_s1(sigma_s1), _r_s1(r_s1), _nas(nas)
                {  };
                
                void operator()()
                {
                    execute();
                }

                task* execute()
                {
                    /* Lower Boundary Condition*/
                    _u_grid[0] = exp(-_r_s1 * _t_inc) * _fd_grid[0];
        
                    /* Upper Boundary Condition */
                    {
                        const T s       = static_cast<T>(_nas) * _s1_inc;
                        const T delta   = (_fd_grid[_nas - 1] - _fd_grid[_nas - 2]) * _s1_inc_inv;
                        const T rv      = _r_s1 * _fd_grid[_nas - 1];
                        const T theta   = (_r_s1 * s * delta) - rv;
            
                        _u_grid[_nas - 1] = (theta * _t_inc) + _fd_grid[_nas - 1];
                    }
        

                    /* Update u (upward sweep) */
                    for (int j = 1; j < _nas - 1; j++)
                    {
                        const T s1      = static_cast<T>(j) * _s1_inc;
                        const T alpha   = _lambda * (0.5 * _sigma_s1 * _sigma_s1 * s1 * s1);
                        const T beta    = _lambda * _r_s1 * s1;
                        const T gamma   = 1.0 / (1.0 + alpha + (_r_s1 * _t_inc));
            
                        _u_grid[j]  = (_fd_grid[j    ] * (1.0 - alpha)  +
                                       _fd_grid[j + 1] * (alpha + beta) +
                                       _u_grid[j - 1] * (alpha - beta)) * gamma;
                    }
                    
                    return nullptr;
                }
                
            private :
                const T    *_fd_grid;
                      T    *_u_grid;
                const T     _s1_inc;
                const T     _s1_inc_inv;
                const T     _lambda;
                const T     _t_inc;
                const T     _sigma_s1;
                const T     _r_s1;
                const int   _nas;
        };
        
        class down_pass_task : public tbb::task
        {
            public : 
                down_pass_task(const T *fd_grid, T *v_grid, const T s1_inc, const T s1_inc_inv, const T lambda, const T t_inc, const T sigma_s1, const T r_s1, const int nas) 
                    : _fd_grid(fd_grid), _v_grid(v_grid), _s1_inc(s1_inc), _s1_inc_inv(s1_inc_inv), _lambda(lambda), _t_inc(t_inc), _sigma_s1(sigma_s1), _r_s1(r_s1), _nas(nas)
                {  };
                
                void operator()()
                {
                    execute();
                }

                task* execute()
                {
                    /* Lower Boundary Condition*/
                    this->_v_grid[0] = exp(-_r_s1 * _t_inc) * _fd_grid[0];
        
                    /* Upper Boundary Condition */
                    {
                        const T s       = static_cast<T>(_nas) * _s1_inc;
                        const T delta   = (_fd_grid[_nas - 1] - _fd_grid[_nas - 2]) * _s1_inc_inv;
                        const T rv      = _r_s1 * _fd_grid[_nas - 1];
                        const T theta   = (_r_s1 * s * delta) - rv;
            
                        _v_grid[_nas - 1] = (theta * _t_inc) + _fd_grid[_nas - 1];
                    }


                    /* Update v (downward sweep) */
                    for (int j = _nas - 2; j > 0; j--)
                    {
                        const T s1      = static_cast<T>(j) * _s1_inc;
                        const T alpha   = _lambda * (0.5 * _sigma_s1 * _sigma_s1 * s1 * s1);
                        const T beta    = _lambda * _r_s1 * s1;
                        const T gamma   = 1.0 / (1.0 + alpha + (_r_s1 * _t_inc));
            
                        _v_grid[j] = (_fd_grid[j    ] * (1.0 - alpha)  +
                                      _fd_grid[j - 1] * (alpha - beta) +
                                      _v_grid[j + 1] * (alpha + beta)) * gamma;
                    }
                    
                    return nullptr;
                }
                
            private :
                const T    *_fd_grid;
                      T    *_v_grid;
                const T     _s1_inc;
                const T     _s1_inc_inv;
                const T     _lambda;
                const T     _t_inc;
                const T     _sigma_s1;
                const T     _r_s1;
                const int   _nas;
        };
        
        class update_task : public tbb::task
        {
            public :
                update_task(const T *fd_grid, T *u_grid, T *v_grid, const T s1_inc, const T s1_inc_inv, const T lambda, const T t_inc, const T sigma_s1, const T r_s1, const int nas1) 
                    : _fd_grid(fd_grid), _u_grid(u_grid), _v_grid(v_grid), _s1_inc(s1_inc), _s1_inc_inv(s1_inc_inv), _lambda(lambda), _t_inc(t_inc), _sigma_s1(sigma_s1), _r_s1(r_s1), _nas1(nas1)
                {  };

                task* execute()
                {
                    up_pass_task *a     = new(tbb::task::allocate_child()) up_pass_task(_fd_grid, _u_grid, _s1_inc, _s1_inc_inv, _lambda, _t_inc, _sigma_s1, _r_s1, _nas1);
                    down_pass_task *b   = new(tbb::task::allocate_child()) down_pass_task(_fd_grid, _v_grid, _s1_inc, _s1_inc_inv, _lambda, _t_inc, _sigma_s1, _r_s1, _nas1);
                    // Set ref_count to "two children plus one for the wait".
                    set_ref_count(3);
                    // Start b running.
                    spawn(*b);
                    // Start a running and wait 
                    spawn_and_wait_for_all(*a);
                    
                    return nullptr;
                }
            
            private :
                const T    *_fd_grid;
                      T    *_u_grid;
                      T    *_v_grid;
                const T     _s1_inc;
                const T     _s1_inc_inv;
                const T     _lambda;
                const T     _t_inc;
                const T     _sigma_s1;
                const T     _r_s1;
                const int   _nas1;
        };
};

    
template<class T>    
T ade_pde_solver<T>::solve() const
{
    /* Economics */
    const T r_s1        = this->_economics[0]->get_interest_rate();
    const T sigma_s1    = this->_economics[0]->get_volatility();
    
    /* Grid set up */
    const grid<T> &grid = (*(this->_grids[0]));
    const int nas       = grid.size();
    T *fd_grid          = new T[2 * nas];
    
    const T s1_max      = 300.0;
    const T s1_inc      = s1_max / static_cast<T>(nas);
    const T s1_inc_inv  = 1.0 / s1_inc;

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
        const T t       = (*this->_t_steps)[i - 1];
        const T t_inc   = (*this->_t_steps)[i] - t;
        const T lambda  = t_inc / (s1_inc * s1_inc);
#if 0
#if 0
        update_task *a = new(tbb::task::allocate_root()) update_task(&fd_grid[tp1_offset], this->_u_grid, this->_v_grid, s1_inc, s1_inc_inv, lambda, t_inc, sigma_s1, r_s1, nas);
        tbb::task::spawn_root_and_wait(*a);
#else

        up_pass_task u_task(&fd_grid[tp1_offset], this->_u_grid, s1_inc, s1_inc_inv, lambda, t_inc, sigma_s1, r_s1, nas);
        down_pass_task d_task(&fd_grid[tp1_offset], this->_v_grid, s1_inc, s1_inc_inv, lambda, t_inc, sigma_s1, r_s1, nas);
        
        boost::thread up(boost::ref(u_task));
        boost::thread dn(boost::ref(d_task));
        
        up.join();
        dn.join();
#endif
#else        
        /* Lower Boundary Condition*/
        this->_u_grid[0] = exp(-r_s1 * t_inc) * fd_grid[tp1_offset];
        this->_v_grid[0] = this->_u_grid[0];
        
        /* Upper Boundary Condition */
        {
            const T s       = static_cast<T>(nas) * s1_inc;
            const T delta   = (fd_grid[tp1_offset + nas - 1] - fd_grid[tp1_offset + nas - 2]) * s1_inc_inv;
            const T rv      = r_s1 * fd_grid[tp1_offset + nas - 1];
            const T theta   = (r_s1 * s * delta) - rv;
            
            this->_u_grid[nas - 1] = (theta * t_inc) + fd_grid[tp1_offset + nas - 1];
            this->_v_grid[nas - 1] = this->_u_grid[nas - 1];
        }
        

        /* Update u (upward sweep) */
        for (int j = 1; j < nas - 1; j++)
        {
            const T s1      = static_cast<T>(j) * s1_inc;
            const T alpha   = lambda * (0.5 * sigma_s1 * sigma_s1 * s1 * s1);
            const T beta    = lambda * r_s1 * s1;
            const T gamma   = 1.0 / (1.0 + alpha + (r_s1 * t_inc));
            
            this->_u_grid[j] = (fd_grid[tp1_offset + j    ] * (1.0 - alpha)  +
                                fd_grid[tp1_offset + j + 1] * (alpha + beta) +
                                this->_u_grid[j - 1]        * (alpha - beta)) * gamma;
        }

        /* Update v (downward sweep) */
        for (int j = nas - 2; j > 0; j--)
        {
            const T s1      = static_cast<T>(j) * s1_inc;
            const T alpha   = lambda * (0.5 * sigma_s1 * sigma_s1 * s1 * s1);
            const T beta    = lambda * r_s1 * s1;
            const T gamma   = 1.0 / (1.0 + alpha + (r_s1 * t_inc));
            
            this->_v_grid[j] = (fd_grid[tp1_offset + j    ] * (1.0 - alpha)  +
                                fd_grid[tp1_offset + j - 1] * (alpha - beta) +
                                this->_v_grid[j + 1]        * (alpha + beta)) * gamma;
        }
#endif

        /* Update fd grid (average) */
        for (int j = 0; j < this->_nas1; j++)
        {
            fd_grid[t_offset + j] = 0.5 * (this->_v_grid[j] + this->_u_grid[j]);
        }
                
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
    
    const T pos_s  = this->_economics[0]->get_spot() / s1_max;
    const T result = fd_grid[tp1_offset + static_cast<int>(pos_s * static_cast<T>(nas))];
    
    /* Clean up */
    delete [] fd_grid;
    
    return result;
}

#endif
