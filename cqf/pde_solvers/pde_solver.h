#ifndef __PDE_SOLVER_H__
#define __PDE_SOLVER_H__

#include <algorithm>
#include <map>
#include <vector>

#include "boost/noncopyable.hpp"

#include "boundary_condition.h"

#include "utility.h"
#include "scenario_point.h"

template<class T>
class cashflow;

template<class T>
class equity_economics;

template<class T>
class risk;

//template<class T>
//class boundary_condition;

template<class T>
class pde_solver : private boost::noncopyable
{
    public :
        /* CTOR */
        /*pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<grid<T>*> &grids,
            const boundary_condition<T> &bc,
            const std::string &log_file, const int nts, 
            const int nas1, const T t, bool log_final_only)
            : _cashflows(cashflows), _economics(economics), _grids(grids), 
              _t_steps(nullptr), _bc(bc), 
              _log_file(log_file), _t(t), _nts(nts), _nas1(nas1), 
              _log_final_only(log_final_only) { };*/
              
        pde_solver(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::vector<risk<T>*> &risks,
            const boundary_condition<T> &bc,
            const std::string &log_file, 
            const T max_dt, bool log_final_only)
            : _cashflows(cashflows), _economics(economics), _risks(risks), 
              _t_steps(build_timesteps(max_dt)),
              _bc(bc), 
              _log_file(log_file),
              _log_final_only(log_final_only)
            { };
        
        /* Virtual DTOR */    
        virtual ~pde_solver()
        {
            if (_t_steps != nullptr)
            {
                delete _t_steps;
            }
        };
        
        void solve()
        {
            /* Build scenario points */
            scenario_point<T> null_bump(this->_t_steps->back(), this->_economics[0]->get_spot(), 
                this->_economics[0]->get_interest_rate(), this->_economics[0]->get_volatility());
            std::map<scenario_point *, pingpong_buffer *> scenarios;
            for (int i = 0; i < this->_risks.size(); i++)
            {
                this->_risks[i]->bump(scenarios, null_bump);
            }
            
            /* Set the end condition (ie payoff) */
            bool disc_cashflow = false;
            for (scenario_map::iterator i = scenarios.begin(); i != scenarios.end(); i++)
            {
                T* initial_cond = (*i).second->write_buffer();
                for (unsigned int j = 0; j < this->_cashflows.size(); j++)
                {
                    if (this->_cashflows[j]->exists(this->_t_steps->back()))
                    {
                        disc_cashflow |= this->_cashflows[j]->add(initial_cond, this->_grids, this->_t_steps->back());
                    }
                }
                (*i).second->flip();
            }
            
            /* Work backwards in time */
            for (int i = this->_t_steps->size() - 1; i > 0; i--)
            {
                const T t       = (*this->_t_steps)[i - 1];
                const T t_inc   = (*this->_t_steps)[i] - t;
                
                /* Time step all the scenario points */
                for (scenario_map::iterator j = scenarios.begin(); j != scenarios.end(); j++)
                {
                    this->iteration((*i).second->write_buffer(), (*i).second->read_buffer(), (*i).first, t_inc);

                    /* Process cashflows */
                    disc_cashflow = false;
                    for (unsigned int k = 0; k < this->_cashflows.size(); k++)
                    {
                        if (this->_cashflows[k]->exists(t))
                        {
                            disc_cashflow |= this->_cashflows[k]->add((*j).second->write_buffer(), this->_grids, t);
                        }
                    }
                
                    /* Swap read and write locations */
                    (*j).second->flip();
                }
                
                /* Aggregate all the risk */
                for (int j = 0; j < this->_risks.size(); j++)
                {
                    this->_risks[j]->aggregate(scenarios);
                }
            }
        }
        
    
    protected :
        const std::vector<cashflow<T> *>           &_cashflows;
        const std::vector<equity_economics<T>*>    &_economics;
        const std::vector<risk<T>*>                &_risks;
        std::vector<T>                             *_t_steps;
        const boundary_condition<T>                &_bc;
        const std::string                          &_log_file;
        const bool                                  _log_final_only;
        
        /* Template method to plug in one iteration of some PDE solver algorithm */        
        virtual void iteration(T *const t_offset, const T *const tp1_offset, const T t_inc, bool disc_cashflow) const = 0;

    private :
        std::vector<T> * build_timesteps(const T max_dt)
        {
            /* Build time steps */
            std::vector<T> required_dates;
            required_dates.push_back(0.0);
            std::back_insert_iterator<std::vector<T> > dates_inserter(required_dates);
            for (unsigned int i = 0; i < _cashflows.size(); i++)
            {
                _cashflows[i]->required_date(dates_inserter);
            }
    
            /* Clean up any duplicates */
            std::sort(required_dates.begin(), required_dates.end());
            typename std::vector<T>::iterator dates_iter = std::unique(required_dates.begin(), required_dates.end());
            required_dates.resize(dates_iter - required_dates.begin());
    
            /* Pad the dates for the minimum step size */
            std::vector<T> *padded_dates = new std::vector<T>();
            pad_data<T>(padded_dates, required_dates, max_dt);
            
            return padded_dates;
        }
};

template<class T>
class pde_solver_2d : public pde_solver<T>
{
    public :
        pde_solver_2d(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const T corr, const int nts, 
            const int nas1, const int nas2, const T t, bool log_final_only)
            : pde_solver<T>(cashflows, economics, std::vector<grid<T>*>(), const_ds_boundary_condition<T>(0.0, 300.0), log_file, nts, nas1, t, log_final_only), 
            _corr(corr), _nas2(nas2)
            {
                assert(this->_economics.size() == 2);
            };
            
        virtual ~pde_solver_2d() {  };
        
    protected :
        const T     _corr;
        const int   _nas2;
};

#endif
