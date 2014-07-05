#ifndef __CASHFLOW_H__
#define __CASHFLOW_H__

#include "boost/noncopyable.hpp"
#include "boost/bind.hpp"

#include "barrier_cashflow.h"

#include "extrapolation.h"
#include "interpolation.h"


template<class T>
class grid;


/* Abstract base for the application of a cash flow that may be a fucntion of s and t */
/* Final conditions are also cash flows (of the size of the derivative payoff */
/* Should return true if the cashflow has a dicontinuity in s */
template<class T>
class cashflow : private boost::noncopyable
{
    public :
        cashflow(const T t, const T s0_in, const bool always_eval = false) 
            : _t(1, t), _s0_in(s0_in), _always_eval(always_eval) {  };
        
        cashflow(const std::vector<T> t, const T s0_in, const bool always_eval = false)
            : _t(t), _s0_in(s0_in), _always_eval(always_eval) {  };
            
        /* Virtual DTOR, this class will be used as a base class */
        virtual ~cashflow() {  };
        
        
        /* Common functions */
        bool required_at(const T t) const
        {
            return _always_eval || (std::find(_t.begin(), _t.end(), t) != _t.end());
        }
        
        void required_date(std::back_insert_iterator<std::vector<T>> iter) const
        {
            std::copy(_t.begin(), _t.end(), iter);
        }
        
        /* Function to discover the dimensionality of the cashflows */
        virtual void dimensions(std::vector<bool> *const d) const = 0;
        
        /* Virutal function to ensure strikes are on the PDE grid */
        virtual void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const = 0;
        
        /* Virtual function to allow application of boundary conditions */
        virtual boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const = 0;
        
        /* Virtual function for processing a cashflow */
        virtual bool add(grid<T> &g, const T t) const = 0;
        
    protected :
        unsigned int time_index_of(const T t) const
        {
            return std::distance(_t.begin(), std::find(_t.begin(), _t.end(), t));
        }
        
        unsigned int number_of_eval_dates() const
        {
            return _t.size();
        }
        
        T get_s0_in() const { return _s0_in; }
    
    private :
        const std::vector<T>    _t;
        const T                 _s0_in;
        const bool              _always_eval;
};


/* Single asset based cashflows */
template<class T, class S>
class one_asset_payoff_cashflow : public cashflow<T>
{
    public :
        one_asset_payoff_cashflow(const S func, const T t, const T s0_in)
            : cashflow<T>(t, s0_in, func.always_eval), _func(func)
            {  };

        void dimensions(std::vector<bool> *const d) const
        {
            /* Request 1 connected dimension */
            assert(d->empty() || !(*d)[0] || !"Error: Required dimension is already registered as disconnected.");
            if (d->empty())
            {
                d->push_back(false);
            }
            
            return;
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            _func.required_stock_price(iter);
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _func.boundary_conditions(grid, value);
        }
        
        bool add(grid<T> &g, const T t) const
        {
            assert((g.dimensions() == 1) || ((g.dimensions() == 2) && g.outer_grid_disconnected()));
        
            /* Get the grid sizes */
            const int dim_m1        = g.dimensions() - 1;
            const int nr_of_grids   = (dim_m1 == 0) ? 1 : g.size(0);
            
            /* For each disconnected grid process the cashflow */
            int payoff_offset       = 0;
            const std::vector<T> &s = g.s_grid(dim_m1);
            T* v_grid               = g.write_v_grid();
            for (int i = 0; i < nr_of_grids; i++)
            {
                for (unsigned int j = 0; j < s.size(); j++)
                {
                    v_grid[payoff_offset] = _func(v_grid[payoff_offset], get_s0_in(), s[j]);
                    ++payoff_offset;
                }
            }
            
            return _func.discontinious;
        }
                
    private :
        using cashflow<T>::get_s0_in;
        
        const S _func;
};


/* Single asset based cashflows */
template<class T, class S>
class one_asset_asian_payoff_cashflow : public cashflow<T>
{
    public :
        one_asset_asian_payoff_cashflow(const S func, const std::vector<T> t, const T s0_in, const T asian_in)
            : cashflow<T>(t, s0_in, func.always_eval), _func(func), 
            _cross_asian_v(new std::vector<T>()), _asian_tmp(new std::vector<T>()), 
            _tmp(new std::vector<T>()), _asian_in(asian_in)
            {  };
            
        ~one_asset_asian_payoff_cashflow()
        {
            delete _cross_asian_v;
            delete _asian_tmp;
            delete _tmp;
        }

        void dimensions(std::vector<bool> *const d) const
        {
            /* Request 1 connected dimension and 1 disconnected dimension */
            assert(d->empty() || !(*d)[0] || !"Error: Required dimension is already registered as disconnected.");
            if (d->empty())
            {
                d->push_back(false);
            }
            
            assert((d->size() < 2) || (*d)[1] || !"Error: Required dimension is already registered as connected.");
            if (d->size() < 2)
            {
                d->push_back(true);
            }
            
            return;
        }

        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            _func.required_stock_price(iter);
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _func.boundary_conditions(grid, value);
        }

        bool add(grid<T> &g, const T t) const
        {
            assert((g.dimensions() == 2) && g.outer_grid_disconnected());
        
            /* For each disconnected grid process the cashflow */
            const std::vector<T> &a_grid = g.s_grid(0);
            const std::vector<T> &s_grid = g.s_grid(1);
            T* v_grid = g.write_v_grid();
            
            const unsigned int t_idx    = this->time_index_of(t) + 1;
            const unsigned int t_idx_m1 = t_idx - 1;
            const T flt_t_idx           = static_cast<T>(t_idx);
            const T t_idx_inv           = 1.0 / flt_t_idx;
            const T t_idx_m1_inv        = 1.0 / (flt_t_idx - 1.0);
            const T a_fac               = (flt_t_idx - 1.0) * t_idx_inv;
            if (t_idx == 1)
            {
                /* Asian grid reduced to 1 node */
                _asian_tmp->resize(1);
                (*_asian_tmp)[0] = 0.0;
                
                permute_asian_values(a_grid, s_grid, v_grid, a_fac, t_idx_inv, t);
                g.replace_s_grid((*_asian_tmp), 0);
            }
            else if (t_idx == this->number_of_eval_dates())
            {
                /* Final payoff */
                int payoff_offset = 0;
                if (_asian_in)
                {
                    for (unsigned int i = 0; i < a_grid.size(); i++)
                    {
                        for (unsigned int j = 0; j < s_grid.size(); j++)
                        {
                            v_grid[payoff_offset] = _func(v_grid[payoff_offset], a_grid[i], s_grid[j]);
                            ++payoff_offset;
                        }
                    }
                }
                else
                {
                    for (unsigned int i = 0; i < a_grid.size(); i++)
                    {
                        for (unsigned int j = 0; j < s_grid.size(); j++)
                        {
                            const T avg = (a_grid[i] * a_fac) + (s_grid[j] * t_idx_inv);
                            v_grid[payoff_offset] = _func(v_grid[payoff_offset], get_s0_in(), avg);
                            ++payoff_offset;
                        }
                    }
                }
            }
            else
            {
                /* Asianing update */
                /* Rebuild uniform asian grid */
                /* This may be slightly too big, but not too much to worry about */
                _asian_tmp->resize((a_grid.size() * a_fac) + 1);
                
                /* This loop will get all, but the last t_idx - 1 or less points */
                const unsigned int whole_iters = (a_grid.size() - 1) / t_idx;
                for (unsigned int i = 0; i < whole_iters; i++)
                {
                    const unsigned int l_idx = i * t_idx;
                    const unsigned int h_idx = l_idx + t_idx;
                    const T step = (a_grid[h_idx] - a_grid[l_idx]) * t_idx_m1_inv;
                    for (unsigned int j = 0; j < t_idx_m1; j++)
                    {
                        (*_asian_tmp)[(i * t_idx_m1) + j] = a_grid[l_idx] + (step * static_cast<T>(j));
                    }
                }
                
                /* This loop will get the remaining points */
                const unsigned int part_iters = a_grid.size() - (whole_iters * t_idx);
                const unsigned int l_idx = a_grid.size() - part_iters - 1;
                const T step = (a_grid.back() - a_grid[l_idx]) * t_idx_m1_inv;
                for (unsigned int i = 0; i < part_iters; i++)
                {
                    (*_asian_tmp)[(whole_iters * t_idx_m1) + i] = a_grid[l_idx] + (step * static_cast<T>(i));
                }
                _asian_tmp->back() = a_grid.back();
                
                permute_asian_values(a_grid, s_grid, v_grid, a_fac, t_idx_inv, t);
                g.replace_s_grid((*_asian_tmp), 0);
            }
            
            return _func.discontinious;
        }
                
    private :
        using cashflow<T>::get_s0_in;
        
        const S                 _func;
        std::vector<T>   *const _cross_asian_v;
        std::vector<T>   *const _asian_tmp;
        std::vector<T>   *const _tmp;
        const bool              _asian_in;
        
        void permute_asian_values(const std::vector<T> &a_grid, const std::vector<T> &s_grid,
            T *const v_grid, const T a_fac, const T s_fac, const T t) const
        {
            /* Permute values onto the new asian grid */
            _cross_asian_v->resize(a_grid.size());
            _tmp->resize(s_grid.size() * _asian_tmp->size());
            for (unsigned int i = 0; i < s_grid.size(); i++)
            {
                /* Copy v for this s_value from across all asian values */
                for (unsigned int j = 0; j < a_grid.size(); j++)
                {
                    (*_cross_asian_v)[j] = v_grid[(j * s_grid.size()) + i];
                }
                cubic_spline_interpolator<T> ci(a_grid.data(), _cross_asian_v->data(), a_grid.size(), true, 0);
                
                /* I would fix at s_value */
                const T s_value = s_grid[i];
                for (unsigned int j = 0; j < _asian_tmp->size(); j++)
                {
                    /* Before fixing I was at prefix_asian_value */
                    const T prefix_asian_value = (*_asian_tmp)[j];
                
                    /* After fixing I was at postfix_asian_value */
                    const T postfix_asian_value = (prefix_asian_value * a_fac) + (s_value * s_fac);
                    
                    /* Find and interpolate around the postfix_asian_value */
                    T inter;
                    if (postfix_asian_value > a_grid.back())
                    {
                        inter = linear_extrapolation(a_grid.data(), _cross_asian_v->data(), a_grid.size(), postfix_asian_value);
                    }
                    else
                    {
                        inter = ci.interpolate(postfix_asian_value);
                    }
                
                    /* Permute from old to new asian index, same s index */
                    (*_tmp)[(j * s_grid.size()) + i] = inter;
                }
            }
            
            /* Copy results over */
            memcpy(v_grid, &_tmp->front(), _tmp->size() * sizeof(T));
            _tmp->clear();
        }
};


/* Single asset based cashflows */
template<class T, class S>
class one_asset_dividend_cashflow : public cashflow<T>
{
    public :
        one_asset_dividend_cashflow(const S func, const T t, const T s0_in)
            : cashflow<T>(t, s0_in, func.always_eval), _func(func), _tmp(new std::vector<T>())
            {  };
            
        ~one_asset_dividend_cashflow()
        {
            delete _tmp;
        }

        void dimensions(std::vector<bool> *const d) const
        {
            /* Request 1 connected dimension */
            assert(d->empty() || !(*d)[0] || !"Error: Required dimension is already registered as disconnected.");
            if (d->empty())
            {
                d->push_back(false);
            }
            
            return;
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            /* No stock prices are required for dividends, just a dense enough grid */
            return;
        }
        
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return NO_BOUNDARY_CONDITION;
        }

        bool add(grid<T> &g, const T t) const
        {
            assert((g.dimensions() == 1) || ((g.dimensions() == 2) && g.outer_grid_disconnected()));
        
            /* Get the grid sizes */
            const int dim_m1        = g.dimensions() - 1;
            const int nr_of_grids   = (dim_m1 == 0) ? 1 : g.size(0);
            
            /* For each disconnected grid process the cashflow */
            const std::vector<T> &s = g.s_grid(dim_m1);
            T* v_grid = g.write_v_grid();
            for (int i = 0; i < nr_of_grids; i++)
            {
                cubic_spline_interpolator<T> ci(s.data(), v_grid, s.size(), true, 0);
                
                for (unsigned int j = 0; j < s.size(); j++)
                {
                    _tmp->push_back(_func(s, ci, j));
                }
            }
            
            /* Commit results, could save memory by doing this in the outer loop */
            memcpy(v_grid, &_tmp->front(), _tmp->size() * sizeof(T));
            _tmp->clear();
            
            
            return _func.discontinious;
        }
                
    private :
        const S                     _func;
        std::vector<T>   *const     _tmp;
};


/* Dividends */
template<class T>
class fixed_dividend
{
    public :
        fixed_dividend(const T div) : _div(div) {  };
        
        T operator()(const std::vector<T> &s, cubic_spline_interpolator<T> &ci, const int i) const
        {
            const T div_s = std::max(s[i] - _div, 0.0);
            return ci.interpolate(div_s);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const T _div;  
};


template<class T>
class percent_dividend
{
    public :
        percent_dividend(const T div_fac) : _div_fac(div_fac) {  };
        
        T operator()(const std::vector<T> &s, cubic_spline_interpolator<T> &ci, const int i) const
        {
            const T div_s = std::max(s[i] - (_div_fac * s[i]), 0.0);
            return ci.interpolate(div_s);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const T _div_fac;  
};

/* Two asset based cashflows */
/*template<class T, class S>
class two_asset_payoff_cashflow : public cashflow<T>
{
    public :
        two_asset_payoff_cashflow(const S &func, const T k_s1, const T k_s2, 
            const T s1_inc, const T s2_inc, const int nas1, const int nas2)
            : _func(func), _k_s1(k_s1), _k_s2(k_s2), _s1_inc(s1_inc), _s2_inc(s2_inc), (nas1), _nas2(nas2)
            {  };

        void dimensions(std::vector<bool> *const d) const
        {
*/            /* Request 2 connected dimension */
/*            assert(d->empty() || !(*d)[0] || !"Error: Required dimension is already registered as disconnected.");
            if (d->empty())
            {
                d->push_back(false);
            }
            
            assert((d->size() < 2) || !(*d)[1] || !"Error: Required dimension is already registered as disconnected.");
            if (d->size() < 2)
            {
                d->push_back(false);
            }
            
            return;
        }
            
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
        
        }        
        
        bool boundary_conditions(std::vector<T> &grid, T *const value)
        {
        
        }

        bool add(grid<T> &g, const T t) const
        {
            assert(g.size() == 2);
            
            int payoff_offset = 0;
            const std::vector<T> &s1 = g[0]->s_grid();
            const std::vector<T> &s2 = g[1]->s_grid();
            for (int i = 0; i < _nas1; i++)
            {
                for (int j = 0; j < _nas2; j++)
                {
                    fd_grid[payoff_offset++] = _func(s1[i], s2[j], _k_s1, _k_s2);
                }
            }
            
            return _func.discontinious;
        }
                
    private :
        const S    &_func;
        const T     _k_s1;
        const T     _k_s2;
        const T     _s1_inc;
        const T     _s2_inc;
        const int   _nas1;
        const int   _nas2;
};*/


/* Call option paying off on the best of two assets */
template<class T>
class bo_call_cashflow_func
{
    public :
        T operator()(const T v, const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return v + std::max(std::max(s1 - k_s1, 0.0), std::max(s2 - k_s2, 0.0));
        }
        
        const static bool discontinious = false;
};


/* Call option paying off on the worst of two assets */
template<class T>
class wo_call_cashflow_func
{
    public :
        T operator()(const T v, const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return v + std::min(std::max(s1 - k_s1, 0.0), std::max(s2 - k_s2, 0.0));
        }
        
        const static bool discontinious = false;
};


/* Binary call paying off if both asset above their strike */
template<class T>
class wo_binary_call_cashflow_func
{
    public :
        T operator()(const T v, const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return v + (((s1 > k_s1) && (s2 > k_s2)) ? 1.0 : 0.0);
        }
        
        const static bool discontinious = true;
};

#endif
