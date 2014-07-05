#ifndef __GRID_MAPPER_H__
#define __GRID_MAPPER_H__

#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>

#include "boost/bind.hpp"
#include "boost/noncopyable.hpp"
#include "boost/shared_ptr.hpp"

#include "utility.h"

#include "cashflow.h"


/* Mapper abstract base class */
template<class T>
class grid_mapper
{
    public :
        /* CTOR */
        grid_mapper(const T min_fac, const T max_fac) : _min_fac(min_fac), _max_fac(max_fac) {  };
        
        /* DTOR */
        virtual ~grid_mapper() {  };
    
        /* Map the grid points to a new grid */
        virtual std::vector<T>* map(const std::vector<cashflow<T> *> &cashflows, const T s,
            std::pair<bool, T> *const lb, std::pair<bool, T> *const ub) const = 0;
        
        /* Helper function */
        void build_required_cashflows(std::vector<T> &required, 
            const std::vector<cashflow<T> *> &cashflows, const T s, 
            std::pair<bool, T> *const lb, std::pair<bool, T> *const ub, 
            const bool inc_minmax = true, const bool inc_spot = true) const 
        {
            /* Get all prices that must be calculated */
            std::back_insert_iterator<std::vector<T>> required_inserter(required);
            for(unsigned int i = 0; i < cashflows.size(); i++)
            {
                cashflows[i]->required_stock_price(required_inserter);
            }
            
            /* The stock price is also required so as not to interpolate the results */
            /* This could possibly be relaxed */
            if (inc_spot)
            {
                required_inserter = s;
            }
            
            /* Sort */
            std::sort(required.begin(), required.end());
            
            /* s_min and s_max are required so that values between are filled in */
            if (inc_minmax)
            {
                /* ooo, push_front */
                required.insert(required.begin(), _min_fac * required.front());
                required.push_back(_max_fac * required.back());
            }
            
            /* Clean up any duplicates */
            auto required_iter = std::unique(required.begin(), required.end());
            required.resize(required_iter - required.begin());
            
            /* Apply any boundary conditions */
            lb->first = false;
            ub->first = false;
            T boundary_value;
            bool applied_below = false;
            bool applied_above = false;
            for(unsigned int i = 0; i < cashflows.size(); i++)
            {
                const boundary_condition_t cond = cashflows[i]->boundary_conditions(required, &boundary_value);
                if (cond == LOWER_BOUNDARY_CONDITION)
                {
                    assert(!applied_below || !"Error: Cannot apply multi boundary conditons to the same boundary.");
                    lb->first = true;
                    lb->second = boundary_value;
                }
                else if (cond == UPPER_BOUNDARY_CONDITION)
                {
                    assert(!applied_above || !"Error: Cannot apply multi boundary conditons to the same boundary.");
                    ub->first = true;
                    ub->second = boundary_value;
                }
            }
        }
        
    protected :
        const T _min_fac;
        const T _max_fac;
};

/* Mapper to assign a uniformed grid */
template<class T>
class uniformed_grid_mapper : public grid_mapper<T>
{
    public : 
        /* CTOR */
        uniformed_grid_mapper(const T max_fac, const T min_fac, const T max_inc)
         : grid_mapper<T>(min_fac, max_fac), _max_inc(max_inc) { }
        
        std::vector<T>* map(const std::vector<cashflow<T> *> &cashflows, const T s,
            std::pair<bool, T> *const lb, std::pair<bool, T> *const ub) const 
        {
            std::vector<T> required;
            this->build_required_cashflows(required, cashflows, s, lb, ub);
    
            /* Pad the stock prices for the minimum step size */
            std::vector<T> *values = new std::vector<T>();
            pad_data<T>(values, required, _max_inc);
            
            return values;
        }

    private :
        const T _max_inc;
};

/* Mapper to assign a non-uniformed grid */
template<class T>
class sinh_grid_mapper : public grid_mapper<T>
{
    public :
        /* CTOR */
        /* k is a point of interest */
        /* c is the proportion of grid points that lie in the region of k */
        sinh_grid_mapper(const T k, const T c, const T max_fac, const T min_fac, const T avg_inc)
         : grid_mapper<T>(min_fac, max_fac), _k(k), _c(c), _avg_inc(avg_inc) { }
               
        std::vector<T>* map(const std::vector<cashflow<T> *> &cashflows, const T s,
            std::pair<bool, T> *const lb, std::pair<bool, T> *const ub) const 
        {
            std::vector<T> required;
            this->build_required_cashflows(required, cashflows, s, lb, ub);
            
            auto values = new std::vector<T>(static_cast<int>((required.back() - required.front()) / _avg_inc) + 1);
            const T nas_t = static_cast<T>(values->size());
            for (unsigned int i = 0; i < values->size(); i++)
            {
                const T uni = std::asinh(-_k / _c) + ((i / nas_t) * (std::asinh((required.back() - _k) / _c) - std::asinh(-_k / _c)));
                (*values)[i] = _k + (_c * std::sinh(uni));
            }

            return values;
        }
        
    private :
        const T _k;
        const T _c;
        const T _avg_inc;
};


#endif
