#ifndef __SINGLE_ASSET_CASHFLOW_H__
#define __SINGLE_ASSET_CASHFLOW_H__

#include "barrier_cashflow.h"


/* Payments */
template<class T>
class fixed_payment
{
    public :
        fixed_payment(const barrier_cashflow<T> barr, const T pay)
            : _barr(barr), _pay(pay) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr(std::max(v + _pay, 0.0), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            /* No stock prices are required for a fixed payment */
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const barrier_cashflow<T>   _barr;
        const T                     _pay;
};


template<class T>
class percent_payment
{
    public :
        percent_payment(const barrier_cashflow<T> barr, const T pay_fac) 
            : _barr(barr), _pay_fac(pay_fac) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr(std::max(v + (s1_out * _pay_fac), 0.0), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            /* No stock prices are required for a percent payment just a dense enough grid */
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const barrier_cashflow<T>   _barr;
        const T                     _pay_fac;  
};


/* Call option */
template<class T>
class euro_call_cashflow_func
{
    public :
        euro_call_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr((v + std::max(s1_out - (s1_in * _k), 0.0)), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = false;
    
    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};


template<class T>
class amer_call_cashflow_func
{
    public :
        amer_call_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr(v + std::max(std::max(s1_out - (s1_in * _k), v), 0.0), s1_out);
        }        
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = true;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};


/* Put option */
template<class T>
class euro_put_cashflow_func
{
    public :
        euro_put_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr(v + std::max((s1_in * _k) - s1_out, 0.0), s1_out);
        }
        
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = false;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};

template<class T>
class amer_put_cashflow_func
{
    public :
        amer_put_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            return _barr(v + std::max(std::max((s1_in * _k) - s1_out, v), 0.0), s1_out);
        }        
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = true;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};


/* Binary call */
template<class T>
class euro_binary_call_cashflow_func
{
    public :
        euro_binary_call_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };

        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            /* Valid for American and European options */
            return _barr(v + ((s1_out > (s1_in * _k)) ? 1.0 : 0.0), s1_out);
        }        
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = false;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};

template<class T>
class amer_binary_call_cashflow_func
{
    public :
        amer_binary_call_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };

        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            /* Valid for American and European options */
            return _barr(v + ((s1_out > (s1_in * _k)) ? 1.0 : 0.0), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = true;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};


/* Binary put */
template<class T>
class euro_binary_put_cashflow_func
{
    public :
        euro_binary_put_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            /* Valid for American and European options */
            return _barr(v + ((s1_out < (s1_in * _k)) ? 1.0 : 0.0), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = false;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};

template<class T>
class amer_binary_put_cashflow_func
{
    public :
        amer_binary_put_cashflow_func(const barrier_cashflow<T> barr, const T k) 
            : _barr(barr), _k(k) {  };
        
        T operator()(const T v, const T s1_in, const T s1_out) const
        {
            /* Valid for American and European options */
            return _barr(v + ((s1_out < (s1_in * _k)) ? 1.0 : 0.0), s1_out);
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            iter = _k;
            _barr.required_stock_price(iter);
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            return _barr.boundary_conditions(grid, value);
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = true;

    private :
        const barrier_cashflow<T>   _barr;
        const T                     _k;
};

#endif
