#ifndef __CASHFLOW_H__
#define __CASHFLOW_H__

#include "boost/noncopyable.hpp"
#include "boost/bind.hpp"

#include "interpolation.h"


/* Abstract base for the application of a cash flow that may be a fucntion of s and t */
/* Final conditions are also cash flows (of the size of the derivative payoff */
/* Should return true if the cashflow has a dicontinuity in s */
template<class T>
class cashflow : private boost::noncopyable
{
    public :
        cashflow(const T t = 1.0, const bool always_eval = false) 
            : _t(t), _always_eval(always_eval) {  };
            
        /* Virtual DTOR, this class will be used as a base class */
        virtual ~cashflow() {  };
        
        
        /* Common functions */
        bool exists(const T t) const
        {
            return _always_eval || (_t == t);
        }
        
        void required_date(std::back_insert_iterator<std::vector<T> > iter) const
        {
            (*iter) = _t;
        }
        
        virtual bool add(T *const fd_grid, const T t) const = 0;
        virtual bool add(T *const fd_grid, const std::vector<grid<T>*> &g, const T t) const = 0;
    
    private :
        const T     _t;
        const bool  _always_eval;
};


/* Single asset based cashflows */
template<class T, class S>
class one_asset_payoff_cashflow : public cashflow<T>
{
    public :
        one_asset_payoff_cashflow(const S &func, const T t, const T s1_inc, const int nas1)
            : cashflow<T>(t, func.always_eval), _func(func), _s1_inc(s1_inc), _nas1(nas1)
            {  };

        bool add(T *const fd_grid, const T t) const
        {
            int payoff_offset = 0;
            for (int i = 0; i < _nas1; i++)
            {
                const T s1 = static_cast<T>(i) * _s1_inc;
                fd_grid[payoff_offset] = _func(fd_grid[payoff_offset], s1);
                ++payoff_offset;
            }
            
            return _func.discontinious;
        }
        
        bool add(T *const fd_grid, const std::vector<grid<T>*> &g, const T t) const
        {
            assert(g.size() == 1);
            
            int payoff_offset = 0;
            const std::vector<T> &s = g[0]->points();
            for (unsigned int i = 0; i < s.size(); i++)
            {
                fd_grid[payoff_offset] = _func(fd_grid[payoff_offset], s[i]);
                ++payoff_offset;
            }
            
            return _func.discontinious;
        }
                
    private :
        const S    &_func;
        const T     _s1_inc;
        const int   _nas1;
};


/* Payments */
template<class T>
class fixed_payment
{
    public :
        fixed_payment(const T pay) : _pay(pay) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(v + _pay, 0.0);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const T _pay;  
};


template<class T>
class percent_payment
{
    public :
        percent_payment(const T pay_fac) : _pay_fac(pay_fac) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(v + (s1 * _pay_fac), 0.0);
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const T _pay_fac;  
};


/* Call option */
template<class T>
class euro_call_cashflow_func
{
    public :
        euro_call_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(s1 - _k, 0.0);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = false;
    
    private :
        const T _k;
};

template<class T>
class amer_call_cashflow_func
{
    public :
        amer_call_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(std::max(s1 - _k, v), 0.0);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = true;

    private :
        const T _k;
};


/* Put option */
template<class T>
class euro_put_cashflow_func
{
    public :
        euro_put_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(_k - s1, 0.0);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = false;

    private :
        const T _k;
};

template<class T>
class amer_put_cashflow_func
{
    public :
        amer_put_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            return std::max(std::max(_k - s1, v), 0.0);
        }
        
        const static bool discontinious = false;
        const static bool always_eval   = true;

    private :
        const T _k;
};


/* Binary call */
template<class T>
class euro_binary_call_cashflow_func
{
    public :
        euro_binary_call_cashflow_func(const T k) : _k(k) {  };

        T operator()(const T v, const T s1) const
        {
            /* Valid for American and European options */
            return (s1 > _k) ? 1.0 : 0.0;
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = false;

    private :
        const T _k;
};

template<class T>
class amer_binary_call_cashflow_func
{
    public :
        amer_binary_call_cashflow_func(const T k) : _k(k) {  };

        T operator()(const T v, const T s1) const
        {
            /* Valid for American and European options */
            return (s1 > _k) ? 1.0 : 0.0;
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = true;

    private :
        const T _k;
};


/* Binary put */
template<class T>
class euro_binary_put_cashflow_func
{
    public :
        euro_binary_put_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            /* Valid for American and European options */
            return (s1 < _k) ? 1.0 : 0.0;
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = false;

    private :
        const T _k;
};

template<class T>
class amer_binary_put_cashflow_func
{
    public :
        amer_binary_put_cashflow_func(const T k) : _k(k) {  };
        
        T operator()(const T v, const T s1) const
        {
            /* Valid for American and European options */
            return (s1 < _k) ? 1.0 : 0.0;
        }
        
        const static bool discontinious = true;
        const static bool always_eval   = true;

    private :
        const T _k;
};


/* Single asset based cashflows */
template<class T, class S>
class one_asset_dividend_cashflow : public cashflow<T>
{
    public :
        one_asset_dividend_cashflow(const S &func, const T t)
            : cashflow<T>(t, func.always_eval), _func(func), _tmp(new std::vector<T>())
            {  };
            
        ~one_asset_dividend_cashflow()
        {
            delete _tmp;
        }

        /* Depricated will not implement */
        bool add(T *const fd_grid, const T t) const
        {
            assert(false);
        }
        
        bool add(T *const fd_grid, const std::vector<grid<T>*> &g, const T t) const
        {
            assert(g.size() == 1);
            
            const std::vector<T> &s = g[0]->points();
            for (unsigned int i = 0; i < s.size(); i++)
            {
                _tmp->push_back(_func(s, fd_grid, i));
            }
            
            memcpy(&fd_grid[0], &_tmp->front(), _tmp->size() * sizeof(T));
            _tmp->clear(); 
            return _func.discontinious;
        }
                
    private :
        const S                &_func;
        std::vector<T>   *const _tmp;
};


/* Dividends */
template<class T>
class fixed_dividend
{
    public :
        fixed_dividend(const T div) : _div(div) {  };
        
        T operator()(const std::vector<T> &s, T *const v, const int i) const
        {
            const T div_s = std::max(s[i] - _div, 0.0);
            typename std::vector<T>::const_iterator iter = std::find_if(s.begin(), s.end(), boost::bind(std::less_equal<T>(), div_s, _1));
            const int v_idx = std::distance(s.begin(), iter);
            
            /* Mid section */
            if ((v_idx > 1) && (v_idx < static_cast<int>(s.size() - 1)))
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cubic_interpolation(v[v_idx - 2], v[v_idx - 1], v[v_idx], v[v_idx + 1], mu);
            }
            /* Lower points */
            else if (v_idx == 1)
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cosine_interpolation(v[0], v[1], mu);
            }
            else if (v_idx == 0)
            {
                return v[0];
            }
            /* Upper points */
            else
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cosine_interpolation(v[v_idx - 1], v[v_idx], mu);
            }
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
        
        T operator()(const std::vector<T> &s, T *const v, const int i) const
        {
            const T div_s = std::max(s[i] - (_div_fac * s[i]), 0.0);
            typename std::vector<T>::const_iterator iter = std::find_if(s.begin(), s.end(), boost::bind(std::less_equal<T>(), div_s, _1));
            const int v_idx = std::distance(s.begin(), iter);
            
            /* Mid section */
            if ((v_idx > 1) && (v_idx < static_cast<int>(s.size() - 1)))
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cubic_interpolation(v[v_idx - 2], v[v_idx - 1], v[v_idx], v[v_idx + 1], 
                    s[v_idx - 2], s[v_idx - 1], s[v_idx], s[v_idx + 1], mu);
            }
            /* Lower points */
            else if (v_idx == 1)
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cosine_interpolation(v[0], v[1], mu);
            }
            else if (v_idx == 0)
            {
                return v[0];
            }
            /* Upper points */
            else
            {
                const T mu = (div_s - s[v_idx - 1]) / (s[v_idx] - s[v_idx - 1]);
                return cosine_interpolation(v[v_idx - 1], v[v_idx], mu);
            }
        }
        
        const static bool discontinious = false;  
        const static bool always_eval   = false;
        
    private :
        const T _div_fac;  
};

/* Two asset based cashflows */
template<class T, class S>
class two_asset_payoff_cashflow : public cashflow<T>
{
    public :
        two_asset_payoff_cashflow(const S &func, const T k_s1, const T k_s2, 
            const T s1_inc, const T s2_inc, const int nas1, const int nas2)
            : _func(func), _k_s1(k_s1), _k_s2(k_s2), _s1_inc(s1_inc), _s2_inc(s2_inc), _nas1(nas1), _nas2(nas2)
            {  };

        bool add(T *const fd_grid, const T t) const
        {
            int payoff_offset = 0;
            for (int i = 0; i < _nas1; i++)
            {
                const T s1 = static_cast<T>(i) * _s1_inc;
                for (int j = 0; j < _nas2; j++)
                {
                    const T s2 = static_cast<T>(j) * _s2_inc;
                    fd_grid[payoff_offset++] = _func(s1, s2, _k_s1, _k_s2);
                }
            }
            
            return _func.discontinious;
        }
        
        bool add(T *const fd_grid, const std::vector<grid<T>*> &g, const T t) const
        {
            assert(g.size() == 2);
            
            int payoff_offset = 0;
            const std::vector<T> &s1 = g[0]->points();
            const std::vector<T> &s2 = g[1]->points();
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
};


/* Call option paying off on the best of two assets */
template<class T>
class bo_call_cashflow_func
{
    public :
        T operator()(const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return std::max(std::max(s1 - k_s1, 0.0), std::max(s2 - k_s2, 0.0));
        }
        
        const static bool discontinious = false;
};


/* Call option paying off on the worst of two assets */
template<class T>
class wo_call_cashflow_func
{
    public :
        T operator()(const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return std::min(std::max(s1 - k_s1, 0.0), std::max(s2 - k_s2, 0.0));
        }
        
        const static bool discontinious = false;
};


/* Binary call paying off if both asset above their strike */
template<class T>
class wo_binary_call_cashflow_func
{
    public :
        T operator()(const T s1, const T s2, const T k_s1, const T k_s2) const
        {
            return ((s1 > k_s1) && (s2 > k_s2)) ? 1.0 : 0.0;
        }
        
        const static bool discontinious = true;
};

#endif
