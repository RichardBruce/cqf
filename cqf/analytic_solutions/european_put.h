#ifndef __EUROPEAN_PUT_H__
#define __EUROPEAN_PUT_H__

#include <cmath>

#include "european_call_put.h"
#include "standard_normal.h"

using std::exp;
using std::log;
using std::sqrt;


template<class T>
T european_put_value(const T s, const T k, const T v, const T r, const T t)
{
    const T v_sqrt_t = v * sqrt(t);
    const T d1 = (log(s / k) + (r + (v * v * 0.5)) * t) / v_sqrt_t;
    const T d2 = d1 - v_sqrt_t;
    
    return (s * -standard_normal_cdf(-d1)) + (k * exp(-r * t) * standard_normal_cdf(-d2));
}


template<class T>
T european_put_delta(const T s, const T k, const T v, const T r, const T t)
{
    const T d1 = european_call_put_d1(s, k, v, r, t);
    return -standard_normal_cdf(-d1);
}


template<class T>
T european_put_theta(const T s, const T k, const T v, const T r, const T t)
{
    const T v_sqrt_t = v * sqrt(t);
    const T d1 = (log(s / k) + (r + (v * v * 0.5)) * t) / v_sqrt_t;
    const T d2 = d1 - v_sqrt_t;
    
    const T pt1 = -(s * v * standard_normal_pdf(d1)) / (2.0 * sqrt(t));
    const T pt2 = r * k * exp(-r * t) * standard_normal_cdf(-d2);
    
    return pt1 + pt2;  
}


template<class T>
T european_put_rho(const T s, const T k, const T v, const T r, const T t)
{
    const T d2 = european_call_put_d2(s, k, v, r, t);
    return -k * t * exp(-r * t) * standard_normal_pdf(-d2);
}


#endif
