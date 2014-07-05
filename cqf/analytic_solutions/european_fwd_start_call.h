#ifndef __EUROPEAN_FWD_START_CALL_H__
#define __EUROPEAN_FWD_START_CALL_H__

#include <cmath>

#include "standard_normal.h"

template<class T>
T european_fwd_start_call_value(const T s, const T k, const T v, const T r, const T b, const T t_a, const T t_b)
{
    const T tenor = t_b - t_a;
    const T mu = b - r;
    const T v_sqrt_t = v * std::sqrt(tenor);
    const T d1 = (std::log(1.0 / k) + ((b + (v * v * 0.5)) * tenor)) / v_sqrt_t;
    const T d2 = d1 - v_sqrt_t;
    
    return ((s * std::exp(mu * t_a)) * 
            ((std::exp(mu * tenor) * standard_normal_cdf(d1)) -
             (std::exp(-r * tenor) * k * standard_normal_cdf(d2))));
}

#endif
