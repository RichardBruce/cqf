#ifndef __EUROPEAN_CASH_BINARY_CALL_H__
#define __EUROPEAN_CASH_BINARY_CALL_H__

#include <cmath>

#include "standard_normal.h"

using std::exp;
using std::log;
using std::sqrt;


template<class T>
T european_cash_binary_call_value(const T s, const T k, const T v, const T r, const T t)
{
    const T v_sqrt_t = v * sqrt(t);
    const T d = (log(s / k) + (r - (v * v * 0.5)) * t) / v_sqrt_t;
    
    return (k * exp(-r * t) * standard_normal_cdf(d));
}

#endif
