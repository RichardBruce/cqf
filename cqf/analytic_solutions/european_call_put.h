#ifndef __EUROPEAN_CALL_PUT_H__
#define __EUROPEAN_CALL_PUT_H__

#include <cmath>

#include "standard_normal.h"

using std::log;
using std::sqrt;


template<class T>
T european_call_put_d1(const T s, const T k, const T v, const T r, const T t)
{
    return (log(s / k) + (r + (v * v * 0.5f)) * t) / (v * sqrt(t));
}


template<class T>
T european_call_put_d2(const T s, const T k, const T v, const T r, const T t)
{
    return (log(s / k) + (r - (v * v * 0.5f)) * t) / (v * sqrt(t));
}


template<class T>
T european_call_put_gamma(const T s, const T k, const T v, const T r, const T t)
{
    const T d1 = european_call_put_d1(s, k, v, r, t);
    return standard_normal_pdf(d1) / (s * v * sqrt(t));
}


template<class T>
T european_call_put_vega(const T s, const T k, const T v, const T r, const T t)
{
    const T d1 = european_call_put_d1(s, k, v, r, t);
    return s * sqrt(t) * standard_normal_pdf(d1);
}

#endif
