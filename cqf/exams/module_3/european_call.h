#ifndef __EUROPEAN_CALL_H__
#define __EUROPEAN_CALL_H__

#include <cmath>

#include "european_call_put.h"
#include "standard_normal.h"

using std::exp;
using std::log;
using std::sqrt;


double european_call_value(const double s, const double k, const double v, const double r, const double t)
{
    const double v_sqrt_t = v * sqrt(t);
    const double d1 = (log(s / k) + (r + (v * v * 0.5)) * t) / v_sqrt_t;
    const double d2 = d1 - v_sqrt_t;
    
    return (s * standard_normal_cdf(d1)) - (k * exp(-r * t) * standard_normal_cdf(d2));
}


#endif
