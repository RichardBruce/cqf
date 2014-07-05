#ifndef __EUROPEAN_MIN_OF_TWO_CALL_H__
#define __EUROPEAN_MIN_OF_TWO_CALL_H__

#include <cmath>

#include "bivariate_standard_normal.h"

template<class T>
T european_min_of_two_call_value(
    const T s1, const T v1,  const T b1,
    const T s2, const T v2,  const T b2, 
    const T r,  const T t,  const T rho, const T k)
{
    const T sqrt_t      = std::sqrt(t);
    const T v1_sqrt_t   = v1 * sqrt_t;
    const T v2_sqrt_t   = v2 * sqrt_t;
    const T v1_sq       = v1 * v1;
    const T v2_sq       = v2 * v2;
    
    const T v           = std::sqrt(v1_sq + v2_sq - (2.0 * rho * v1 * v2));
    const T v_sqrt_t    = v * sqrt_t;
    
    const T d = (log(s1 / s2) + (b1 - b2 + (0.5 * v * v)) * t) / v_sqrt_t;
    
    const T rho1    = (v1 - (rho * v2)) / v;
    const T rho2    = (v2 - (rho * v1)) / v;
    const T y1      = (log(s1 / k) + (b1 + (v1_sq * 0.5)) * t) / v1_sqrt_t;
    const T y2      = (log(s2 / k) + (b2 + (v2_sq * 0.5)) * t) / v2_sqrt_t;

//    const T d1_adj  = s1_call ? d1 : -d1;
//    const T d2_adj  = s2_call ? d2 : -d2;
//    const T rho_adj = (s1_call == s2_call) ? rho : -rho;
        
    return (s1 * exp((b1 - r) * t) * bivariate_standard_normal_cdf<T>(y1,            -d,              -rho1)) +
           (s2 * exp((b2 - r) * t) * bivariate_standard_normal_cdf<T>(y2,             d  - v_sqrt_t,  -rho2)) +
           (k  * exp(-r * t)       * bivariate_standard_normal_cdf<T>(y1 - v1_sqrt_t, y2 - v2_sqrt_t,  rho));
}

#endif
