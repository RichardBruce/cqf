#ifndef __EUROPEAN_TWO_ASSET_CASH_BINARY_H__
#define __EUROPEAN_TWO_ASSET_CASH_BINARY_H__

#include <cmath>

#include "bivariate_standard_normal.h"

template<class T>
T european_two_asset_cash_binary_value(
    const T s1, const T k1, const T v1,  const T b1,
    const T s2, const T k2, const T v2,  const T b2, 
    const T r,  const T t,  const T rho, const T c, 
    bool s1_call, bool s2_call)
{
    const T sqrt_t = std::sqrt(t);
    const T v1_sqrt_t = v1 * sqrt_t;
    const T v2_sqrt_t = v2 * sqrt_t;
    
    const T d1 = (log(s1 / k1) + (b1 - (v1 * v1 * 0.5)) * t) / v1_sqrt_t;
    const T d2 = (log(s2 / k2) + (b2 - (v2 * v2 * 0.5)) * t) / v2_sqrt_t;

    const T d1_adj  = s1_call ? d1 : -d1;
    const T d2_adj  = s2_call ? d2 : -d2;
    const T rho_adj = (s1_call == s2_call) ? rho : -rho;
        
    return (c * exp(-r * t) * bivariate_standard_normal_cdf<T>(d1_adj, d2_adj, rho_adj));
}

#endif
