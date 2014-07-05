#ifndef __EUROPEAN_ASSET_BINARY_PUT_H__
#define __EUROPEAN_ASSET_BINARY_PUT_H__

#include <cmath>

#include "standard_normal.h"

using std::exp;
using std::log;
using std::sqrt;


float european_asset_binary_put_value(const float s, const float k, const float v, const float r, const float t)
{
    const float v_sqrt_t = v * sqrt(t);
    const float d = (log(s / k) + (r + (v * v * 0.5f)) * t) / v_sqrt_t;
    
    return (s * standard_normal_cdf(-d));
}

#endif
