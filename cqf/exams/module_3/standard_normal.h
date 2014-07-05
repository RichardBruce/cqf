#ifndef __STANDARD_NORMAL_H__
#define __STANDARD_NORMAL_H__

#include <cmath>

#include "constants.h"


double standard_normal_pdf(double x)
{
    return SQRT_2PI_INV * std::exp(-0.5 * x * x);
}

/* Abramowiz and Stegun's approximation */
double standard_normal_cdf(double x)
{
    /* Too small and too big */
    if (x < -6.0)
    {
        return 0.0;
    }
    
    if (x > 6.0)
    {
        return 1.0;
    }
    
    /* Co-efficients */
    const double b0 = 0.31938153;
    const double b1 = -0.356563782;
    const double b2 = 1.781477937;
    const double b3 = -1.821255978;
    const double b4 = 1.330274429;
    const double p = 0.2316419;
    const double c = 0.3989423;
    
    const double a = fabs(x);
    const double t = 1.0 / (1.0 + (a * p));
    const double b = c * std::exp(-x * (x * 0.5));
    const double n = (((((b4 * t) + b3) * t + b2) * t + b1) * t + b0) * t;
    const double half_cdf = 1.0 - b * n;
    if (x < 0.0)
    {
        return 1.0 - half_cdf;
    }
    
    return half_cdf;
}

#endif
