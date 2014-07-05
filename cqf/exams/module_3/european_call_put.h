#ifndef __EUROPEAN_CALL_PUT_H__
#define __EUROPEAN_CALL_PUT_H__

#include <cmath>

#include "standard_normal.h"

using std::log;
using std::sqrt;


double european_call_put_d1(const double s, const double k, const double v, const double r, const double t)
{
    return (log(s / k) + (r + (v * v * 0.5)) * t) / (v * sqrt(t));
}


double european_call_put_vega(const double s, const double k, const double v, const double r, const double t)
{
    const double d1 = european_call_put_d1(s, k, v, r, t);
    return s * sqrt(t) * standard_normal_pdf(d1);
}

#endif
