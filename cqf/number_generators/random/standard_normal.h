#ifndef __STANDARD_NORMAL_H__
#define __STANDARD_NORMAL_H__

#include <cmath>

#include <assert.h>

#include "constants.h"


template<class T>
T standard_normal_pdf(const T x)
{
    return (1.0f / std::sqrt(2.0f * PI)) * std::exp(-0.5f * x * x);
}

/* Abramowiz and Stegun's approximation */
template<class T>
T standard_normal_cdf(const T x)
{
    /* Too small and too big */
    if (x < -6.0f)
    {
        return 0.0f;
    }
    
    if (x > 6.0f)
    {
        return 1.0f;
    }
    
    /* Co-efficients */
    static const T b0 = 0.31938153f;
    static const T b1 = -0.356563782;
    static const T b2 = 1.781477937;
    static const T b3 = -1.821255978;
    static const T b4 = 1.330274429;
    static const T p = 0.2316419;
    static const T c = 0.3989423;
    
    const T a = fabs(x);
    const T t = 1.0f / (1.0f + (a * p));
    const T b = c * std::exp(-x * (x * 0.5f));
    const T n = (((((b4 * t) + b3) * t + b2) * t + b1) * t + b0) * t;
    const T half_cdf = 1.0f - b * n;
    if (x < 0.0)
    {
        return 1.0f - half_cdf;
    }
    
    return half_cdf;
}

/* Error function approximation for cdf */
//template<class T>
//T stdnormal_cdf(const T u)
//{
//    static const T a0 = 0.01161110663653770;
//    static const T a1 = 0.3951404679838207;
//    static const T a2 = 28.46603853776254;
//    static const T a3 = 188.7426188426510;
//    static const T a4 = 3209.377589138469;
//
//    static const T b0 = 0.1767766952966369;
//    static const T b1 = 8.344316438579620;
//    static const T b2 = 172.5514762600375;
//    static const T b3 = 1813.893686502485;
//    static const T b4 = 8044.716608901563;
//
//    static const T c0 = 0.0000000215311535474403846;
//    static const T c1 = 0.564188496988670089;
//    static const T c2 = 8.88314979438837594;
//    static const T c3 = 66.1191906371416295;
//    static const T c4 = 298.635138197400131;
//    static const T c5 = 881.952221241769090;
//    static const T c6 = 1712.04761263407058;
//    static const T c7 = 2051.07837782607147;
//    static const T c8 = 1230.33935479799725;
//
//    static const T d0 = 1.0;
//    static const T d1 = 15.7449261107098347;
//    static const T d2 = 117.693950891312499;
//    static const T d3 = 537.181101862009858;
//    static const T d4 = 1621.38957456669019;
//    static const T d5 = 3290.79923573345963;
//    static const T d6 = 4362.61909014324716;
//    static const T d7 = 3439.36767414372164;
//    static const T d8 = 1230.33935480374942;
//
//    static const T p0 = 0.0163153871373020978;
//    static const T p1 = 0.305326634961232344;
//    static const T p2 = 0.360344899949804439;
//    static const T p3 = 0.125781726111229246;
//    static const T p4 = 0.0160837851487422766;
//    static const T p5 = 0.000658749161529837803;
//
//    static const T q0 = 1.0;
//    static const T q1 = 2.56852019228982242;
//    static const T q2 = 1.87295284992346047;
//    static const T q3 = 0.527905102951428412;
//    static const T q4 = 0.0605183413124413191;
//    static const T q5 = 0.00233520497626869185;
//
//    T y = fabs(u);
//    if (y <= (0.46875 * std::sqrt(2.0)))
//    {
//        /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
//        T z = y * y;
//        y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4) / ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
//        return 0.5 + y;
//    }
//    
//    T z = exp(-y * y * 0.5) * 0.5;
//    if (y <= 4.0)
//    {
//        /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
//        y = y * (1.0 / std::sqrt(2.0));
//        y = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8)
//            / ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);
//    
//        y = z * y;
//    }
//    else
//    {
//        /* evaluate erfc() for |u| > sqrt(2)*4.0 */
//        z = z * std::sqrt(2.0) / y;
//        y = 2.0 / (y * y);
//        y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
//       
//        y = z * (M_1_SQRTPI - y);
//    }
//    
//    return (u < 0.0) ? y : (1.0 - y);
//};

/* Acklam's approximation */
template<class T>
T inverse_standard_normal_cdf(const T p, const bool use_halley = false)
{
    /* Valid range for a probability */
    assert(p > 0.0);
    assert(p < 1.0);

    /* Constant coefficients */    
    static const T a0 = -39.69683028665376;
    static const T a1 = 220.9460984245205;
    static const T a2 = -275.9285104469687;
    static const T a3 = 138.3577518672690;
    static const T a4 = -30.66479806614716;
    static const T a5 = 2.506628277459239;

    static const T b0 = -54.47609879822406;
    static const T b1 = 161.5858368580409;
    static const T b2 = -155.6989798598866;
    static const T b3 = 66.80131188771972;
    static const T b4 = -13.28068155288572;
 
    static const T c0 = -0.007784894002430293;
    static const T c1 = -0.3223964580411365;
    static const T c2 = -2.400758277161838;
    static const T c3 = -2.549732539343734;
    static const T c4 = 4.374664141464968;
    static const T c5 = 2.938163982698783;
    
    static const T d0 = 0.007784695709041462;
    static const T d1 = 0.3224671290700398;
    static const T d2 = 2.445134137142996;
    static const T d3 = 3.754408661907416;

    /* Rational approximation */
    T t, u;
    const T q = std::min(p, 1 - p);
    if (q > 0.02425) 
    {
        u = q - 0.5;
        t = u * u;
        u = u * (((((a0 * t + a1) * t + a2) * t + a3) * t + a4) * t + a5) / 
            (((((b0 * t + b1) * t + b2) * t + b3) * t + b4) * t + 1);
    } 
    else 
    {
        t = std::sqrt(-2.0 * std::log(q));
        u = (((((c0 * t + c1) * t + c2) * t + c3) * t + c4) * t + c5) /
            ((((d0 * t + d1) * t + d2) * t + d3) * t + 1);
    }

    /* Opyional one iteration of Halley refinement */
    if (use_halley)
    {
        const T error = standard_normal_cdf<T>(u) - q;
        const T correction = error * std::sqrt(2.0 * PI) * std::exp(u * u * 0.5);
        u = u - correction / (1.0 + (u * t * 0.5));
    }
    
    return (p > 0.5 ? -u : u);
};

#endif
