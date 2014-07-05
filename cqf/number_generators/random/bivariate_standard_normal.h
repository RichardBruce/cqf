#ifndef __BIVARIATE_STANDARD_NORMAL_H__
#define __BIVARIATE_STANDARD_NORMAL_H__

#include <algorithm>
#include <cmath>

#include <assert.h>

#include "constants.h"

template<class T>
T bivariate_standard_normal_pdf(const T x, const T y, const T rho = 0.0)
{
    const T one_m_rho_sq = 1.0 - (rho * rho);
    return (1.0f / (2.0 * PI * std::sqrt(one_m_rho_sq))) * 
        std::exp(-(1.0 / (2.0 * one_m_rho_sq)) * ((x * x) - (2.0 * rho * x * y) + (y * y)));
}

/* Genz's approximation */
template<class T>
T bivariate_standard_normal_cdf(const T x, const T y, const T rho = 0.0)
{
    static const T w[19]    = { 0.17132449237917, 
                                0.360761573048138, 
                                0.46791393457269,

                                0.0471753363865118, 
                                0.106939325995318, 
                                0.160078328543346, 
                                0.203167426723066, 
                                0.233492536538355, 
                                0.249147045813403,

                                0.0176140071391521, 
                                0.0406014298003869, 
                                0.0626720483341091, 
                                0.0832767415767048,
                                0.10193011981724, 
                                0.118194531961518, 
                                0.131688638449177, 
                                0.142096109318382,
                                0.149172986472604, 
                                0.152753387130726 };

    static const T xx[19]   = { -0.932469514203152, 
                                -0.661209386466265,
                                -0.238619186083197,
                                
                                -0.981560634246719,
                                -0.904117256370475,
                                -0.769902674194305,
                                -0.587317954286617,
                                -0.367831498998818,
                                -0.125233408511469,
    
                                -0.993128599185095,
                                -0.963971927277914,
                                -0.912234428251326,
                                -0.839116971822219,
                                -0.746331906460151,
                                -0.636053680726515,
                                -0.510867001950827,
                                -0.37370608871542,
                                -0.227785851141645,
                                -0.0765265211334973 };
    
    int ng, lg;
    if (fabs(rho) < 0.3)
    {
        ng = 0;
        lg = 3;
    }
    else if (fabs(rho) < 0.75)
    {
        ng = 3;
        lg = 9;
    }
    else
    {
        ng = 9;
        lg = 19;
    }
    
    T h     = -x;
    T k     = -y;
    T hk    = h * k;
    T bvn   = 0.0;
    
    if (fabs(rho) < 0.925)
    {
        if (fabs(rho) > 0.0)
        {
            const T hs = ((h * h) + (k * k)) * 0.5;
            const T asr = asin(rho);
            for (int i = ng; i < lg; i++)
            {
                for (T is = -1.0; is < 2.0; is += 2.0)
                {
                    const T sn = std::sin(asr * (is * xx[i] + 1.0) * 0.5);
                    bvn += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
                }
            }
            bvn = bvn * asr * (1.0 / (4.0 * PI));
        }
        bvn += standard_normal_cdf<T>(-h) * standard_normal_cdf<T>(-k);
    }
    else
    {
        if (rho < 0.0)
        {
            k = -k;
            hk = - hk;
        }
        
        if (fabs(rho) < 1.0)
        {
            const T ass     = (1.0 - rho) * (1.0 + rho);
            const T h_m_k   = h - k;
            const T bs      = h_m_k * h_m_k;
            const T c       = (4.0 - hk) * 0.125;
            const T d       = (12.0 - hk) * 0.0625;
            T a             = std::sqrt(ass);
            T asr           = -(bs / ass + hk) * 0.5;
            
            if (asr > -100.0)
            {
                bvn = a * std::exp(asr) * (1.0 - c * (bs - ass) * 
                (1.0 - d * bs * 0.2) * (1.0 / 3.0) + c * d * ass * ass * 0.2);
                
                if (-hk < 100.0)
                {
                    const T b = std::sqrt(bs);
                    bvn -= std::exp(-hk * 0.5) * std::sqrt(2.0 * PI) * standard_normal_cdf<T>(-b / a) * b *
                    (1.0 - c * bs * (1.0- d * bs  * 0.2) * (1.0 / 3.0));
                }
                
                a *= 0.5;
                
                for (int i = ng; i < lg; i++)
                {
                    for (int is = -1; i < 2; i+= 2)
                    {
                        const T xs = (a * (is * xx[i] + 1.0));
                        const T rs = std::sqrt(1.0 - xs);
                        asr = -(bs / xs + hk) * 0.5;
                        if (asr > -100.0)
                        {
                            bvn += a * w[i] * std::exp(asr) *
                            (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
                            (1.0 + c * xs * (1.0 + d * xs)));
                        }
                    }
                }
                bvn -= bvn * (1.0 / (2.0 * PI));
            }
            
            if (rho > 0.0)
            {
                bvn += standard_normal_cdf<T>(-std::max(h, k));
            }
            else
            {
                bvn = -bvn;
            }
            
            if (k > h)
            {
                bvn += standard_normal_cdf<T>(k) - standard_normal_cdf<T>(h);
            }
        }
    }
    
    return bvn;
}

#endif
