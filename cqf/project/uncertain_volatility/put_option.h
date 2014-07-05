#ifndef __PUT_OPTION_H__
#define __PUT_OPTION_H__

#include <algorithm>
#include <iostream>


class put_option : public option
{
    public :
        /* CTOR */
        put_option(const double expiry, const double strike, const double notional, const double bid, const double ask) 
         : option(expiry, strike, notional, bid, ask)
        {
        
        }

        /* Base class virtual functions */        
        virtual void payoff(double *p, const double s_inc, const int nas)
        {
            for (int i = 0; i < nas; i++)
            {
                const double s = static_cast<double>(i) * s_inc;
                p[i] = std::max(strike - s, 0.0) * notional;
            }
        }
        
        virtual option* clone()
        {
            return new put_option(expiry, strike, notional, bid, ask);
        }
};

#endif
