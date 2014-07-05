#ifndef __CALL_OPTION_H__
#define __CALL_OPTION_H__

#include <algorithm>


class call_option : public option
{
    public :
        /* CTOR */
        call_option(const double expiry, const double strike, const double notional, const double bid, const double ask) 
         : option(expiry, strike, notional, bid, ask)
        {
        
        }
        
        /* Base class virtual functions */        
        virtual void payoff(double *p, const double s_inc, const int nas)
        {
            for (int i = 0; i < nas; i++)
            {
                const double s = static_cast<double>(i) * s_inc;
                p[i] = std::max(s - strike, 0.0) * notional;
            }
        }
        
        virtual option* clone()
        {
            return new call_option(expiry, strike, notional, bid, ask);
        }
};

#endif
