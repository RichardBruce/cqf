#ifndef __BINARY_PUT_OPTION_H__
#define __BINARY_PUT_OPTION_H__

class binary_put_option : public option
{
    public :
        /* CTOR */
        binary_put_option(const double expiry, const double strike, const double notional, const double bid, const double ask) 
         : option(expiry, strike, notional, bid, ask)
        {
        
        }
        
        /* Base class virtual functions */        
        virtual void payoff(double *p, const double s_inc, const int nas)
        {
            for (int i = 0; i < nas; i++)
            {
                const double s = static_cast<double>(i) * s_inc;
                if (s < strike)
                {
                    p[i] = notional;
                }
                else
                {
                    p[i] = 0.0;
                }
            }
        }
        
        virtual option* clone()
        {
            return new binary_put_option(expiry, strike, notional, bid, ask);
        }
};

#endif
