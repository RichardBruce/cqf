#ifndef __OPTION_H__
#define __OPTION_H__

#include <iostream>


class option
{
    public :
        /* CTOR */
        option(const double expiry, const double strike, const double notional, const double bid, const double ask) 
            : expiry(expiry), strike(strike), notional(notional), bid(bid), ask(ask)
            {
        
            }
            
        /* Virtual DTOR for inheritance */
        virtual ~option() { };
        
        option& set_notional(const double n)
        {
            notional = n;
            return (*this);
        }
        
        double get_notional() const
        {
            return notional;
        }
        
        double get_expiry() const
        {
            return expiry;
        }
        
        double get_strike() const
        {
            return strike;
        }
        
        double get_price() const
        {
            /* Sell */
            if (notional < 0.0)
            {
                return notional * bid;
            }
            /* Buy */
            else
            {
                return notional * ask;
            }
        }
        
        void dump()
        {
            std::cout << "expiry: " << expiry << std::endl;
            std::cout << "strike: " << strike << std::endl;
            std::cout << "notional: " << notional << std::endl;
            std::cout << "bid: " << bid << std::endl;
            std::cout << "ask: " << ask << std::endl;
        }
        
        /* Get the payoff of an options (option specific) */
        virtual void payoff(double *p, const double s_inc, const int nas) = 0;
        
        /* Allow the option to be copied */
        virtual option* clone() = 0;
    
    protected :
        const double    expiry;
        const double    strike;
        double          notional;
        const double    bid;
        const double    ask;
};


#endif
