#ifndef __IMPLIED_VOLATILITY_H__
#define __IMPLIED_VOLATILITY_H__

#include "european_call.h"
#include "european_call_put.h"


class implied_volatility
{
    public :
        implied_volatility(const double target_error, const int iteration_limit) : 
        _target_error(target_error), _iteration_limit(iteration_limit) {  }

        /* Functions to solve for the volatility of individual options */
        double back_out_from_call(const double s, const double k, const double r, const double t, const double v, double guess = 0.2)
        {
            /* Newton Raphson solver for volatility */
            int i = 0;
            double err = _target_error + 1.0f;
            while ((i++ < _iteration_limit) && (err > _target_error))
            {
                const double guess_v = european_call_value(s, k, guess, r, t);
                const double vega = european_call_put_vega(s, k, guess, r, t);
                
                err = guess_v - v;
                guess = guess - (err / vega);
            }
            
            return guess;
        }
        
    private :
        const double _target_error;      /* Target error when solving for implied vol                    */
        const int    _iteration_limit;   /* Maximum number of iterations when solving for implied vol    */
};

#endif
