#ifndef __OPTION_PAYOFF_H__
#define __OPTION_PAYOFF_H__

#include "payoff.h"

class option_payoff : public payoff
{
    public :
        /* Note the underlyings must all have the same start date which is also the 
           options evaluation date */
        option_payoff(const vector<payoff*> *u, const double d, const double p, const double s, const bool c = true) 
            : underlyings(u), start_date(d), pay_date(p), strike(s), call(c) {  }
        
        ~option_payoff()
        {
            /* Clean up */
            delete underlyings;
        }
        
        void required_rate(std::back_insert_iterator<vector<double> > &r) const
        {
            /* The options required rates */
            r = 0.0;
            
            /* The underlyings required rates */
            for (unsigned int i = 0; i < underlyings->size(); i++)
            {
                (*underlyings)[i]->required_rate(r);
            }
        }
        
        void required_date(std::back_insert_iterator<vector<double> > &d) const
        {
            /* The options required dates */
            d = this->start_date;

            /* The underlyings required dates */
            for (unsigned int i = 0; i < underlyings->size(); i++)
            {
                (*underlyings)[i]->required_date(d);
            }
        }

        double evaluate(const vector<double> &rates, const vector<double> &dates, const double *const mc_points) const
        {
            /* Value the underlyings */
            double undl_val = 0.0;
            for (unsigned int i = 0; i < underlyings->size(); i++)
            {
                undl_val += (*underlyings)[i]->evaluate(rates, dates, mc_points);
            }
            
            /* Find the date to discount from in the mc points */
            const std::vector<double>::const_iterator dsc_ptr = std::find(dates.begin(), dates.end(), this->pay_date);
            const int dsc_idx = dsc_ptr - dates.begin();
            assert(dsc_ptr != dates.end());
            
            /* Find the start date in the mc points */
            const std::vector<double>::const_iterator start_ptr = std::find(dates.begin(), dates.end(), this->start_date);
            const int start_idx = start_ptr - dates.begin();
            assert(start_ptr != dates.end());
            
            /* Find the risk free rate */
            const std::vector<double>::const_iterator rf_rate_ptr = std::find(rates.begin(), rates.end(), 0.0);
            const int rf_rate_idx = rf_rate_ptr - rates.begin();
            assert(rf_rate_ptr != rates.end());
            
            /* Evaluate the option payoff */
            const double value  = call ? (undl_val - this->strike) : (this->strike - undl_val);
            const double payoff = std::max(value, 0.0);
            
            /* Discount the pay off to the start date at the spot rate */
            const double pv_dr = integrate_rate(rates, dates, mc_points, start_idx, dsc_idx, rf_rate_idx);
            const double pv_df = std::exp(-pv_dr);
            
            return payoff * pv_df;
        }
  
    private :
      const vector<payoff*>    *underlyings;   /* The options underlyings      */
      const double              start_date;     /* The date to present value to */
      const double              pay_date;       /* The payment date             */
      const double              strike;         /* the strike                   */
      const bool                call;           /* Call or put option           */
};

#endif
