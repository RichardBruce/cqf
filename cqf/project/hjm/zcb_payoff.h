#ifndef __ZCB_PAYOFF_H__
#define __ZCB_PAYOFF_H__

#include "payoff.h"

/* At the start date integrate the forward curve up to the tenor */
class zcb_payoff : public payoff
{
    public :
        zcb_payoff(const double s, const double t, const double n)
        : start(s), tenor(t), notional(n) { }
        
        void required_rate(std::back_insert_iterator<std::vector<double> > &r) const
        {
            r = 0.0;
        }
        
        void required_date(std::back_insert_iterator<std::vector<double> > &d) const
        {
            d = this->start;
            d = this->start + this->tenor;
        }
        
        double evaluate(const vector<double> &rates, const vector<double> &dates, const double *const mc_points) const
        {
            /* Find start date in the mc points */
            const std::vector<double>::const_iterator start_ptr = std::find(dates.begin(), dates.end(), this->start);
            const int start_idx = start_ptr - dates.begin();
            assert(start_ptr != dates.end());

             /* Find the tenor date in the mc points */
            const std::vector<double>::const_iterator tenor_ptr = std::find(dates.begin(), dates.end(), this->start + this->tenor);
            const int tenor_idx = tenor_ptr - dates.begin();
            assert(tenor_ptr != dates.end());
            
            /* Find the risk free rate */
            const std::vector<double>::const_iterator rate_ptr = std::find(rates.begin(), rates.end(), 0.0);
            const int rate_idx = rate_ptr - rates.begin();
            assert(rate_ptr != rates.end());
            
            /* Discount at the risk free rate the start date */
            const double pv_dr = integrate_rate(rates, dates, mc_points, start_idx, tenor_idx, rate_idx);
            const double pv_df = std::exp(-pv_dr);

            return pv_df * this->notional;
        }

    private :
        const double    start;      /* Time to the start of the ZCB     */
        const double    tenor;      /* Time in the period of the ZCB    */
        const double    notional;   /* The notional of the ZCB          */
};

#endif
