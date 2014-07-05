#ifndef __KO_FLOORLET_PAYOFF_H__
#define __KO_FLOORLET_PAYOFF_H__

#include "payoff.h"

class ko_floorlet_payoff : public payoff
{
    public :
        ko_floorlet_payoff(const double d, const double e, const double p, const double t, const double n, 
            const double s, const double r, const double b, const double c = 12.0) 
            : start_date(d), eval_date(e), pay_date(p), tenor(t), notional(n), strike(s), rate(r), barrier(b), comp_freq(c) {  }

        void required_rate(std::back_insert_iterator<std::vector<double> > &r) const
        {
            r = 0.0;
            r = this->rate;
        }
        
        void required_date(std::back_insert_iterator<std::vector<double> > &d) const
        {
            d = this->start_date;
            d = this->eval_date;
            d = this->pay_date;
        }

        double evaluate(const vector<double> &rates, const vector<double> &dates, const double *const mc_points) const
        {
            /* Find the date to discount from in the mc points */
            const std::vector<double>::const_iterator dsc_ptr = std::find(dates.begin(), dates.end(), this->pay_date);
            const int dsc_idx = dsc_ptr - dates.begin();
            assert(dsc_ptr != dates.end());
            
            /* Find the evaluation date in the mc points */
            const std::vector<double>::const_iterator tenor_ptr = std::find(dates.begin(), dates.end(), this->eval_date);
            const int tenor_idx = tenor_ptr - dates.begin();
            assert(tenor_ptr != dates.end());
            
            /* Find the start date in the mc points */
            const std::vector<double>::const_iterator start_ptr = std::find(dates.begin(), dates.end(), this->start_date);
            const int start_idx = start_ptr - dates.begin();
            assert(start_ptr != dates.end());
            
            /* Find the risk free rate */
            const std::vector<double>::const_iterator rf_rate_ptr = std::find(rates.begin(), rates.end(), 0.0);
            const int rf_rate_idx = rf_rate_ptr - rates.begin();
            assert(rf_rate_ptr != rates.end());
            
            /* Find the strike rate */            
            const std::vector<double>::const_iterator rate_ptr = std::find(rates.begin(), rates.end(), this->rate);
            const int rate_idx = rate_ptr - rates.begin();
            assert(rate_ptr != rates.end());

            /* Find the rate in the mc points */
            const double mc_rate = integrate_curve(rates, mc_points, rf_rate_idx, rate_idx, tenor_idx) / this->rate;
            
            /* Convert to a discrete rate */
            const double discrete_rate = this->comp_freq * (std::exp(mc_rate / this->comp_freq) - 1.0);
            
            /* Calculate the payoff */
            const double payoff = (mc_rate < barrier) ? 0.0 : (this->tenor * std::max(this->strike - discrete_rate, 0.0));
            
            /* Discount the pay off to the start date at the spot rate */
            const double pv_dr = integrate_rate(rates, dates, mc_points, start_idx, dsc_idx, rf_rate_idx);
            const double pv_df = std::exp(-pv_dr);

            return payoff * pv_df * this->notional;
        }

    private :
        const double    start_date; /* The date to present value to         */
        const double    eval_date;  /* The date to evaluate the payoff      */
        const double    pay_date;   /* The date to pay the payoff           */
        const double    tenor;      /* The period the interest is paid over */
        const double    notional;   /* The notional of the floorlet         */
        const double    strike;     /* Strike of the floorlet               */
        const double    rate;       /* The rate the floorlet is struck on   */
        const double    barrier;    /* KO barrier level                     */
        const double    comp_freq;  /* Compounding frequency                */
};

#endif
