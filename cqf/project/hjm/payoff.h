#ifndef __PAYOFF_H__
#define __PAYOFF_H__

#include <iostream>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

using std::vector;


class payoff
{
    public :
        /* Virtual DTOR for overloading */
        virtual ~payoff() { };
        
        /* Function to build the required mc points */
        virtual void required_rate(std::back_insert_iterator<vector<double> > &r) const = 0;
        virtual void required_date(std::back_insert_iterator<vector<double> > &d) const = 0;

        /* Function to evaluate the payoff after mc simulation */
        virtual double evaluate(const vector<double> &rates, const vector<double> &dates, const double *const mc_points) const = 0;
    
    protected :
        /* Utility to integrate a rate from time t0 to time t1 */
        double integrate_rate(const vector<double> &rates, const vector<double> &dates, const double *const mc_points, const int t0_idx, const int t1_idx, const int rate_idx) const
        {
            /* t1 must be after t0 */
            assert(t0_idx <= t1_idx);
            if (t0_idx == t1_idx)
            {
                return 0.0;
            }
            
            /* Integrate from t0 to t1 */
            double integ = 0.0;
            for (int i = t0_idx; i < t1_idx - 1;)
            {
                const int mc_idx = (i * rates.size()) + rate_idx;
                
                /* Use Simpsons rule if possible */
                if ((i + 1) < t1_idx)
                {
                    const double dd = dates[i + 2] - dates[i];
                    integ += (mc_points[mc_idx] + 4.0 * mc_points[mc_idx + 1] + mc_points[mc_idx + 2]) * dd * (1.0 / 6.0);
                    i += 2;
                }
                /* Else use the trapezium rule */
                else
                { 
                    const double dd = dates[i + 1] - dates[i];
                    integ += (mc_points[mc_idx] + mc_points[mc_idx + 1]) * dd * 0.5;
                    ++i;
                }
            }
            
            return integ;
        }

        /*  Utility to integrate the rate curve from rate r0 to rate r1 on a given day */
        double integrate_curve(const vector<double> &rates, const double *const mc_points, const int r0_idx, const int r1_idx, const int date_idx) const
        {
            /* r1 must be greater than r0 */
            assert(r0_idx <= r1_idx);
            if (r0_idx == r1_idx)
            {
                return 0.0;
            }
            
            /* Integrate from r0 to r1 using the trapezium rule */
            double integ = 0.0;
            for (int i = r0_idx; i < r1_idx - 1;)
            {
                const int mc_idx = (date_idx * rates.size()) + i;
                
                /* Use Simpsons rule if possible */
                if ((i + 1) < r1_idx)
                {
                    const double dr = rates[i + 2] - rates[i];
                    integ += (mc_points[mc_idx] + 4.0 * mc_points[mc_idx + 1] + mc_points[mc_idx + 2]) * dr * (1.0 / 6.0);
                    i += 2;
                }
                /* Else use the trapezium rule */
                else
                {
                    const double dr = rates[i + 1] - rates[i];
                    integ += (mc_points[mc_idx] + mc_points[mc_idx + 1]) * dr * 0.5;
                    ++i;
                }
            }
            
            return integ;
        }
};

#endif
