#ifndef __HJM_SIMULATOR_H__
#define __HJM_SIMULATOR_H__

#include <string>
#include <vector>


#include "blocked_range.h"
#include "atomic.h"

#include "mersenne_twister.h"
#include "box_muller.h"
#include "sobol_numbers.h"
#include "standard_normal.h"
#include "brownian_bridge.h"

#include "payoff.h"


class hjm_simulator
{
    public :
        /* CTOR */
        hjm_simulator(const std::string &dir_nums, const std::vector<payoff *> &payoffs, 
        double *init_curve, const std::vector<double> &rates, 
        const std::vector<double> &dates, const double *const vols, 
        const double *const drifts, const int total_sims, const int factors) :
          payoffs(payoffs), rates(rates), dates(dates), vols(vols), drifts(drifts), 
          nr_of_paths_inv(1.0 / static_cast<double>(total_sims)), average_value(0.0), 
          sobol_rand(NULL), bridge(NULL), factors(factors)
        {
            /* Take a copy of the initial curve */
            this->fwd_curve = new double [rates.size() * dates.size()];
            memcpy(this->fwd_curve, init_curve, rates.size() * sizeof(double));

            /* Build the corect random number generator */
            /* This is pretty hacky, but theres no way im paying for a virtual function
               call when the code should be inlined (sobol: all ~1 XOR worth) */
            if (!dir_nums.empty())
            {
                /* Build sobol number generator */
                sobol_rand = new sobol_numbers<double>(dir_nums.c_str(), total_sims + 1, ((dates.size() - 1) * factors));
                
                /* Build Brownian Bridge */
                bridge = new brownian_bridge<double>(dates.size() - 1);
            }
            else
            {
                /* Build normally distributed random number generators */
                const int idx = hjm_simulator::seed_idx.fetch_and_add(1);
                mt_rand = new mersenne_twister(5489 * (idx + 1));
            }
        }
        
        
        /* TBB splitting CTOR */
        hjm_simulator(const hjm_simulator &h, tbb::split) :
          payoffs(h.payoffs), rates(h.rates), dates(h.dates), vols(h.vols), drifts(h.drifts), 
          nr_of_paths_inv(h.nr_of_paths_inv), average_value(0.0), factors(h.factors)
        {
            /* Take a copy of the initial curve */
            this->fwd_curve = new double [rates.size() * dates.size()];
            memcpy(this->fwd_curve, &h.fwd_curve[0], rates.size() * sizeof(double));
            
            /* Hacky see above */
            if (h.sobol_rand != NULL)
            {
                sobol_rand  = new sobol_numbers<double>((*h.sobol_rand));
                bridge      = new brownian_bridge<double>((*h.bridge));
            }
            else
            {
                const int idx = hjm_simulator::seed_idx.fetch_and_add(1);
                mt_rand = new mersenne_twister(5489 * (idx + 1));
            }
        }

        
        /* DTOR */
        ~hjm_simulator()
        {
            delete [] fwd_curve;
            
            /* Hacky again see above */
            if (sobol_rand != NULL)
            {
                delete sobol_rand;
                delete bridge;
            }
            else
            {
                delete mt_rand;
            }
        }
        
        void operator() (const tbb::blocked_range<size_t> &r)
        {
            /* Precalculate the time steps */
            std::vector<double> dt(dates.size());
            std::adjacent_difference(dates.begin(), dates.end(), dt.begin());
            std::vector<double> sqrt_dt(dt);
            for (unsigned int i = 0; i < sqrt_dt.size(); i++)
            {
                sqrt_dt[i] = std::sqrt(sqrt_dt[i]);
            }


            /* Precalculate the rate steps */
            std::vector<double> dr_inv(rates.size());
            std::adjacent_difference(rates.begin(), rates.end(), dr_inv.begin());
            for (unsigned int i = 0; i < dr_inv.size(); i++)
            {
                dr_inv[i] = 1.0 / dr_inv[i];
            }

            /* Prepare the random number generator */        
            double *gen_nums    = new double [(dt.size() - 1) * factors];
            double *normal_nums = NULL;
            if (sobol_rand != NULL)
            {
                normal_nums = new double [(dt.size() - 1) * factors];
                
                /* Skip to the numbers at the start of the sim */
                if (r.begin() != 0)
                {
                    sobol_rand->skip(r.begin() - 1);
                }
                else
                {
                    /* Skipped because 0 doesnt do well in an inverse cumulative normal */
                    sobol_rand->get_next(gen_nums);
                }
            }
            
            /*  Run lots of sims */
            for (size_t sim = r.begin(); sim < r.end(); sim++)
            {
                /* Random numbers for this iteration */
                if (sobol_rand != NULL)
                {
                    sobol_rand->get_next(&gen_nums[0]);
                    for (unsigned int i = 0; i < ((dt.size() - 1) * factors); i++)
                    {
                        normal_nums[i] = inverse_standard_normal_cdf(gen_nums[i]);
                    }
                    
                    for (int j = 0; j < factors; j++)
                    {
                        bridge->map_randoms(&normal_nums[j * (dt.size() - 1)], &gen_nums[j * (dt.size() - 1)]);
                    }
                }
                else
                {
                    for (unsigned int i = 0; i < ((dt.size() - 1) * factors); i++)
                    {
                        gen_nums[i] = inverse_standard_normal_cdf(mt_rand->get_next_prob());
                    }
                }

                /* Foreach time step */
                unsigned int t_idx = dr_inv.size();
                unsigned int last_t_idx = 0;
                for (unsigned int i = 1; i < dt.size(); i++)
                {
                    /* Build the forward curve a point at a time */
                    unsigned int rand_idx = (i - 1) * factors;
                    double musiela = 0.0;
                    for (unsigned int j = 0; j < dr_inv.size() - 1; j++)
                    {
                        const double last_t = fwd_curve[last_t_idx++];
                        
                        double diff = gen_nums[rand_idx] * vols[j];
                        for (int k = 1; k < factors; k++)
                        {
                            diff += gen_nums[i - 1 + (k * (dt.size() - 1))] * vols[(k * dr_inv.size()) + j];
                        }
                        diff *= sqrt_dt[i];
                    
                        musiela = ((fwd_curve[last_t_idx] - last_t) * dr_inv[j + 1]);
                        fwd_curve[t_idx++] = last_t + diff + (drifts[j] + musiela) * dt[i];
                    }
        
                    /* Last point is built slightly differently */
                    const double last_t = fwd_curve[last_t_idx++];
                
                    double diff = gen_nums[rand_idx] * vols[dr_inv.size() - 1];
                    for (int j = 1; j < factors; j++)
                    {
                        diff += gen_nums[i - 1 + (j * (dt.size() - 1))] * vols[(j * dr_inv.size()) + dr_inv.size() - 1];
                    }
                    diff *= sqrt_dt[i];
                
                    fwd_curve[t_idx++] = last_t + diff + (drifts[dr_inv.size() - 1] + musiela) * dt[i];
                }

                /* Calculate the total payoff of this sim */
                double v = 0.0;
                for (unsigned int j = 0; j < payoffs.size(); j++)
                {
                    const double payoff_comp = payoffs[j]->evaluate(rates, dates, fwd_curve);
                    v += payoff_comp;
                }
                average_value += (v * nr_of_paths_inv);
            }

            /* Clean up */
            delete [] gen_nums;
            
            if (sobol_rand != NULL)
            {
                delete [] normal_nums;
            }
        }
        
        void join(const hjm_simulator &h)
        {
            average_value += h.average_value;
        }
        
        double value() const
        {
            return average_value;
        }
        
    private :
        /* This CTOR was never tested dont use it */
        hjm_simulator(const hjm_simulator &h) :
          payoffs(h.payoffs), rates(h.rates), dates(h.dates), vols(h.vols), drifts(h.drifts), 
          nr_of_paths_inv(h.nr_of_paths_inv), average_value(0.0), mt_rand(h.mt_rand), factors(h.factors)
          {
            if (h.sobol_rand != NULL)
            {
                sobol_rand  = new sobol_numbers<double>(*(h.sobol_rand));
                bridge      = new brownian_bridge<double>(*(h.bridge));
            }
          }
        
        /* This operator was never tested dont use it */
        hjm_simulator& operator=(const hjm_simulator &)
        {
            return *this;
        }
        
        const std::vector<payoff *>         &payoffs;
        const std::vector<double>           &rates;
        const std::vector<double>           &dates; 
        const double                 *const vols;
        const double                 *const drifts;
        const double                        nr_of_paths_inv;
        double                       *      fwd_curve;
        double                              average_value;
        mersenne_twister             *      mt_rand;
        sobol_numbers<double>        *      sobol_rand;
        brownian_bridge<double>      *      bridge;
        const int                           factors;
        static tbb::atomic<int>             seed_idx;
};

tbb::atomic<int> hjm_simulator::seed_idx;   /* Initialised 0 */


#endif
