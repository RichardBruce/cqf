#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <assert.h>

#include "metropolis.h"

#include "option.h"
#include "fd_solver.h"


double metropolis::perturb(double *values, std::vector<option*> *nodes, int n, int low_idx, double low_val, double d, const double s, const int iter)
{
    const double flt_n          = static_cast<double>(n);
    const double flt_n_p1_inv   = 1.0 / static_cast<double>(n + 1);
    const double int_max_inv    = 1.0 / static_cast<double>(std::numeric_limits<int>::max());

    /* Annealing schedule */
    const double max_iter_flt   = static_cast<double>(max_iter);  
    const double iter_flt       = static_cast<double>(iter);  
    const double anlg_max_iter  = static_cast<double>(max_iter_flt * anlg_period);
    double rand_scale           = 0.0;
    if (iter_flt < anlg_max_iter)
    {
        rand_scale = (anlg_max_iter - iter_flt) / anlg_max_iter;
    }
    
    /* Perturb the nodes of the simplex */
    for (int i = 0; i < n; i++)
    {
        double center = 0.0;
        for (int j = 0; j <= n; j++)
        {
            center += nodes[j][i]->get_notional() * flt_n_p1_inv;
        }

        const double notional = d * nodes[low_idx][i]->get_notional() + (1.0 - d) * ((flt_n + 1.0) / flt_n * center - nodes[low_idx][i]->get_notional() / flt_n);
        const double dist = notional - nodes[low_idx][i]->get_notional();
        const double rand = static_cast<double>(rand_gen->get_next()) * int_max_inv;
        nodes[low_idx][i]->set_notional(notional + (rand * rand_scale * dist));
    }
    
    /* See if the result is better */
    const double trial_val = solver.solve(nodes[low_idx], s);
    if (trial_val > low_val)
    {
        values[low_idx] = trial_val;
        return trial_val;
    }
    
    return low_val;
}


double metropolis::maximise(const std::vector<option*> &hedge, option &target, const double s)
{
    /* Value the initial guess */
//    std::cout << "Beginning Initial Solve" << std::endl;
    double *values = new double [hedge.size() + 1];
    std::vector<option*> *nodes = new std::vector<option*> [hedge.size() + 1];
    for (unsigned int i = 0; i <= hedge.size(); i++)
    {
        nodes[i].resize(hedge.size() + 1);
        for (unsigned int j = 0; j < hedge.size(); j++)
        {
            nodes[i][j] = hedge[j]->clone();
            if (i == j)
            {
                nodes[i][j]->set_notional(nodes[i][j]->get_notional() * 2.0);
            }
        }
        nodes[i][hedge.size()] = target.clone();
        values[i] = solver.solve(nodes[i], s);
    }
//    std::cout << "Completed Initial Solve" << std::endl;
    
//    std::cout << "Beginning Optimisation" << std::endl;   
    int hi_idx;
    int iter = 0;
    double rtol = 2.0 * tolerance;
    while ((rtol > tolerance) && (iter < max_iter))
    {
        /* Find the two lowest and the highest nodes */
        unsigned int low_idx  = 0;
        unsigned int slw_idx  = 0;
        unsigned int hi_idx   = 0;
        double low_val  = values[low_idx];
        double slw_val  = values[slw_idx];
        double hi_val   = values[hi_idx];
        for (unsigned int i = 1; i <= hedge.size(); i++)
        {
            if (values[i] < low_val)
            {
                slw_idx = low_idx;
                slw_val = low_val;
                low_idx = i;
                low_val = values[i];
            }
            else if (values[i] < slw_val)
            {
                slw_idx = i;
                slw_val = values[i];
            }
            else if (values[i] > hi_val)
            {
                hi_idx = i;
                hi_val = values[i];
            }
        }
        
        /* Reflect the lowest value */
        double trial_val = perturb(values, nodes, hedge.size(), low_idx, low_val, -1.0, s, iter);
        
        /* If it worked keep going */
        if (trial_val >= hi_val)
        {
            trial_val = perturb(values, nodes, hedge.size(), low_idx, trial_val, 2.0, s, iter);
        }
        /* Else if its still the lowest try a short distance */
        else if (trial_val <= slw_val)
        {
            const double tmp = values[low_idx];
            trial_val = perturb(values, nodes, hedge.size(), low_idx, trial_val, 0.5, s, iter);
            if (trial_val <= tmp)
            {
                for (unsigned int i = 0; i <= hedge.size(); i++)
                {
                    if (i != hi_idx)
                    {
                        for (unsigned int j = 0; j < hedge.size(); j++)
                        {
                            nodes[i][j]->set_notional(0.5 * (nodes[i][j]->get_notional() + nodes[hi_idx][j]->get_notional()));
                        }
                    }
                    values[i] = solver.solve(nodes[i], s);
                }
            }
        }
        
        ++iter;
        rtol = 2.0 * fabs(hi_val - trial_val) / (fabs(hi_val) + fabs(low_val));
    }

//    std::cout << "Optimisation Complete" << std::endl;
    double max_val = values[0];
    hi_idx = 0;
    for (unsigned int i = 1; i <= hedge.size(); i++)
    {
        if (values[i] > max_val)
        {
            hi_idx = i;
            max_val = values[i];
        }
    }
    
//    std::cout << "Max value is " << values[hi_idx] << " index " << hi_idx << std::endl;
    std::cout << max_val;
    for (unsigned int i = 0; i < hedge.size(); i++)
    {
        std::cout << " " << nodes[hi_idx][i]->get_notional();
    }
    std::cout << std::endl;
    
    /* Clean up */
    for (unsigned int i = 0; i <= hedge.size(); i++)
    {
        for (unsigned int j = 0; j <= hedge.size(); j++)
        {
            delete nodes[i][j];
        }
    }
    delete [] nodes;
    delete [] values;
    
    return max_val;
}
