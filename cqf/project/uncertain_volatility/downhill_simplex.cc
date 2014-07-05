#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <assert.h>

#include "downhill_simplex.h"

#include "option.h"
#include "fd_solver.h"


void downhill_simplex::deform(double *values, std::vector<option*> *nodes, int n, int low_idx, double *low_val, double d, const double s)
{
    const double flt_n = static_cast<double>(n);
    const double flt_n_p1_inv = 1.0 / static_cast<double>(n + 1);
    for (int i = 0; i < n; i++)
    {
        double center = 0.0;
        for (int j = 0; j <= n; j++)
        {
            center += nodes[j][i]->get_notional() * flt_n_p1_inv;
        }

        const double notional = d * nodes[low_idx][i]->get_notional() + (1.0 - d) * ((flt_n + 1.0) / flt_n * center - nodes[low_idx][i]->get_notional() / flt_n);
        nodes[low_idx][i]->set_notional(notional);
    }
    
    const double trial_val = solver.solve(nodes[low_idx], s);
    if (trial_val > (*low_val))
    {
        values[low_idx] = trial_val;
        (*low_val) = trial_val;
    }
}


double downhill_simplex::maximise(const std::vector<option*> &hedge, option &target, const double s)
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
        unsigned int low_idx = 0;
        unsigned int slw_idx = 0;
        unsigned int hi_idx  = 0;
        double  low_val = values[low_idx];
        double  slw_val = values[slw_idx];
        double  hi_val  = values[hi_idx];
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
        
        deform(values, nodes, hedge.size(), low_idx, &low_val, -1.0, s);
        double trial_val = low_val;
        if (trial_val >= hi_val)
        {
            deform(values, nodes, hedge.size(), low_idx, &low_val, 2.0, s);
            trial_val = low_val;
        }
        else if (trial_val <= slw_val)
        {
            const double tmp = values[low_idx];
            deform(values, nodes, hedge.size(), low_idx, &low_val, 0.5, s);
            trial_val = low_val;
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
        rtol = 2.0 * fabs(hi_val - low_val) / (fabs(hi_val) + fabs(low_val));
    }

////    std::cout << "Optimisation Complete" << std::endl;
////    double* res = new double [hedge.size() * (hedge.size() + 2)];
//    for(int i = 0; i <= hedge.size(); i++)
//    {
//        for (int j = 0; j < hedge.size(); j++)
//        {
////            res[(i * (hedge.size() + 1)) + j] = nodes[i][j]->get_notional();
//            std::cout << "Notional " << j << " = " << nodes[i][j]->get_notional() << ", ";
//        }
////        res[(i * (hedge.size() + 1)) + hedge.size()] = values[i];
//        std::cout << "Value " << values[i] << std::endl;
//    }
    
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
    std::cout << -max_val;
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
