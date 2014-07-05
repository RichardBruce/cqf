#ifndef __METROPOLIS_H__
#define __METROPOLIS_H__

#include "optimiser.h"

#include "mersenne_twister.h"


class metropolis : public optimiser
{
    public :
        /* CTOR */
        metropolis(fd_solver &solver, const double tolerance = 0.0001, const double anlg_period = 0.5, const int max_iter = 500)
            : optimiser(solver), tolerance(tolerance), anlg_period(anlg_period), max_iter(max_iter)
            {
                rand_gen = new mersenne_twister();
            }
        
        /* DTOR */
        ~metropolis()
        {
            delete rand_gen;
        }
        
        /* Abstract base member */
        virtual double maximise(const std::vector<option*> &hedge, option &target, const double s);
    
    private :
        /* Copy CTOR */
        metropolis(const metropolis &m) 
            : optimiser(m.solver), tolerance(m.tolerance), anlg_period(m.anlg_period), max_iter(m.max_iter)
            {
                rand_gen = NULL;
            }
        
        /* Assignment operator */
        metropolis& operator=(const metropolis &) { return *this; }
        
        double perturb(double *values, std::vector<option*> *nodes, int n, int low_idx, double low_val, double d, const double s, const int iter);
        
        mersenne_twister    *rand_gen;
        const double        tolerance;
        const double        anlg_period;
        const int           max_iter;
};

#endif
