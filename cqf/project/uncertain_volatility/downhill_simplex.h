#ifndef __DOWNHILL_SIMPLEX_H__
#define __DOWNHILL_SIMPLEX_H__

#include "optimiser.h"


class downhill_simplex : public optimiser
{
    public :
        downhill_simplex(fd_solver &solver, const double tolerance = 0.0001, const int max_iter = 500)
            : optimiser(solver), tolerance(tolerance), max_iter(max_iter)
            {
        
            }
        
        virtual double maximise(const std::vector<option*> &hedge, option &target, const double s);
    
    private :
        void deform(double *values, std::vector<option*>*nodes, int n, int low_idx, double *low_val, double d, const double s);
        
        const double    tolerance;
        const int       max_iter;
};

#endif
