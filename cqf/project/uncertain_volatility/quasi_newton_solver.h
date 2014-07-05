#ifndef __QUASI_NEWTON_SOLVER_H__
#define __QUASI_NEWTON_SOLVER_H__

#include "optimiser.h"


class quasi_newton_solver : public optimiser
{
    public :
        quasi_newton_solver(fd_solver &solver, const double tolerance = 0.0001, const int max_iter = 500)
            : optimiser(solver), finite_diff_step(0.00001), alpha(1.0e-4), max_step(10.0), 
            tolerance(tolerance), grad_tol(0.00001), max_iter(max_iter)
            {
        
            }
        
        virtual double maximise(const std::vector<option*> &hedge, option &target, const double s);
    
    private :
        void direction(const std::vector<option*> &p, double *df, const double f, const double s, const int n);
        
        bool line_search(std::vector<option*> &xold, const double fold, double *g, double *p,
	std::vector<option*> &x, double &f, const double stpmax, int n, double s);
        
        const double    finite_diff_step;
        const double    alpha;
        const double    max_step;
        const double    tolerance;
        const double    grad_tol;    
        const int       max_iter;
};

#endif
