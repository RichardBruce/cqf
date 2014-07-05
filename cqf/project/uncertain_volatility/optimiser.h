#ifndef __OPTIMISER_H__
#define __OPTIMISER_H__

class fd_solver;
class option;

class optimiser
{
    public :
        optimiser(const fd_solver &solver) 
            : solver(solver)
            {
        
            }
            
        /* Virtual DTOR for inheritance */
        virtual ~optimiser() { };
        
        virtual double maximise(const std::vector<option*> &hedge, option &target, const double s) = 0;
    
    protected :
        const fd_solver &solver;
};

#endif
