#ifndef __FD_SOLVER_H__
#define __FD_SOLVER_H__

#include <assert.h>


/* Forward declaration */
class option;

class fd_solver
{
    public :
        fd_solver(const double max_vol, const double min_vol, const double r, const double max_t_inc, const double max_s_inc)
            : max_vol(max_vol), min_vol(min_vol), r(r), max_t_inc(max_t_inc), max_s_inc(max_s_inc)
            {
                assert(max_vol >= min_vol);
            }
            
        double solve(const std::vector<option*> &options, const double s) const;
            
    private :
        const double    max_vol;
        const double    min_vol;
        const double    r;
        const double    max_t_inc;
        const double    max_s_inc;
};

#endif
