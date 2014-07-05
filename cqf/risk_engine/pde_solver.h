#ifndef __PDE_SOLVER_H__
#define __PDE_SOLVER_H__

#include <vector>
#include <algorithm>

#include "boost/noncopyable.hpp"


template<class T>
class scenario_point;


template<class T>
class pde_solver : private boost::noncopyable
{
    public :
        /* CTOR */
        
        /* Virtual DTOR */    
        virtual ~pde_solver() {  };
        
        virtual void solve(scenario_point<T> &p) = 0;
};

/*template<class T>
class pde_solver_2d : public pde_solver<T>
{
    public :
        pde_solver_2d(const std::vector<cashflow<T> *> &cashflows, 
            const std::vector<equity_economics<T>*> &economics,
            const std::string &log_file, const T corr, const int nts, 
            const int nas1, const int nas2, const T t, bool log_final_only)
            : pde_solver<T>(cashflows, economics, std::vector<grid<T>*>(), const_ds_boundary_condition<T>(0.0, 300.0), log_file, nts, nas1, t, log_final_only), 
            _corr(corr), _nas2(nas2)
            {
                assert(this->_economics.size() == 2);
            };
            
        virtual ~pde_solver_2d() {  };
        
    protected :
        const T     _corr;
        const int   _nas2;
};*/

#endif
