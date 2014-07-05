#ifndef __RISK_ENGINE_H__
#define __RISK_ENGINE_H__

#include <vector>
#include <set>

#include "boost/noncopyable.hpp"

#include "risk.h"


template<class T>
class risk_engine : private boost::noncopyable
{
    public :
        /* CTOR */
        risk_engine(std::vector<risk<T> *> risks, scenario_point<T> &null_bump)
            : _risks(risks), _scenarios(build_scenarios(null_bump)) {  }
        
        /* DTOR */
        ~risk_engine()
        {
            for (typename risk<T>::scenario_iter i = _scenarios->begin(); i != _scenarios->end(); i++)
            {
                delete (*i);
            }
                
            delete _scenarios;
        }
        
        /* Calculate risks */
        void run(pde_solver<T> &calc_engine)
        {
            /* Initialise the scenarios (build grids, set up initial condition, possible partial time step) */
            for (typename risk<T>::scenario_iter i = _scenarios->begin(); i != _scenarios->end(); i++)
            {
                (*i)->initialise(calc_engine);
            }
            
            for (unsigned int i = 0; i < this->_risks.size(); i++)
            {
                this->_risks[i]->aggregate();
            }
            
            /* Work back in time until done */
            while (!(*_scenarios->begin())->done())
            {
                /* Iterate each scenario point */
                /* This must be done in the set order defined by the scenarios */
                for (typename risk<T>::scenario_iter i = _scenarios->begin(); i != _scenarios->end(); i++)
                {
                    (*i)->iteration(calc_engine);
                }
                
                /* Aggregate all the risk */
                for (unsigned int i = 0; i < this->_risks.size(); i++)
                {
                    this->_risks[i]->aggregate();
                }
            }
        }
    
    
    private :
        const std::vector<risk<T> *>            _risks;
        const typename risk<T>::scenario_set *  _scenarios;
        
        /* Private member to build scenario points required by risks being calculated */
        typename risk<T>::scenario_set* build_scenarios(scenario_point<T> &null_bump)
        {
            typename risk<T>::scenario_set *scenarios = new typename risk<T>::scenario_set();
            for (unsigned int i = 0; i < this->_risks.size(); i++)
            {
                this->_risks[i]->bump(*scenarios, null_bump);
            }
            
            return scenarios;
        }
};

#endif
