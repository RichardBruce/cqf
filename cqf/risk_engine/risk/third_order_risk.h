#ifndef __THIRD_ORDER_RISK_H__
#define __THIRD_ORDER_RISK_H__

#include "risk.h"


/* Bumper for third order volatility risk */
template<class T>
class dvol3_scenario_bumper : public scenario_bumper<T>
{
    public :
        dvol3_scenario_bumper(const T order, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD) 
            : scenario_bumper<T>(order, bump_size, bump_dir, bump_style)
            {  };
            
        aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const
        {
            T bumped_vol;
            scenario_point<T> *bumped;
            switch (this->_bump_dir)
            {
                case BI :
                    bumped_vol = this->apply_style_and_bump(null_bump.get_volatility(), (-2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_volatility(bumped_vol);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_vol = this->apply_style_and_bump(null_bump.get_volatility(), -this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_volatility(bumped_vol);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_vol = this->apply_style_and_bump(null_bump.get_volatility(),  this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_volatility(bumped_vol);
                    bumps.push_back(this->insert(scenarios, bumped));

                    bumped_vol = this->apply_style_and_bump(null_bump.get_volatility(), (2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_volatility(bumped_vol);
                    bumps.push_back(this->insert(scenarios, bumped));
                    break;
                
                default :
                    assert(!"Error: Unknown bump direction.");
            }
            
            T* s_bump_sizes = new T[4];
            s_bump_sizes[0] = 0.0;
            s_bump_sizes[1] = 0.0;
            s_bump_sizes[2] = 0.0;
            s_bump_sizes[3] = 0.0;
            
            return new aggregator<T>(*bumps[1], s_bump_sizes, this->get_risk_weightings(100.0), 
                this->_order, bumps[0]->get_grid().total_geometric_size());
        }
};


/* Bumper for third order interest rate risk */
template<class T>
class drate3_scenario_bumper : public scenario_bumper<T>
{
    public :
        drate3_scenario_bumper(const T order, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD) 
            : scenario_bumper<T>(order, bump_size, bump_dir, bump_style)
            {  };
            
        aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const
        {
            T bumped_ir;
            scenario_point<T> *bumped;
            switch (this->_bump_dir)
            {
                case BI :
                    bumped_ir = this->apply_style_and_bump(null_bump.get_interest_rate(), (-2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_interest_rate(bumped_ir);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_ir = this->apply_style_and_bump(null_bump.get_interest_rate(), -this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_interest_rate(bumped_ir);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_ir = this->apply_style_and_bump(null_bump.get_interest_rate(),  this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_interest_rate(bumped_ir);
                    bumps.push_back(this->insert(scenarios, bumped));

                    bumped_ir = this->apply_style_and_bump(null_bump.get_interest_rate(), (2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_interest_rate(bumped_ir);
                    bumps.push_back(this->insert(scenarios, bumped));
                    break;
                
                default :
                    assert(!"Error: Unknown bump direction.");
            }
            
            T* s_bump_sizes = new T[4];
            s_bump_sizes[0] = 0.0;
            s_bump_sizes[1] = 0.0;
            s_bump_sizes[2] = 0.0;
            s_bump_sizes[3] = 0.0;
            
            return new aggregator<T>(*bumps[1], s_bump_sizes, this->get_risk_weightings(100.0), 
                this->_order, bumps[0]->get_grid().total_geometric_size());
        }
};


/* Bumper for third order time risk */
template<class T>
class dtime3_scenario_bumper : public scenario_bumper<T>
{
    public :
        dtime3_scenario_bumper(const T order, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD) 
            : scenario_bumper<T>(order, bump_size, bump_dir, bump_style)
            {  };
            
        aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const
        {
            /* The null bump must be simulated to pick up cashflows */
            scenario_point<T> *null_copy = this->insert(scenarios, new scenario_point<T>(null_bump));
            
            T bumped_time;
            scenario_point<T> *bumped;
            switch (this->_bump_dir)
            {
                case BI :
                    bumped_time = this->apply_style_and_bump(null_bump.get_time(), (-2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_time_offset(null_bump.get_time() - bumped_time, null_copy);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_time = this->apply_style_and_bump(null_bump.get_time(), -this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_time_offset(null_bump.get_time() - bumped_time, null_copy);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_time = this->apply_style_and_bump(null_bump.get_time(),  this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_time_offset(null_bump.get_time() - bumped_time, null_copy);
                    bumps.push_back(this->insert(scenarios, bumped));

                    bumped_time = this->apply_style_and_bump(null_bump.get_time(), (2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_time_offset(null_bump.get_time() - bumped_time, null_copy);
                    bumps.push_back(this->insert(scenarios, bumped));
                    break;
                
                default :
                    assert(!"Error: Unknown bump direction.");
            }
            
            T* s_bump_sizes = new T[4];
            s_bump_sizes[0] = 0.0;
            s_bump_sizes[1] = 0.0;
            s_bump_sizes[2] = 0.0;
            s_bump_sizes[3] = 0.0;
            
            return new aggregator<T>(*null_copy, s_bump_sizes, this->get_risk_weightings(100.0), 
                this->_order, bumps[0]->get_grid().total_geometric_size());
        }
};


/* Bumper for third order spot risk */
template<class T>
class dspot3_scenario_bumper : public scenario_bumper<T>
{
    public :
        dspot3_scenario_bumper(const T order, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD) 
            : scenario_bumper<T>(order, bump_size, bump_dir, bump_style)
            {  };
            
        aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const
        {
            assert((null_bump.get_grid().dimensions() == 1) || (null_bump.get_grid().dimensions() == 2 && null_bump.get_grid().outer_grid_disconnected()));

            T bumped_spot;
            scenario_point<T> *bumped;
            const int dim_m1 = null_bump.get_grid().dimensions() - 1;
            switch (this->_bump_dir)
            {
                case BI :
                    bumped_spot = this->apply_style_and_bump(null_bump.get_spot(), (-2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_spot(bumped_spot, dim_m1);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_spot = this->apply_style_and_bump(null_bump.get_spot(), -this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_spot(bumped_spot, dim_m1);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_spot = this->apply_style_and_bump(null_bump.get_spot(),  this->_bump_size);
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_spot(bumped_spot, dim_m1);
                    bumps.push_back(this->insert(scenarios, bumped));
                    
                    bumped_spot = this->apply_style_and_bump(null_bump.get_spot(),  (2.0 * this->_bump_size));
                    bumped = new scenario_point<T>(null_bump);
                    bumped->set_spot(bumped_spot, dim_m1);
                    bumps.push_back(this->insert(scenarios, bumped));
                    break;
                
                default :
                    assert(!"Error: Unknown bump direction.");
            }
                    
            T* s_bump_sizes = new T[4];
            s_bump_sizes[0] = -this->_bump_size * 2.0;
            s_bump_sizes[1] = -this->_bump_size;
            s_bump_sizes[2] =  this->_bump_size;
            s_bump_sizes[3] =  this->_bump_size * 2.0;
            
            /* Create a node that will be stepped, but not calculated */
            scenario_point<T> *zombie = new scenario_point<T>(null_bump);
            zombie->zombify();
            zombie = this->insert(scenarios, zombie);
            
            return new aggregator<T>(*zombie, s_bump_sizes, this->get_risk_weightings(100.0), 
                this->_order, bumps[0]->get_grid().total_geometric_size());
        }
};

#endif
