#ifndef __RISK_H__
#define __RISK_H__

#include <iostream>
#include <fstream>
#include <set>
#include <string>

#include "dumping.h"

#include "grid_mapper.h"
#include "scenario_point.h"
#include "aggregator.h"


/* Forward delcaration for typedefs */
template<class T>
class risk;


/* Enums to describe methods of bumping */
enum bump_dir_t { NO_BUMP_DIR = 0, DOWN = 1, UP = 2, BI = 3 };
enum bump_style_t { NO_BUMP_STYLE = 0, ADD = 1, MULT = 2, OVERRIDE = 3 };


/* Virtual base class for scenario bumping */
template<class T>
class scenario_bumper
{
    public : 
        /* CTOR */
        scenario_bumper(const int order, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD) 
            : _order(order), _bump_size(bump_size), _bump_dir(bump_dir), 
              _bump_style(bump_style)
              {  };
        
	/* Virtual DTOR */
	virtual ~scenario_bumper() {  };

        /* Allow default Copy CTOR */
        
        /* Multiplier when aggregating risk */
        T get_risk_multiplier(const T v) const
        {
            const T bump_size = apply_style_and_bump(v, bump_size_with_direction()) - v;
            return pow(bump_size, static_cast<T>(_order));
        }
        
        /* Virtual member for bumping scenarios */
        virtual aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const = 0;

    protected :        
        const int           _order;
        const T             _bump_size;
        const bump_dir_t    _bump_dir;
        const bump_style_t  _bump_style;
        
        
        T apply_style_and_bump(const T data, const T bump) const 
        {
            /* Bump based on style */
            switch (this->_bump_style)
            {
                case NO_BUMP_STYLE :
                    return data;
                case ADD :
                    return data + bump;
                case MULT :
                    return data + (bump * data);
                case OVERRIDE :
                    return bump;
                default :
                    assert(!"Error: Invalid bump style.");
            }
        }

        
        /* Weightings of each scenario point when aggreating risk */
        T* get_risk_weightings(const T v) const
        {
            T *weightings = new T[_order + 1];
            const T mult_inv = 1.0 / get_risk_multiplier(v);
            switch (_order)
            {
                case 0 :
                    weightings[0] = 1.0;
                    break;
                case 1 :
                    weightings[0] = -mult_inv;
                    weightings[1] =  mult_inv;
                    break;
                case 2 :
                    weightings[0] = mult_inv;
                    weightings[1] = mult_inv * -2.0;
                    weightings[2] = mult_inv;
                    break;
                case 3 :
                    weightings[0] =  mult_inv;
                    weightings[1] =  mult_inv * -2.0;
                    weightings[2] =  mult_inv *  2.0;
                    weightings[3] = -mult_inv;
                    break;
            }
            
            return weightings;
        }
        
        
        /* Clean insert to scenario set */
        scenario_point<T>* insert(typename risk<T>::scenario_set &scenarios, scenario_point<T> *bump) const
        {
            auto res = scenarios.insert(bump);
            
            /* Delete the bump if it wasnt inserted */
            if (!res.second)
            {
                (*res.first)->dezombify(*bump);
                delete bump;
            }
            
            /* Return the bump in the set */
            return *res.first;
        }

        
    private :
        T bump_size_with_direction() const
        {
            if (_order == 1)
            {
                if (_bump_dir == DOWN)
                {
                    return -_bump_size;
                }
                else if (_bump_dir == UP)
                {
                    return _bump_size;
                }
                else
                {
                    return 2.0 * _bump_size;
                }
            }
            else if (_order > 1)
            {
                assert(_bump_dir == BI);
                return _bump_size;
            }
            
            return 0.0;
        }
};


/* Bumper for fair value */
template<class T>
class value_scenario_bumper : public scenario_bumper<T>
{
    public :
        value_scenario_bumper() 
            : scenario_bumper<T>(0, 0.0, NO_BUMP_DIR, NO_BUMP_STYLE) {  };
            
        aggregator<T>* bump(std::vector<scenario_point<T>*> &bumps, typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump) const
        {
            bumps.push_back(this->insert(scenarios, new scenario_point<T>(null_bump)));
            
            T* s_bump_sizes = new T[1];
            s_bump_sizes[0] = 0.0;
            
            return new aggregator<T>(*bumps[0], s_bump_sizes, this->get_risk_weightings(100.0), 
                this->_order,  bumps[0]->get_grid().total_geometric_size());
        }
};

    
/* Virtual base class for risk calculations */
template<class T>
class risk
{
    public :
        /* Typedef to reduce typing */
        typedef typename std::set<scenario_point<T> *, scenario_point_compare<T> >                          scenario_set;
        typedef typename std::set<scenario_point<T> *, scenario_point_compare<T> >::iterator                scenario_iter;
        typedef typename std::set<scenario_point<T> *, scenario_point_compare<T> >::const_iterator          scenario_const_iter;
        
        typedef std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> s_range_pair;

        /* CTOR */
        risk(const scenario_bumper<T> *bumper, const std::string &log_file, 
             const T disc_grid_offset, const bool root_node, const bool grid_risk)
            : _log_file(log_file), _bumper(bumper), _aggregator(nullptr), _results(nullptr), 
              _name(log_file), _disc_grid_offset(disc_grid_offset), _root_node(root_node), 
              _grid_risk(grid_risk)
            {  };
           
        /* Copy CTOR */
        risk(const risk<T> &r)
            : _bumper(r._bumper), _aggregator(nullptr), _results(nullptr), 
              _name(r._name), _disc_grid_offset(r._disc_grid_offset), 
              _root_node(r._root_node), _grid_risk(r._grid_risk)
            {
                /* It would probably be ok if root nodes were copied, but not expected */
                assert(!_root_node || !"Error: Root risk node should not be copied");
            };
            
        /* DTOR */
        virtual ~risk()
        {
            if (_log_file.is_open())
            {
                _log_file.close();
            }

            if (_aggregator != nullptr)
            {
                delete _aggregator;
            }
            
            if (_results != nullptr)
            {
                delete [] _results;
            }
        }
        
        
        /* Virtual members */
        virtual void bump(scenario_set &scenarios, scenario_point<T> &null_bump) = 0;
        virtual risk_result<T> aggregate() = 0;
        virtual risk<T>* deep_copy() = 0;
        
    protected :
        std::ofstream                                       _log_file;
        const boost::shared_ptr<const scenario_bumper<T> >  _bumper;
        aggregator<T>               *                       _aggregator;
        risk_result<T>              *                       _results;
        std::vector<scenario_point<T>*>                     _scenario_points;
        const std::string                                   _name;              /* For output   */
        const T                                             _disc_grid_offset;
        const bool                                          _root_node;
        const bool                                          _grid_risk;
        
    private :
        risk<T> operator=(const risk<T> &) {  };
};

    
/* Virtual base class for risk calculations */
template<class T>
class risk_node : public risk<T>
{
    public :
        /* CTOR */
        risk_node(const scenario_bumper<T> *bumper, const std::string &log_file, 
                  const T disc_grid_offset, const bool root_node, const bool grid_risk)
            : risk<T>(bumper, log_file, disc_grid_offset, root_node, grid_risk)
            {  };
           
        /* Copy CTOR */
        risk_node(const risk_node<T> &r) : risk<T>(r) {  };
        
        /* Default DTOR */
        
        /* Virtual members */
        risk_node<T>* deep_copy()
        {
            risk_node<T> *copy = new risk_node<T>(*this);
            return copy; 
        }
        
        
        void bump(typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump)
        {
            /* Create the bumped scenarios */
            this->_aggregator = this->_bumper->bump(this->_scenario_points, scenarios, null_bump);
            this->_results = new risk_result<T>[this->_scenario_points.size()];
            
            return;
        }

        risk_result<T> aggregate()
        {
            /* Aggregation */
            risk_result<T> result;
            if (this->_grid_risk)
            {
                /* Check if all scenarios have begun */
                bool begun = true;
                for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                {
                    begun &= this->_scenario_points[i]->begun();
                }
                
                if (begun)
                {
                    for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                    {
                        this->_results[i] = this->_scenario_points[i]->get_risk_result();
                    }
                    
                    result = this->_aggregator->aggregate(this->_results, false);
                    assert((result.dimensions() == 1) || ((result.dimensions() == 2) && result.outer_grid_disconnected()));
                    if (this->_root_node)
                    {
                        /* Dump middle disconnected grid */
                        const int dim_m1        = result.dimensions() - 1;
                        const int nr_of_grids   = (dim_m1 == 0) ? 0 : static_cast<int>((result.size(0) - 1) * this->_disc_grid_offset);
                        const int grid_offset   = nr_of_grids * result.size(dim_m1);
                        dump_row_to_gnuplot(this->_log_file, &result.values()[grid_offset], result.spot_grid(dim_m1).first, result.spot_grid(dim_m1).second, this->_scenario_points[0]->get_time());
                    }
                }
            }
            else
            {
                /* Check if all scenarios are complete */
                bool done = true;
                for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                {
                    done &= this->_scenario_points[i]->done();
                }
                
                if (done)
                {
                    for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                    {
                        this->_results[i] = this->_scenario_points[i]->get_risk_result();
                    }
                    
                    result = this->_aggregator->aggregate(this->_results, true);
                    if (this->_root_node)
                    {
                        assert((result.dimensions() == 1) || (result.dimensions() == 2 && result.outer_grid_disconnected()));
                        assert((result.size(0) == 1) || (result.dimensions() == 1) || !"Error: Dimensionality must be reduced to 1 before completion.");
                        std::cout << this->_name << ": " << result.values()[0] << std::endl;
                    }
                }
            }
                    
            return result;
        }
        
    private :
        risk_node<T> operator=(const risk<T> &) {  };
};


/* Cross product risk */
template<class T>
class composite_risk : public risk<T>
{
    public :
        /* CTOR */
        composite_risk(const scenario_bumper<T> *const bumper, const std::string &log_file,
                       risk<T> *const risk_node, const T disc_grid_offset, 
                       const bool root_node, const bool grid_risk)
            : risk<T>(bumper, log_file, disc_grid_offset, root_node, grid_risk), 
              _risk_node(new std::vector<risk<T>*>(1, risk_node))
              {  };
            
        /* Copy CTOR */
        composite_risk(const composite_risk<T> &r)
            : risk<T>(r),
              _risk_node(new std::vector<risk<T>*>(1, (*r._risk_node)[0]))
        {  };
            
        /* DTOR */
        ~composite_risk()
        {
            for (unsigned int i = 0; i < _risk_node->size(); i++)
            {
                delete (*_risk_node)[i];
            }
            
            delete _risk_node;
        }

        composite_risk<T>* deep_copy()
        {
            composite_risk<T> *copy = new composite_risk<T>(*this);
            return copy;
        }
        

        void bump(typename risk<T>::scenario_set &scenarios, scenario_point<T> &null_bump)
        {
            /* Create the bumped scenarios */
            this->_aggregator = this->_bumper->bump(this->_scenario_points, scenarios, null_bump);
            
            /* Create sub node and bump them */
            _risk_node->resize(this->_scenario_points.size());
            for (unsigned int i = 1; i < this->_scenario_points.size(); i++)
            {
                (*_risk_node)[i] = (*_risk_node)[0]->deep_copy();
            }
            
            for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
            {
                (*_risk_node)[i]->bump(scenarios, (*this->_scenario_points[i]));
            }
            
            this->_results = new risk_result<T>[this->_scenario_points.size()];
            
            return;
        }
        
        risk_result<T> aggregate() 
        {
            risk_result<T> result;
            if (this->_grid_risk)
            {
                /* Check if all scenarios have begun */
                bool begun = true;
                for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                {
                    begun &= this->_scenario_points[i]->begun();
                }
                
                if (begun)
                {
                    for (unsigned int i = 0; i < _risk_node->size(); i++)
                    {
                        this->_results[i] = (*_risk_node)[i]->aggregate();
                    }
            
                    result = this->_aggregator->aggregate(this->_results, false);
                    assert((result.dimensions() == 1) || (result.dimensions() == 2 && result.outer_grid_disconnected()));
                    if (this->_root_node)
                    {
                        /* Dump middle disconnected grid */
                        const int dim_m1        = result.dimensions() - 1;
                        const int nr_of_grids   = (dim_m1 == 0) ? 0 : static_cast<int>((result.size(0) - 1) * this->_disc_grid_offset);
                        const int grid_offset   = nr_of_grids * result.size(dim_m1);
                        dump_row_to_gnuplot(this->_log_file, &result.values()[grid_offset], result.spot_grid(dim_m1).first, result.spot_grid(dim_m1).second, this->_scenario_points[0]->get_time());
                    }
                }
            }
            else
            {
                /* Check if all scenarios are complete */
                bool done = true;
                for (unsigned int i = 0; i < this->_scenario_points.size(); i++)
                {
                    done &= this->_scenario_points[i]->done();
                }
                
                if (done)
                {
                    for (unsigned int i = 0; i < _risk_node->size(); i++)
                    {
                        this->_results[i] = (*_risk_node)[i]->aggregate();
                    }
            
                    result = this->_aggregator->aggregate(this->_results, true);
                    if (this->_root_node)
                    {
                        assert((result.dimensions() == 1) || (result.dimensions() == 2 && result.outer_grid_disconnected()));
                        assert((result.size(0) == 1) || (result.dimensions() == 1) || !"Error: Dimensionality must be reduced to 1 before completion.");
                        std::cout << this->_name << ": " << result.values()[0] << std::endl;
                    }
                }
            }
            
            return result;
        }
        
    private :
        composite_risk<T> operator=(const composite_risk<T> &) {  };
        
        std::vector<risk<T>*>       *const  _risk_node;
};


#endif
