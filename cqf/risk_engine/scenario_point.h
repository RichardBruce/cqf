#ifndef __SCENARIO_POINT_H__
#define __SCENARIO_POINT_H__

#include <vector>
#include <algorithm>

#include "boost/shared_ptr.hpp"

#include "utility.h"

#include "grid_mapper.h"

#include "cashflow.h"

#include "pde_solver.h"


template<class S>
class scenario_point_compare
{
    public : 
        bool operator()(const scenario_point<S> *const lhs, const scenario_point<S> *const rhs) const
        {
            return (*lhs) < (*rhs);
        }
};


/* Immutable class to provide a convenient bundle of data needed for risk aggregation */
template<class T>
class risk_result
{
    public :
        typedef std::pair<typename std::vector<T>::const_iterator,typename std::vector<T>::const_iterator> iterator_pair;
        
        /* CTOR */
        risk_result(const std::vector<iterator_pair> *const s, 
                    const T *const v, const T spot, const bool d) 
            : _s(s), _v(v), _spot(spot), _d(d) {  };
        
        /* Default CTOR for use in STL containers */
        risk_result() {  };
        
        /* DTOR */
        ~risk_result() {  };
        
        /* Allow default copy CTOR */
        
        /* Access function */
        const T *const          values()                            const { return _v;                                                  }
        const iterator_pair &   spot_grid(const unsigned int dim)   const { return (*_s)[dim];                                          }
        const T                 spot_price()                        const { return _spot;                                               }
        const unsigned int      size(const unsigned int dim)        const { return std::distance((*_s)[dim].first, (*_s)[dim].second);  }
        const unsigned int      dimensions()                        const { return _s->size();                                          }
        const bool              outer_grid_disconnected()           const { return _d;                                                  }
        
        /* Deep copy of s grid for construction of new results */
        std::vector<iterator_pair> *const deep_copy_s_grid() const
        {
            return new std::vector<iterator_pair>(_s->begin(), _s->end());
        }
    
    private :
        boost::shared_ptr<const std::vector<iterator_pair>>     _s;     /* Iterator marking the spot grid these results are defined on */
        const T                                             *   _v;     /* Option value */
        T                                                       _spot;  /* The spot price the results are defined at for single value results */
        bool                                                    _d;     /* Weather the outer s grid is disconnected */
};


template<class T>
class scenario_point
{
    public : 
        typedef typename risk_result<T>::iterator_pair iterator_pair;
        
        /* CTOR */
        /* This builds the null bump, all other scenarios should be built from the null bump */
        scenario_point(const std::vector<cashflow<T> *> &c, grid<T> &g, const T max_dt, 
                       const T s, const T r, const T v)
            : _cashflows(c), _timesteps(build_timesteps(max_dt)), _grid(g), _null_time(nullptr), 
              _time_off(0.0), _spot(s), _interest_rate(r), _volatility(v), 
              _cur_t(_timesteps->size() - 1), _disc_cashflow(false), _zombie(false) {  }

        /* Copy CTOR */
        scenario_point(scenario_point &s)
            : _cashflows(s._cashflows),  _timesteps(s._timesteps),  _grid(s._grid), 
              _null_time(s._null_time), _time_off(s._time_off), _spot(s._spot),
              _interest_rate(s._interest_rate), _volatility(s._volatility), _cur_t(s._cur_t),
              _disc_cashflow(s._disc_cashflow), _zombie(s._zombie) {  };

    
       /* Comparison for use in STL maps */
        bool operator<(const scenario_point &rhs) const
        {
            if (_spot < rhs._spot)
            {
                return true;
            }
            
            if (_interest_rate < rhs._interest_rate)
            {
                return true;
            }
    
            if (_volatility < rhs._volatility)
            {
                return true;
            }
                    
            return (_time_off < rhs._time_off);
        }
        
        
        void initialise(pde_solver<T> &calc_engine)
        {
            /* Process cashflows */
            _disc_cashflow = process_cashflows();
            
            /* Possible time offset step */
            if (_time_off > 0.0)
            {
                _grid.flip();
                calc_engine.solve(*this);
            }
            
            /* Move to the next time step */
            _grid.flip();
            --_cur_t;
        }
        
        
        /* Calculate a one step move */
        void iteration(pde_solver<T> &calc_engine)
        {

            /* If this scenario point is at an offset time it must be calculated from */
            /* its equivalent non offset time in order to include cashflows without */
            /* interpolations */
            if (_null_time != nullptr)
            {
                const grid<T> &g    = _null_time->_grid;
                const int size      = g.total_geometric_size();
                _disc_cashflow      = _null_time->_disc_cashflow;
                std::copy(&g.read_v_grid()[0], &g.read_v_grid()[size], &_grid.read_v_grid()[0]);
            }
            
            /* Calculate the step */
            assert(_cur_t >= 0);
            if (!_zombie)
            {
                calc_engine.solve(*this);
            }
            
            /* Process cashflows */
            if (_null_time == nullptr)
            {
                _disc_cashflow = process_cashflows();
            }
            
            /* Move to the next time step */
            _grid.flip();
            --_cur_t;
        }
        
        
        /* Setters */
        scenario_point& set_time_offset(const T t, const scenario_point<T> *n)
        { 
            _time_off  = t; 
            _null_time = n;
            return *this;
        }
        
        scenario_point& set_spot(const T s, const unsigned int dim)
        {
            _grid.shift(_cashflows, s, s - _spot, dim);
            _spot = s; 
            return *this;
        }
        
        scenario_point& set_interest_rate(const T r)    { _interest_rate    = r; return *this; }
        scenario_point& set_volatility(const T v)       { _volatility       = v; return *this; }
        
        /* Getters */
        bool done()                 const { return _cur_t == 0;             }
        T get_time()                const { return (*_timesteps)[_cur_t];   }
        T get_time_offset()         const { return _time_off;               }
        T get_spot()                const { return _spot;                   }
        T get_interest_rate()       const { return _interest_rate;          }
        T get_volatility()          const { return _volatility;             }
        bool is_disc_cashflow()     const { return _disc_cashflow;          }
        grid<T>& get_grid()               { return _grid;                   }
        const grid<T>& get_grid()   const { return _grid;                   }
        
        const risk_result<T> get_risk_result() const 
        {
            std::vector<iterator_pair> *s_grid_iters = new std::vector<iterator_pair>();
            for (int i = 0; i < _grid.dimensions(); i++)
            {
                s_grid_iters->push_back(std::make_pair(_grid.s_grid(i).begin(), _grid.s_grid(i).end()));
            }
            
            return risk_result<T>(s_grid_iters, _grid.read_v_grid(), _spot, _grid.outer_grid_disconnected());
        }
            
        bool begun() const 
        {
            return static_cast<unsigned int>(_cur_t) != (_timesteps->size() - 1);
        }

        T get_time_increment() const 
        {
            return (_time_off > 0.0) ? _time_off : ((*_timesteps)[_cur_t] - (*_timesteps)[_cur_t - 1]) + _time_off;
        }
        
        scenario_point& zombify() 
        { 
            _zombie = true; 
            return *this;
        }
        
        scenario_point& dezombify(const scenario_point &s) 
        {
            _zombie &= s._zombie;
            return *this;
        } 

    private :
        std::vector<cashflow<T> *>          _cashflows;     /* Any cashflows including the payoff           */
        boost::shared_ptr<std::vector<T> >  _timesteps;     /* The times that require simulating            */
        grid<T>                             _grid;          /* The PDE grid this scenario is on             */
        const scenario_point<T>            *_null_time;     /* This bump, but with no time offset           */
        T                                   _time_off;      /* The offset from the iteration time           */
        T                                   _spot;          /* Spot price                                   */
        T                                   _interest_rate; /* Interest rate                                */
        T                                   _volatility;    /* Volatility                                   */
        int                                 _cur_t;         /* Current time index                           */
        bool                                _disc_cashflow; /* Weather the last cashflow was discontinious  */
        bool                                _zombie;        /* Weather this scenario should be processed    */
        
        /* Build time steps based on cashflow times and maximum time step */
        std::vector<T> * build_timesteps(const T max_dt)
        {
            /* Build time steps */
            std::vector<T> required_dates;
            required_dates.push_back(0.0);
            std::back_insert_iterator<std::vector<T> > dates_inserter(required_dates);
            for (unsigned int i = 0; i < _cashflows.size(); i++)
            {
                _cashflows[i]->required_date(dates_inserter);
            }
    
            /* Clean up any duplicates */
            std::sort(required_dates.begin(), required_dates.end());
            typename std::vector<T>::iterator dates_iter = std::unique(required_dates.begin(), required_dates.end());
            required_dates.resize(dates_iter - required_dates.begin());
    
            /* Pad the dates for the minimum step size */
            std::vector<T> *padded_dates = new std::vector<T>();
            pad_data<T>(padded_dates, required_dates, max_dt);
                        
            return padded_dates;
        }
        
        /* Calculate any cachflows and add the to the grid */
        bool process_cashflows()
        {
            bool disc = false;
            const T t = (*_timesteps)[_cur_t];
            for (unsigned int i = 0; i < _cashflows.size(); i++)
            {
                if (_cashflows[i]->required_at(t))
                {
                    disc |= _cashflows[i]->add(_grid, t);
                }
            }
            
            return disc;
        }
};

#endif
