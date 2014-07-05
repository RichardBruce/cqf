#ifndef __RISK_H__
#define __RISK_H__


#include <map>

#include "grid_mapper.h"
#include "scenario_point.h"


/* Enums to describe methods of bumping */
enum bump_dir_t { DOWN = 0, UP = 1, BI = 2 };
enum bump_style_t { ADD = 0, MULT = 1, OVERRIDE = 2 };

/* Typedef to reduce typing */
class pingpong_buffer;
typedef map<scenario_point *, pingpong_buffer *> scenario_map;


template<class T>
class pingpong_buffer
{
    public :
        /* CTOR */
        pingpong_buffer(const int length)
        : _buffer(new T[length]), _write_idx(0), _length(length) {  };
        
        /* DTOR */
        ~pingpong_buffer()
        {
            delete [] _buffer;
        }
        
        /* Buffer access member */
        const T* read_buffer()  const { return &_buffer[_write_idx ^ _length];  }
        T* write_buffer()       const { return &_buffer[_write_idx];            }
        int length()            const { return _length;                         }
        
        void flip() { _write_idx ^= _length; }
        
    private :
        T *const    _buffer;
        int         _write_idx;
        const int   _length;
};


/* Virtual base class for risk calculations */
template<class T>
class risk
{
    public :
        /* CTOR */
        risk(const grid &g, const T bump_size = 0.01, 
            const bump_dir_t bump_dir = BI, const bump_style_t bump_style = ADD)
            : _grid(g), _bump_size(bump_size), _bump_dir(bump_dir), _bump_style(bump_style)
            {  }
           
        /* DTOR */
        virtual ~risk() {  }
        
        /* Virtual members */
        virtual void bump(scenario_map &scenarios, const scenario_point &null_bump) const = 0;
        virtual void aggregate(scenario_map &scenarios) const = 0;
        
    protected :
        const grid          _grid;
        const T             _bump_size;
        const bump_dir_t    _bump_dir;
        const bump_style_t  _bump_style;
        
        T apply_style_and_bump(const T data, const T bump) const 
        {
            switch (this->_bump_style)
            {
                case ADD :
                    return data + bump;
                    break;
                    
                case MULT :
                    return data + (bump * data);
                    break;
                    
                case OVERRIDE :
                    return bump;
                    break;
                    
                default :
                    assert(!"Error: Unknown bump style");
            }
        }
}


/* First order volatility risk */
template<class T>
class vega_risk : public risk<T>
{
    public :
        /* CTOR */
        vega_risk(const grid &g, const T bump_size = 0.01, const bump_dir_t bump_dir = BI, 
            const bump_style_t bump_style = ADD)
            : risk(g, bump_size, bump_dir, bump_style) { }
            
        void bump(scenario_map &scenarios, const scenario_point &null_bump) const
        {
            scenario_point *bumped;
            T bumped_vol;
            switch (this->_bump_dir)
            {
                case DOWN :
                    bumped_vol = apply_style_and_bump(null_bump.get_volatility(), -this->_bump_size);
                    bumped = new scenario_point(null_bump)->set_volatility(bumped_vol);
                    scenarios.insert(bumped, new pingpong_buffer(this->_grid.size()));
                    this->_dn = bumped;
                    this->_up = null_bump;
                    break;
                    
                case UP :
                    bumped_vol = apply_style_and_bump(null_bump.get_volatility(),  this->_bump_size);
                    bumped = new scenario_point(null_bump)->set_volatility(bumped_vol);
                    scenarios.insert(bumped, new pingpong_buffer(this->_grid.size()));
                    this->_dn = null_bump;
                    this->_up = bumped;
                    break;
                    
                case BI :
                    bumped_vol = apply_style_and_bump(null_bump.get_volatility(), -this->_bump_size);
                    bumped = new scenario_point(null_bump)->set_volatility(bumped_vol);
                    scenarios.insert(bumped, new pingpong_buffer(this->_grid.size()));
                    this->_dn = bumped;
                    
                    bumped_vol = apply_style_and_bump(null_bump.get_volatility(),  this->_bump_size);
                    bumped = new scenario_point(null_bump)->set_volatility(bumped_vol);
                    scenarios.insert(bumped, new pingpong_buffer(this->_grid.size()));
                    this->_up = bumped;
                    break;
                
                default :
                    assert(!"Error: Unknow bump direction");
            }
            
            return;
        }
        
        void aggregate(scenario_map &scenarios) const 
        {
            T *up;
            T *dn;
            T bump_size;
            switch (this->_bump_dir)
            {
                case DOWN, UP :
                    dn = scenarios[this->_dn]->read_buffer();
                    up = scenarios[this->_up]->read_buffer();
                    bump_size = this->_bump_size;
                    break;
                    
                case BI :
                    dn = scenarios[this->_dn]->read_buffer();
                    up = scenarios[this->_up]->read_buffer();
                    bump_size = 2.0 * this->_bump_size;
                    break;
                
                default :
                    assert(!"Error: Unknow bump direction");
            }
            
            apply_aggregation(up, dn, bump_size);
            return;
        }
    
    private :
        scenario_point _dn; /* Up scenario point used   */
        scenario_point _up; /* Down scenario point used */
        
        void apply_aggregation(const T *const up, const T *dn, const T bump_size) const
        {
            const T bump_size_inv = 1.0 / bump_size;
            for (int i = 0; i < this->_grid.size(); i++)
            {
                const T risk = (up[i] - dn[i]) / bump_size_inv;
            }
        }
}


#endif
