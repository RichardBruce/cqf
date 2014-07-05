#ifndef __SCENARIO_POINT_H__
#define __SCENARIO_POINT_H__


template<class T>
class scenario_point
{
    public : 
        /* CTOR */
        scenario_point(const T t, const T s, const T r, const T v)
            : _t(t), _s(s), _r(r), _v(v) {  }
            
        /* Allow default copy CTOR, assignment operator and DTOR */

    /* Comparison for use in STL maps */
    bool operator<(const scenario_point &rhs) const
    {
        if (this->_t < rhs._t)
        {
            return true;
        }
        
        if (this->_s < rhs._s)
        {
            return true;
        }
        
        if (this->_r < rhs._r)
        {
            return true;
        }
        
        return (this->_v < rhs._v);
    }
    
    
    /* Setters */
    scenario_point& set_time_to_maturity(const T t) { this->_t = t; return *this; }
    scenario_point& set_spot(const T s)             { this->_s = s; return *this; }
    scenario_point& set_interest_rate(const T r)    { this->_r = r; return *this; }
    scenario_point& set_volatility(const T v)       { this->_v = v; return *this; }
    
    /* Getters */
    T get_time_to_maturity()    const { return this->_t; }
    T get_spot()                const { return this->_s; }
    T get_interest_rate()       const { return this->_r; }
    T get_volatility()          const { return this->_v; }
        
    private :
        T _t;   /* The until maturity   */
        T _s;   /* Spot price           */
        T _r;   /* Interest rate        */
        T _v;   /* Volatility           */
};

#endif
