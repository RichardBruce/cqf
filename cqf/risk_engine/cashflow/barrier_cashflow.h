#ifndef __BARRIER_CASHFLOW_H__
#define __BARRIER_CASHFLOW_H__


/* Barrier enums */
enum barrier_dir_t { NO_BARRIER_DIR = 0, UP_BARRIER = 1, DOWN_BARRIER = 2 };
enum barrier_type_t { NO_BARRIER_TYPE = 0, IN_BARRIER = 1, OUT_BARRIER = 2 };

enum boundary_condition_t { NO_BOUNDARY_CONDITION = 0, LOWER_BOUNDARY_CONDITION = 1, UPPER_BOUNDARY_CONDITION = 2 };

template<class T>
class barrier_cashflow
{
    public :
        barrier_cashflow(const T b, const barrier_dir_t dir, const barrier_type_t type, const bool amer) 
            :   /* No barrier always in                                             In if > b, else in if < b */
              _b(((dir == NO_BARRIER_DIR) || (type == NO_BARRIER_TYPE)) ? -1.0 : ((static_cast<int>(dir) == static_cast<int>(type)) ? b : -b)),
              _amer(amer), _neg_s(static_cast<int>(dir) != static_cast<int>(type))
            {  };
        
        T operator()(const T v, T s) const
        {
            /* In if s < b, so negate s and b and apply > */
            s = _neg_s ? -s : s;
            
            /* In if > b, else 0 */
            if (s > _b)
            {
                return v;
            }
            else
            {
                return 0.0;
            }
        }
        
        void required_stock_price(std::back_insert_iterator<std::vector<T> > iter) const
        {
            /* Include the barrier in the grid if present */
            if (_neg_s || (_b >= 0.0))
            {
                iter = _neg_s ? -_b : _b;
            }
            
            return;
        }
        
        boundary_condition_t boundary_conditions(std::vector<T> &grid, T *const value) const
        {
            /* Only need for american barriers */
            if (_amer)
            {
                (*value) = 0.0;
                /* Remove gird points < the barrier */
                typename std::vector<T>::iterator remove_from;
                if (_neg_s)
                {
                    remove_from = std::remove_if(grid.begin(), grid.end(), boost::bind(std::less<T>(), -_b, _1));
                    grid.erase(remove_from, grid.end());
                    return UPPER_BOUNDARY_CONDITION;
                }
                /* Remove grid points > the barrier */
                else
                {
                    remove_from = std::remove_if(grid.begin(), grid.end(), boost::bind(std::greater<T>(), _b, _1));
                    grid.erase(remove_from, grid.end());
                    return LOWER_BOUNDARY_CONDITION;
                }
            }
            else
            {
                return NO_BOUNDARY_CONDITION;
            }
        }
        
    private :
        const T     _b;
        const bool  _amer;
        const bool  _neg_s;
};

#endif
