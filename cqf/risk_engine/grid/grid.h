#ifndef __GRID_H__
#define __GRID_H__


#include "grid_mapper.h"

/* Representation of the FD grid. May be uniform or non-uniform */
template<class T>
class grid
{
    public :
        typedef std::pair<typename std::vector<T>::const_iterator,typename std::vector<T>::const_iterator> iterator_pair;
        
        /* CTOR */
        grid(const std::vector<boost::shared_ptr<const grid_mapper<T>>> mappers, 
             const std::vector<cashflow<T> *> &cashflows, const T spot, const bool disconnect)
            : _mappers(mappers), _lb(std::make_pair(false, 0.0)), _ub(std::make_pair(false, 0.0)), _s_grid(build_s_grids(cashflows, spot)),
              _deltas(calculate_deltas()), _dim_s_size(calculate_total_geometric_size()), 
              _v_grid(new T [_dim_s_size * 2]), _write_idx(0), 
              _disconnect(disconnect)
            {
                memset(_v_grid, 0, _dim_s_size * 2 * sizeof(T));
            };
            
        /* Copy CTOR */
        grid(grid &g)
            : _mappers(g._mappers), _lb(g._lb), _ub(g._ub), _s_grid(g._s_grid), 
              _deltas(g._deltas), _dim_s_size(g._dim_s_size), 
              _v_grid(new T [_dim_s_size * 2]), _write_idx(g._write_idx), 
              _disconnect(g._disconnect)
            {
                memcpy(_v_grid, g._v_grid, _dim_s_size * 2 * sizeof(T));
            }
        
        /* DTOR */
        ~grid()
        {
            delete [] _v_grid;
        };
        
        
        grid& shift(const std::vector<cashflow<T> *> &cashflows, const T spot, const T shift, const unsigned int dim)
        {
            assert((dimensions() == 1) || ((dimensions() == 2) && outer_grid_disconnected()));
            
            /* Deep copy and shift the grid up */
            _s_grid[dim].reset(new std::vector<T>(_s_grid[dim]->begin(), _s_grid[dim]->end()));
            for (unsigned int i = 0; i < _s_grid[dim]->size(); i++)
            {
                (*_s_grid[dim])[i] += shift;
            }
          
            /* Insert required s that didnt move */
            std::vector<T> required;
            _mappers[dim]->build_required_cashflows(required, cashflows, spot, &_lb, &_ub, true, false);
            _s_grid[dim]->insert(_s_grid[dim]->end(), required.begin(), required.end());
          
            /* Sort and clean s */
            typename std::vector<T>::iterator s_iter = std::remove_if(_s_grid[dim]->begin(), _s_grid[dim]->end(), boost::bind(std::less<T>(), _1, 0.0));
            std::sort(_s_grid[dim]->begin(), s_iter);
            s_iter = std::unique(_s_grid[dim]->begin(), s_iter);
            _s_grid[dim]->resize(s_iter - _s_grid[dim]->begin());
            remove_small_steps(_s_grid[dim], required);
            
            /* Rebuild the deltas */
            _deltas[dim].reset(new std::vector<T>(_s_grid[dim]->size()));
            std::adjacent_difference(_s_grid[dim]->begin(), _s_grid[dim]->end(), _deltas[dim]->begin());
            
            /* Recalculate size */
            _dim_s_size = calculate_total_geometric_size();
    
            /* Rebuild and clear v */
            delete [] _v_grid;
            _v_grid = new T [_dim_s_size * 2];
            memset(_v_grid, 0, _dim_s_size * 2 * sizeof(T));
            
            return *this;
        }
        
        
        /* Accesss functions */
        bool                    outer_grid_disconnected()       const   { return _disconnect;                           }
        int                     dimensions()                    const   { return _s_grid.size();                        }
        int                     size(const unsigned int dim)    const   { return _s_grid[dim]->size();                  }
        const std::vector<T>&   s_grid(const unsigned int dim)  const   { return (*_s_grid[dim]);                       }
        const T*                read_v_grid()                   const   { return &_v_grid[_write_idx ^ _dim_s_size];    }
        T*                      read_v_grid()                           { return &_v_grid[_write_idx ^ _dim_s_size];    }

        T*          write_v_grid()          { return &_v_grid[_write_idx];  }
        const T*    write_v_grid()  const   { return &_v_grid[_write_idx];  }
        void        flip()                  { _write_idx ^= _dim_s_size;    }
        
        /* Number of points in s grid */
        int total_geometric_size() const
        {
            return _dim_s_size;
        }
        
        /* Get const all s grids */
        const std::vector<boost::shared_ptr<std::vector<T>>> *const get_s_grids() const
        {
            return &_s_grid;
        }
        
        /* Deep copy of all s grids */
        std::vector<boost::shared_ptr<std::vector<T>>> deep_copy_s_grids(const T* bounds = nullptr) const
        {
            std::vector<boost::shared_ptr<std::vector<T>>> ret;
            for (unsigned int i = 0; i < _s_grid.size(); i++)
            {
                /* Pick optional bound for filtering */
                T bound = -0.1;
                if (bounds != nullptr)
                {
                    bound = bounds[i];
                }

                const typename std::vector<T>::iterator s_start = std::find_if(_s_grid[i]->begin(), _s_grid[i]->end(), boost::bind(std::less_equal<T>(), bound, _1));
                boost::shared_ptr<std::vector<T>> grid(new std::vector<T>(s_start, _s_grid[i]->end()));
                ret.push_back(grid);
            }
            
            return ret;
        }
        
        std::vector<iterator_pair> * get_grid_bounds() const
        {
            std::vector<iterator_pair> *s_grid_iters = new std::vector<iterator_pair>();
            for (int i = 0; i < dimensions(); i++)
            {
                s_grid_iters->push_back(std::make_pair(s_grid(i).begin(), s_grid(i).end()));
            }
            
            return s_grid_iters;
        }
        
        /* Replace an s grid */
        void replace_s_grid(const std::vector<T> &r, const unsigned int dim)
        {
            _s_grid[dim].reset(new std::vector<T>(r.begin(), r.end()));
            return;
        }

        
        /* Interpolate a results from the grid */
//        T interpolate_result(const T a) const
//        {
//            return interpolate(s_grid(), read_v_grid(), size(), a);
//        }
        
        /* PDE Solver coeffs for this grid */
        int get_delta_coeffs_lower(T *const coeffs, const unsigned int dim) const
        {
            const T d1_p_d2 = (*_deltas[dim])[1] + (*_deltas[dim])[2];
          
            coeffs[0] = (-2.0 * (*_deltas[dim])[1] - (*_deltas[dim])[2]) / ((*_deltas[dim])[1] * d1_p_d2);
            coeffs[1] = d1_p_d2 / ((*_deltas[dim])[1] * (*_deltas[dim])[2]);
            coeffs[2] = -(*_deltas[dim])[1] / ((*_deltas[dim])[2] * d1_p_d2);
            
            return 3;
        }
        
        int get_delta_coeffs_middle(T *const coeffs, const unsigned int dim, const unsigned int pos) const
        {
            const T d1_p_d2 = (*_deltas[dim])[pos] + (*_deltas[dim])[pos + 1];
          
            coeffs[0] = -(*_deltas[dim])[pos + 1] / ((*_deltas[dim])[pos] * d1_p_d2);
            coeffs[1] = ((*_deltas[dim])[pos + 1] - (*_deltas[dim])[pos]) / ((*_deltas[dim])[pos] * (*_deltas[dim])[pos + 1]);
            coeffs[2] = (*_deltas[dim])[pos] / ((*_deltas[dim])[pos + 1] * d1_p_d2);
            
            return 3;
        }
        
        int get_delta_coeffs_upper(T *const coeffs, const unsigned int dim) const
        {
            const unsigned int pos = _s_grid[dim]->size() - 1;
            const T d1_p_d2 = (*_deltas[dim])[pos - 1] + (*_deltas[dim])[pos];
          
            coeffs[0] = (*_deltas[dim])[pos] / ((*_deltas[dim])[pos - 1] * d1_p_d2);
            coeffs[1] = (-(*_deltas[dim])[pos - 1] - (*_deltas[dim])[pos]) / ((*_deltas[dim])[pos - 1] * (*_deltas[dim])[pos]);
            coeffs[2] = ((*_deltas[dim])[pos - 1] + 2.0 * (*_deltas[dim])[pos]) / ((*_deltas[dim])[pos] * d1_p_d2);
            
            return 3;
        }
        
        int get_gamma_coeffs_lower(T *const coeffs) const
        {
            coeffs[0] = 0.0;
            coeffs[1] = 0.0;
            coeffs[2] = 0.0;
            
            return 0;
        }
        
        int get_gamma_coeffs_middle(T *const coeffs, const unsigned int dim, const unsigned int pos) const
        {
            const T d1_p_d2 = (*_deltas[dim])[pos] + (*_deltas[dim])[pos + 1];
          
            coeffs[0] =  2.0 / ((*_deltas[dim])[pos    ] * d1_p_d2);
            coeffs[1] = -2.0 / ((*_deltas[dim])[pos    ] * (*_deltas[dim])[pos + 1]);
            coeffs[2] =  2.0 / ((*_deltas[dim])[pos + 1] * d1_p_d2);
          
            return 3;
        }
        
        int get_gamma_coeffs_upper(T *const coeffs) const
        {
            coeffs[0] = 0.0;
            coeffs[1] = 0.0;
            coeffs[2] = 0.0;
            
            return 0;
        }
        
        const std::pair<bool, T>& get_value_coeffs_lower() const
        {
            return _lb;
        }
        
        const std::pair<bool, T>& get_value_coeffs_upper() const
        {
            return _ub;
        }
        
        
        /* Debug function */
        void dump()
        {
            for (unsigned int i = 0; i < _s_grid.size(); i++)
            {
                for (unsigned int j = 0; j < size(); j++)
                {
                    std::cout << (*_s_grid[i])[j] << ", " << (*_deltas[i])[j] << ", " << read_v_grid()[(_s_grid[i]->size() * i) + j] << std::endl;
                }
                std::cout << std::endl;
            }
        }

    private :
        const std::vector<boost::shared_ptr<const grid_mapper<T>>>  _mappers;
        std::pair<bool, T>                                          _lb;
        std::pair<bool, T>                                          _ub;
        std::vector<boost::shared_ptr<std::vector<T>>>              _s_grid;
        std::vector<boost::shared_ptr<std::vector<T>>>              _deltas;
        int                                                         _dim_s_size;
        T                                                       *   _v_grid;
        int                                                         _write_idx;
        bool                                                        _disconnect;
        
        /* Prevent assignment */
        grid& operator=(grid &) {  };
        
        /* Build the s grids */
        std::vector<boost::shared_ptr<std::vector<T>>> build_s_grids(const std::vector<cashflow<T> *> &cashflows, const T spot)
        {
            std::vector<boost::shared_ptr<std::vector<T>>> s_grids;
            for (unsigned int i = 0; i < _mappers.size(); i++)
            {
                boost::shared_ptr<std::vector<T>> grid(_mappers[i]->map(cashflows, spot, &_lb, &_ub));
                s_grids.push_back(grid);
            }
            
            return s_grids;
        }
        
        /* Remove small grid steps */
        void remove_small_steps(boost::shared_ptr<std::vector<T>> s_grid, const std::vector<T> &required) const
        {
            const T space_threshold = 0.01;
            
            /* Copy down points with big enough spacings */
            /* The highest and lowest points must always be included */
            unsigned int write_idx = 0;
            unsigned int required_idx = 0;
            unsigned int candidate_idx = 0;
            while (candidate_idx < s_grid->size())
            {
                while (candidate_idx < s_grid->size())
                {
                    /* The gaps is big enough overwrite the next point */
                    if (((*s_grid)[candidate_idx] - (*s_grid)[write_idx]) > space_threshold)
                    {
                        (*s_grid)[++write_idx] = (*s_grid)[candidate_idx++];
                        break;
                    }
                    /* Candidate is required overwrite the write_idx */
                    else if (fabs(((*s_grid)[candidate_idx] - required[required_idx])) < 0.0001)
                    {
                        ++required_idx;
                        (*s_grid)[write_idx] = (*s_grid)[candidate_idx++];
                        break;
                    }
                    /* No good, next candidate */
                    else
                    {
                        ++candidate_idx;
                    }
                }
            }
            
            /* Drop unused points */
            s_grid->resize(write_idx + 1);
            
            return;
        }
        
        /* Calculate step size of s grids */
        std::vector<boost::shared_ptr<std::vector<T>>> calculate_deltas()
        {
            std::vector<boost::shared_ptr<std::vector<T>>> deltas;
            for (unsigned int i = 0; i < _s_grid.size(); i++)
            {
                boost::shared_ptr<std::vector<T>> delta(new std::vector<T>(_s_grid[i]->size()));
                std::adjacent_difference(_s_grid[i]->begin(), _s_grid[i]->end(), delta->begin());
                deltas.push_back(delta);
            }
            
            return deltas;
        }
        
        /* Product of the size of the grid accross all dimensions */
        int calculate_total_geometric_size() const
        {
            int total = 1;
            for (unsigned int i = 0; i < _s_grid.size(); i++)
            {
                total *= _s_grid[i]->size();
            }
            
            return total;
        }
};

#endif

