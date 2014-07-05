#ifndef __AGGREGATOR_H__
#define __AGGREGATOR_H__

#include "scenario_point.h"


template<class T>
class aggregator
{
    public :
        typedef typename risk_result<T>::iterator_pair iterator_pair;
        
        /* CTOR */
        aggregator(const scenario_point<T> &defining_risk, const T *s_bump_sizes, 
                   const T *weightings, const unsigned int order, const unsigned int size)
            : _defining_risk(defining_risk), _s_bump_sizes(s_bump_sizes), 
              _weightings(weightings), _risk(new std::vector<T>()),
              _inter(new std::vector<linear_interpolator<T>>(order + 1)),
              _order_p1(order + 1), _size(size)
            {
                _risk->reserve(_size);
            };
              
        /* Copy CTOR */
        aggregator(const aggregator<T> &a)
            : _defining_risk(a._defining_risk), _s_bump_sizes(a._s_bump_sizes), 
              _weightings(a._weightings), _risk(new std::vector<T>()), _order_p1(a._order_p1), _size(a._size)
            {
                _risk->reserve(_size);
            };
              
        /* DTOR */
        ~aggregator()
        {
            delete [] _s_bump_sizes;
            delete [] _weightings;
            delete _risk;
            delete _inter;
        }
        
        const risk_result<T> aggregate(const risk_result<T> *const results, const bool single_value)
        {
            /* Check results dimensions */
            for (unsigned int i = 0; i < _order_p1; i++)
            {
                assert((results[i].dimensions() == 1) || (results[i].dimensions() == 2 && results[i].outer_grid_disconnected()));
            }
            
            /* Get the grid sizes */
            const unsigned int dim_m1       = _defining_risk.get_grid().dimensions() - 1;
            const unsigned int nr_of_grids  = (dim_m1 == 0) ? 1 : _defining_risk.get_grid().size(0);
            
            /* Get a copy of the grids */
            std::vector<iterator_pair> *const def_s_grid = _defining_risk.get_grid().get_grid_bounds();
            --((*def_s_grid)[dim_m1].second);
            
            /* Clear results */
            const unsigned int def_s_grid_size = std::distance((*def_s_grid)[dim_m1].first, (*def_s_grid)[dim_m1].second);
            if (_risk->size() < def_s_grid_size)
            {
                _risk->resize(def_s_grid_size);
            }
    
            /* Aggregation */
            if (!single_value)
            {
                /* Move up the defining risk s grid until all s grids are > 0.0 */
                /* Move down the defining risk s grid until all s grid < max */
                for (unsigned int i = 0; i < _order_p1; i++)
                {
                    const iterator_pair &res_s_grid = results[i].spot_grid(dim_m1);
                    while ((*res_s_grid.first) > (*((*def_s_grid)[dim_m1].first) + _s_bump_sizes[i]))
                    {
                        ++((*def_s_grid)[dim_m1].first);
                    }

                    while ((*(res_s_grid.second - 1)) < (*((*def_s_grid)[dim_m1].second) + _s_bump_sizes[i]))
                    {
                        --((*def_s_grid)[dim_m1].second);
                    }
                }
                
                /* Process disconnected grids one at a time */
                int result_idx = 0;
                for (unsigned int i = 0; i < nr_of_grids; i++)
                {
                    /* Build interpolators */
                    build_interpolator(results, i, dim_m1);
            
                    /* Foreach defining risk with valid risks either side */
                    for (typename std::vector<T>::const_iterator j = (*def_s_grid)[dim_m1].first; j != (*def_s_grid)[dim_m1].second; j++)
                    {
                        /* Find the FV index related to the current s */
                        const T target_spot = (*j) + _s_bump_sizes[0];
//                        _risk->push_back(interpolate(results[0].spot_grid(dim_m1).first, results[0].spot_grid(dim_m1).second, &results[0].values()[i * results[0].size(dim_m1)], target_spot) * _weightings[0]);
                        (*_risk)[result_idx] = (*_inter)[0].interpolate(target_spot) * _weightings[0];
                        for (unsigned int k = 1; k < _order_p1; k++)
                        {
                            const T target_spot = (*j) + _s_bump_sizes[k];
                            (*_risk)[result_idx] += (*_inter)[k].interpolate(target_spot) * _weightings[k];
//                            _risk->back() += interpolate(results[k].spot_grid(dim_m1).first, results[k].spot_grid(dim_m1).second, &results[k].values()[i * results[k].size(dim_m1)], target_spot) * _weightings[k];
                        }
                        ++result_idx;
                    }
                }
            }
            else
            {
                /* Move the defining grid to the spot price */
                (*def_s_grid)[dim_m1].first = std::find((*def_s_grid)[dim_m1].first, (*def_s_grid)[dim_m1].second, _defining_risk.get_spot());
                (*def_s_grid)[dim_m1].second = (*def_s_grid)[dim_m1].first;

                /* Process disconnected grids one at a time */
                for (unsigned int i = 0; i < nr_of_grids; i++)
                {
                    /* Build interpolators */
                    build_interpolator(results, i, dim_m1);
            
                    /* Dot product of FV and weightings */
//                    _risk->push_back(interpolate(results[0].spot_grid(dim_m1).first, results[0].spot_grid(dim_m1).second, &results[0].values()[i * results[0].size(dim_m1)], results[0].spot_price()) * _weightings[0]);
                    (*_risk)[i] = (*_inter)[0].interpolate(results[0].spot_price()) * _weightings[0];
                    for (unsigned int j = 1; j < _order_p1; j++)
                    {
//                        _risk->back() += interpolate(results[j].spot_grid(dim_m1).first, results[j].spot_grid(dim_m1).second, &results[j].values()[i * results[j].size(dim_m1)], results[j].spot_price()) * _weightings[j];
                        (*_risk)[i] += (*_inter)[j].interpolate(results[j].spot_price()) * _weightings[j];
                    }
                }
            }
            
            return risk_result<T>(def_s_grid, _risk->data(), _defining_risk.get_spot(), _defining_risk.get_grid().outer_grid_disconnected());
        }
        
    private :
        aggregator<T>& operator=(const aggregator<T>&) {  };
        
        /* Set up an interpolator per result */
        void build_interpolator(const risk_result<T> *const r, const int g, const int d)
        {
            for (unsigned int i = 0; i < _order_p1; i++)
            {
                const int size = std::distance(r[i].spot_grid(d).first, r[i].spot_grid(d).second);
                (*_inter)[i] = linear_interpolator<T>(&(*r[i].spot_grid(d).first), &r[i].values()[g * r[i].size(d)], size, true, 0);
            }
        }
        
        const scenario_point<T> &               _defining_risk;
        const T                 *               _s_bump_sizes;
        const T                 *               _weightings;
        std::vector<T>          *               _risk;
        std::vector<linear_interpolator<T>>   * _inter;     /* Using linear interpolation because all values should be very close to exact */
        const unsigned int                      _order_p1;
        const unsigned int                      _size;
};

#endif
