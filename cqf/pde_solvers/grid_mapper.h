#ifndef __GRID_MAPPER_H__
#define __GRID_MAPPER_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>
#include <vector>

#include "boost/bind.hpp"
#include "boost/noncopyable.hpp"


/* Mapper abstract base class */
template<class T>
class grid_mapper
{
    public :
        /* DTOR */
        virtual ~grid_mapper() {  };
    
        /* Map the grid points to a new grid */
        virtual std::vector<T>* map() const = 0;
};

/* Mapper to assign a uniformed grid */
template<class T>
class uniformed_grid_mapper : public grid_mapper<T>
{
    public : 
        /* CTOR */
        uniformed_grid_mapper(const T max, const T min, const int nas)
         : _inc((max - min) / static_cast<T>(nas)), _min(min), _nas(nas) { }
        
        std::vector<T>* map() const 
        {
            std::vector<T>* values = new std::vector<T>(_nas);
            for (int i = 0; i < _nas; i++)
            {
                (*values)[i] = _min + (i * _inc);
            }

            return values;
        }

    private :
        const T     _inc;
        const T     _min;
        const int   _nas;
};

/* Mapper to assign a non-uniformed grid */
template<class T>
class sinh_grid_mapper : public grid_mapper<T>
{
    public :
        /* CTOR */
        /* k is a point of interest */
        /* c is the proportion of grid points that lie in the region of k */
        sinh_grid_mapper(const T k, const T c, const T max, const T min, const int nas)
         : _k(k), _c(c), _max(max), _min(min), _nas(nas) { }
               
        std::vector<T>* map() const 
        {
            std::vector<T>* values = new std::vector<T>(_nas);
            const T nas_t = static_cast<T>(_nas);
            for (int i = 0; i < _nas; i++)
            {
                const T uni = asinh(-_k / _c) + ((i / nas_t) * (asinh((_max - _k) / _c) - asinh(-_k / _c)));
                (*values)[i] = _k + (_c * std::sinh(uni));
            }

            return values;
        }
        
    private :
        const T     _k;
        const T     _c;
        const T     _max;
        const T     _min;
        const int   _nas;
};


/* Representation of the FD grid. May be uniform or non-uniform */
template<class T>
class grid : private boost::noncopyable
{
    public :
        /* CTOR */
        grid(const grid_mapper<T> *const mapper) :
            _values(mapper->map()), _deltas(new std::vector<T>(_values->size()))
            {
                std::adjacent_difference(_values->begin(), _values->end(), _deltas->begin());
            };
        
        /* DTOR */
        ~grid()
        {
            delete _values;
            delete _deltas;
        };
        
        
        /* Accesss functions */
        int                     size()   const  { return _values->size();  }
        const std::vector<T>&   points() const  { return (*_values);       }
        
        
        T interpolate_result(const T *const v, const T a) const
        {
            typename std::vector<T>::const_iterator iter = std::find_if(_values->begin(), _values->end(), boost::bind(std::less_equal<T>(), a, _1));
            std::cout << "Found: " << (*iter) << std::endl;

            const int v_idx = std::distance(_values->begin(), iter);
            if ((*iter) == a)
            {
                return v[v_idx];
            }
            else
            {
                const T v3 = (*(iter + 1));
                const T v2 = (*iter--);
                const T v1 = (*iter--);
                const T v0 = (*iter);
                
                const T d = a - v1;
                const T d_sq = d * d;
                const T c0 = v3 - v2 - v0 + v1;
                const T c1 = v0 - v1 - c0;
                const T c2 = v2 - v0;
                const T c3 = v1;
                return (c0 * d * d_sq) + (c1 * d_sq) + (c2 * d) + c3;
                
//                const T d = u_s0 - l_s0;
//                const T w0 = (u_s0 - a) / d;
//                const T w1 = (a - l_s0) / d;
//                return (v[v_idx - 1] * w0) + (v[v_idx] * w1);
            }
        }
        
        int get_delta_coeffs(T *const coeffs, const unsigned int pos) const
        {
            /* Lower Boundary */
            if (pos == 0)
            {
                const T d1_p_d2 = (*_deltas)[1] + (*_deltas)[2];
                
                coeffs[0] = (-2.0 * (*_deltas)[1] - (*_deltas)[2]) / ((*_deltas)[1] * d1_p_d2);
                coeffs[1] = d1_p_d2 / ((*_deltas)[1] * (*_deltas)[2]);
                coeffs[2] = -(*_deltas)[1] / ((*_deltas)[2] * d1_p_d2);
            }
            /* Upper Boundary */
            else if (pos == (_values->size() - 1))
            {
                const T d1_p_d2 = (*_deltas)[pos - 1] + (*_deltas)[pos];
                
                coeffs[0] = (*_deltas)[pos] / ((*_deltas)[pos - 1] * d1_p_d2);
                coeffs[1] = (-(*_deltas)[pos - 1] - (*_deltas)[pos]) / ((*_deltas)[pos - 1] * (*_deltas)[pos]);
                coeffs[2] = ((*_deltas)[pos - 1] + 2.0 * (*_deltas)[pos]) / ((*_deltas)[pos] * d1_p_d2);
            }
            /* Middle */
            else
            {
                const T d1_p_d2 = (*_deltas)[pos] + (*_deltas)[pos + 1];
                
                coeffs[0] = -(*_deltas)[pos + 1] / ((*_deltas)[pos] * d1_p_d2);
                coeffs[1] = ((*_deltas)[pos + 1] - (*_deltas)[pos]) / ((*_deltas)[pos] * (*_deltas)[pos + 1]);
                coeffs[2] = (*_deltas)[pos] / ((*_deltas)[pos + 1] * d1_p_d2);
            }
            
            return 3;
        }
        
        int get_gamma_coeffs(T *const coeffs, const unsigned int pos) const
        {
            if ((pos != 0) && (pos != (_values->size() - 1)))
            {
                const T d1_p_d2 = (*_deltas)[pos] + (*_deltas)[pos + 1];
                
                coeffs[0] =  2.0 / ((*_deltas)[pos    ] * d1_p_d2);
                coeffs[1] = -2.0 / ((*_deltas)[pos    ] * (*_deltas)[pos + 1]);
                coeffs[2] =  2.0 / ((*_deltas)[pos + 1] * d1_p_d2);
                
                return 3;
            }
            else
            {
                coeffs[0] = 0.0;
                coeffs[1] = 0.0;
                coeffs[2] = 0.0;
                return 0;
            }
            
        }
        
        /* Debug function */
        void dump()
        {
            for (unsigned int i = 0; i < _values->size(); i++)
            {
                std::cout << (*_values)[i] << ", " << (*_deltas)[i] << std::endl;
            }
        }
        
    private :
        const std::vector<T>    *const _values;
        std::vector<T>          *const _deltas;
};

#endif
