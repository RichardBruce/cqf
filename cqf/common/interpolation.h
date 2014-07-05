#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <limits>


/* Base class for interpolation */
template<class T>
class interpolator
{
    public :
        /* CTOR */
        /* Doesnt take ownership or copy x and y */
        interpolator(const T *x, const T *y, const int s, const bool c, const int i = 0)
            : _x(x), _y(y), _s(s), _i(i), _c(c) { };
        
        /* Accept default copy CTOR and assignment operator */
        
        /* Virtual DTOR, class is intended for overriding */
        virtual ~interpolator() { };
        
        /* Interpolation pure virtual function */
        virtual T interpolate(const T a) = 0;
        
        /* Rebuild (re-construct) virtual function */
        virtual interpolator& rebuild(const T *x, const T *y, const int s, const bool c, const int i = 0)
        {
            _x = x;
            _y = y;
            _s = s;
            _c = c;
            _i = i;
            
            return *this;
        }
        
    protected :
        /* Return the first index not less than a */
        int search(const T a)
        {
            /* Sanity check */
            assert(_x != nullptr);
            assert(_y != nullptr);
            
            /* Values are not correlate, begin binary search */
            if (!_c)
            {
                _i = std::distance(&_x[0], std::lower_bound(&_x[0], &_x[_s], a));
            }
            /* Values are correlates and up, begin linear search upwards */
            else if (_x[_i] < a)
            {
                while ((_i < _s) && (_x[_i] < a))
                {
                    ++_i;
                }
            }
            /* Values are correlates and down, begin linear search downwards */
            else
            {
                while ((_i > 0) && (_x[_i - 1] >= a))
                {
                    --_i;
                }
            }
            
            return _i;
        }
        
        /* Access functions */
        T get_x(const int i) const { return _x[i]; }
        T get_y(const int i) const { return _y[i]; }
        
        int get_size() const { return _s; }
                
    private :
        const T     *   _x; /* X data       */
        const T     *   _y; /* Y data       */
        int             _s; /* Size         */
        int             _i; /* Index        */
        bool            _c; /* Correlated   */
};


/* Linear interpolation class */
template<class T>
class linear_interpolator : public interpolator<T>
{
    public :
        /* Empty CTOR for use in STL containers */
        linear_interpolator() : interpolator<T>(nullptr, nullptr, 0, false, 0) {  };
        
        /* CTOR */
        linear_interpolator(const T *x, const T *y, const int s, const bool c, const int i = 0)
            : interpolator<T>(x, y, s, c, i) {  };
            
        /* Accept default copy CTOR, DTOR and assignment operator */
        
        T interpolate(const T a)
        {
            /* Find index of first x not less than a */
            const int idx = this->search(a);
                        
            /* Exact match check */
            if (fabs(get_x(idx) - a) < (100.0 * std::numeric_limits<T>::epsilon()))
            {
                return get_y(idx);
            }

            /* Sanity check */
            assert((idx > 0) || !"Error: This means extrapolation not interpolation.");
            assert((idx < get_size()) || !"Error: This means extrapolation not interpolation.");
            
            /* Linear interpolation */
            return get_y(idx - 1) + 
                ((get_y(idx) - get_y(idx - 1)) * ((a - get_x(idx - 1)) / (get_x(idx) - get_x(idx - 1))));
        }
        
        /* Accept base rebuild */
    
    private :
        using interpolator<T>::get_x;
        using interpolator<T>::get_y;
        using interpolator<T>::get_size;
};


/* Linear interpolation class */
template<class T>
class cosine_interpolator : public interpolator<T>
{
    public :
        /* Empty CTOR for use in STL containers */
        cosine_interpolator() : interpolator<T>(nullptr, nullptr, 0, false, 0) {  };
        
        /* CTOR */
        cosine_interpolator(const T *x, const T *y, const int s, const bool c, const int i = 0)
            : interpolator<T>(x, y, s, c, i) {  };
            
        /* Accept default copy CTOR, DTOR and assignment operator */
        
        T interpolate(const T a)
        {
            /* Find index of first x not less than a */
            const int idx = this->search(a);

            /* Exact match check */
            if (fabs(get_x(idx) - a) < (100.0 * std::numeric_limits<T>::epsilon()))
            {
                return get_y(idx);
            }

            /* Sanity check */
            assert((idx > 0) || !"Error: This means extrapolation not interpolation.");
            assert((idx < get_size()) || !"Error: This means extrapolation not interpolation.");
            
            /* Cosine interpolation */
            const T mu = (a - get_x(idx - 1)) / (get_x(idx) - get_x(idx - 1));
            const T mu_sq = (1.0 - std::cos(mu * PI)) * 0.5;
            return ((get_y(idx - 1) * (1.0 - mu_sq)) + (get_y(idx) * mu_sq));
        }
    
        /* Accept base rebuild */
    
    private :
        using interpolator<T>::get_x;
        using interpolator<T>::get_y;
        using interpolator<T>::get_size;
};


/* Cubic polynomial interpolation class */
template<class T>
class cubic_poly_interpolator : public interpolator<T>
{
    public :
        /* Empty CTOR for use in STL containers */
        cubic_poly_interpolator() : interpolator<T>(nullptr, nullptr, 0, false, 0) {  };
        
        /* CTOR */
        cubic_poly_interpolator(const T *x, const T *y, const int s, const bool c, const int i = 0)
            : interpolator<T>(x, y, s, c, i) {  };
            
        /* Accept default copy CTOR, DTOR and assignment operator */
        
        T interpolate(const T a)
        {
            /* Find index of first x not less than a */
            const int idx = this->search(a);
            
            /* Exact match check */
            if (fabs(get_x(idx) - a) < (100.0 * std::numeric_limits<T>::epsilon()))
            {
                return get_y(idx);
            }

            /* Sanity check */
            assert((idx > 0) || !"Error: This means extrapolation not interpolation.");
            assert((idx < get_size()) || !"Error: This means extrapolation not interpolation.");
            
            /* Cubic polynomial interpolation */
            /* Result near the edge, linearly interpolate the result */
            if ((idx < 2) || (idx == (get_size() - 1)))
            {
                return get_y(idx - 1) + 
                    (get_y(idx) * ((a - get_x(idx - 1)) / (get_x(idx) - get_x(idx - 1))));
            }
            /* Result in the middle, cubic spline interpolate the result */
            else
            {
                return cubic_interpolation(
                       get_y(idx - 2), get_y(idx - 1), get_y(idx), get_y(idx + 1), 
                       get_x(idx - 2), get_x(idx - 1), get_x(idx), get_x(idx + 1), 
                       (a - get_x(idx - 1)) / (get_x(idx) - get_x(idx - 1)));
            }
        }
    
        /* Accept base rebuild */
    
    private :
        using interpolator<T>::get_x;
        using interpolator<T>::get_y;
        using interpolator<T>::get_size;
        
        T cubic_interpolation(const T y0, const T y1, const T y2, const T y3, 
            const T x0, const T x1, const T x2, const T x3, const T mu)
        {
            const T mu_sq = mu * mu;
    
            const T d1 = (y2 - y0) / (x2 - x0);
            const T d2 = (y3 - y1) / (x3 - x1);
    
            const T yy = y2 - y1;
            const T a0 = -(2.0 * yy) + d1 + d2;
            const T a1 =  (3.0 * yy) - (2.0 * d1) - d2;
            const T a2 = d1;
            const T a3 = y1;

           return ((a0 * mu * mu_sq) + (a1 * mu_sq) + (a2 * mu) + a3);
        }
};


/* Cubic spline interpolation class */
template<class T>
class cubic_spline_interpolator : public interpolator<T>
{
    public :
        /* Empty CTOR for use in STL containers */
        cubic_spline_interpolator() : interpolator<T>(nullptr, nullptr, 0, false, 0), _y2(nullptr) {  };
        
        /* CTOR */
        /* Uses custom matrix solving to reduce s^2 memory to s */
        cubic_spline_interpolator(const T *x, const T *y, const int s, const bool c, const int i = 0)
            : interpolator<T>(x, y, s, c, i), _y2(new T[s << 1])
            {
                calc_second_derivatives(x, y, s);
            };
          
        /* DTOR */
        ~cubic_spline_interpolator()
        {
            if (_y2 != nullptr)
            {
                delete [] _y2;
            }
        }
        
        /* Copy CTOR */
        cubic_spline_interpolator(const cubic_spline_interpolator &c)
            : interpolator<T>(c)
        {
            delete [] _y2;
            _y2 = new T[c._s << 1];
            memcpy(_y2, c._y2, c._s * sizeof(T));
        }
        
        /* Move CTOR */
        cubic_spline_interpolator(const cubic_spline_interpolator &&c)
            : interpolator<T>(c)
        {
            delete [] _y2;
            _y2 = c._y2;
            c._y2 = nullptr;
        }
        
        /* Assignment operator */
        cubic_spline_interpolator& operator=(const cubic_spline_interpolator &c)
        {
            /* Self assignment is benign */
            interpolator<T>::operator=(c);
            
            /* Safe copy y2 */
            delete [] _y2;
            _y2 = new T[c._s << 1];
            memcpy(_y2, c._y2, c._s * sizeof(T));
            
            return *this;
        };
        
        /* Move operator */
        cubic_spline_interpolator& operator=(cubic_spline_interpolator &&c)
        {
            /* Prevent self assignment */
            if (this == &c)
            {
                return *this;
            }
            
            interpolator<T>::operator=(c);
            
            /* Pilege y2 */
            if (_y2 != nullptr)
            {
                delete [] _y2;
            }
            _y2 = c._y2;
            c._y2 = nullptr;
            
            return *this;
        };
        
        T interpolate(const T a)
        {
            assert(_y2 != nullptr);
            
            /* Constants */
            static const T one_sixth = 1.0 / 6.0;
               
            /* Find index of first x not less than a */
            const int idx = this->search(a);
            
            /* Exact match check */
            if (fabs(get_x(idx) - a) < (100.0 * std::numeric_limits<T>::epsilon()))
            {
                return get_y(idx);
            }

            /* Sanity check */
            assert((idx > 0) || !"Error: This means extrapolation not interpolation.");
            assert((idx < get_size()) || !"Error: This means extrapolation not interpolation.");
            
            /* Cubic spline interpolation */
            const T d = get_x(idx) - get_x(idx - 1);
            const T c0 = (get_x(idx) - a) / d;
            const T c1 = (a - get_x(idx - 1)) / d;
            const T c2 = (c0 * c0 * c0) - c0;
            const T c3 = (c1 * c1 * c1) - c1;
            return (c0 * get_y(idx - 1)) + (c1 * get_y(idx)) +
                (((c2 * _y2[idx - 1]) + (c3 * _y2[idx])) * d * d * one_sixth);
        }
        
        cubic_spline_interpolator& rebuild(const T *x, const T *y, const int s, const bool c, const int i = 0)
        {
            /* Resize y2 */
            if (get_size() < s)
            {
                if (_y2 != nullptr)
                {
                    delete [] _y2;
                }
                
                _y2 = new T[s << 1];
            }
            
            /* Call base rebuild */
            interpolator<T>::rebuild(x, y, s, c, i);
            
            /* Solve the second derviatives of y */
            calc_second_derivatives(x, y, s);
            
            /* Return a reference to this */
            return *this;
        }
        
    private :
        using interpolator<T>::get_x;
        using interpolator<T>::get_y;
        using interpolator<T>::get_size;
                
        void calc_second_derivatives(const T *x, const T *y, const int s)
        {
            /* Lower boundary condition */
            _y2[0] = 0.0;
            _y2[s] = 0.0;
          
            /* Tridiagonal fill and elimination */
            for (int i = 1; i < s - 1; i++)
            {
                const T x_ratio     = (x[i] - x[i - 1]) / (x[i + 1] - x[i]);
                const T rho         = (x_ratio * _y2[i - 1]) + 2.0;
                const T grad_diff   = ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) - ((y[i] - y[i - 1]) / (x[i] - x[i - 1]));
                _y2[i]              = (x_ratio - 1.0) / rho;
                _y2[s + i]          = (6.0 * grad_diff / (x[i + 1] - x[i - 1]) - x_ratio * _y2[s + i - 1]) / rho;
            }
          
            /* Upper boundary condition */
            _y2[s - 1] = 0.0;
          
            /* Back substitution */
            for (int i = s - 2; i >=0; i--)
            {
                _y2[i] = (_y2[i] * _y2[i + 1]) + _y2[s + i];
            }
        }


        T   *_y2;   /* Second derivative of y */
};

#endif
