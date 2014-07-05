#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <assert.h>
#include <limits>
#include <initializer_list>

#include "constants.h"


/* TridiagonalPlus mean Tridiagonal with 3 elements on the top and bottom rows */
enum matrix_optimisation_flags_t { None = 0, Diagonal = 1, Tridiagonal = 2, TridiagonalPlus = 3, UpperTriangular = 4, LowerTriangular = 5, Symetric = 6 };


template<class T> class matrix
{
    public :
        /* CTOR */
        /* "Expert CTOR", for performance it may be better to release control
            of the "data" input and let the matrix delete it. This means that
            data must be a heap not a stack variable. Do this by setting 
            "delete_data" to true. It may also be better for performance if data
            is not deleted by the matrix class, in which case set "share_data" to true*/
        /* If you use this input "data" must be declares with new and in calling
           this CTOR you release control of "data" ie/ you must never read or 
           write data again. */
        matrix(T *const data, const int x, const int y, bool delete_data = false, bool share_data = false) :
         _x(x), _y(y), _share_data(share_data)
        {
            if (delete_data || share_data)
            {
                _data = data;
            }
            else
            {
                _data = new T [_x * _y];
                memcpy(_data, data, (_x * _y * sizeof(T)));
            }
        }
        
        matrix(std::initializer_list<std::initializer_list<T> > data)
        {
            typedef typename std::initializer_list<std::initializer_list<T> >::const_iterator outer_iter;
            typedef typename std::initializer_list<T>::const_iterator inner_iter;

            _data = new T [data.size() * data.begin()->size()];
            _x = data.size();
            _y = data.begin()->size();

            int idx = 0;
            for (outer_iter i = data.begin(); i != data.end(); i++)
            {
                assert(i->size() == static_cast<unsigned int>(_y));
                for (inner_iter j = i->begin(); j != i->end(); j++)
                {
                    _data[idx++] = (*j);
                }
            }
        }
        
        /* Copy CTOR */
        matrix(const matrix<T> &m) : _data(new T[m._x * m._y]), _x(m._x), _y(m._y)
        {
            memcpy(_data, m._data, (_x * _y * sizeof(T)));
        }
        
        /* Static factory to get an identity matrix */
        static matrix<T>* identity_matrix(const int n)
        {
            matrix<T> *identity = new matrix<T>(n, n);
            memset(identity->_data, 0, (n * n * sizeof(T)));
            for (int i = 0; i < n; i++)
            {
                identity->_data[(i * n) + i] = 1.0;
            }
            
            return identity;
        }
        
        /* Static factory to get a diagonal matrix from an array */
        static matrix<T> diagonal_matrix(const T *const data, const int x, const int y)
        {
            matrix<T> diagonal(x, y);
            memset(diagonal._data, 0, (x * y * sizeof(T)));
            for (int i = 0; i < std::min(x, y); i++)
            {
                diagonal._data[(i * y) + i] = data[i];
            }
            
            return diagonal;
        }
        
        /* Assignment operator */
        matrix& operator=(const matrix<T> &m)
        {
            _x = m._x;
            _y = m._y;
            
            delete [] _data;
            _data = new T[_x * _y];
            memcpy(_data, m._data, _x * _y * sizeof(T));
            
            return *this;
        }
        
        /* DTOR */
        ~matrix()
        {
            if (!_share_data)
            {
                delete [] _data;
            }
        }
        
        /* Standard operators */
        matrix operator+(const matrix<T> &m) const
        {
            assert(m._x == _x);
            assert(m._y == _y);
            
            matrix res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] + m._data[i];
            }
            
            return res;
        }
        
        matrix operator+(const T s) const
        {
            matrix res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] + s;
            }
            
            return res;
        }
        
        matrix& operator+=(const matrix<T> &m)
        {
            assert(m._x == _x);
            assert(m._y == _y);
            
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] += m._data[i];
            }
            
            return *this;
        }
        
        matrix& operator+=(const T s)
        {
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] += s;
            }
            
            return *this;
        }
        
        matrix operator-(const matrix<T> &m) const
        {
            assert(m._x == _x);
            assert(m._y == _y);
            
            matrix res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] - m._data[i];
            }
            
            return res;
        }
        
        matrix operator-(const T s) const
        {
            matrix res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] - s;
            }
            
            return res;
        }

        matrix& operator-=(const matrix<T> &m)
        {
            assert(m._x == _x);
            assert(m._y == _y);
            
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] -= m._data[i];
            }
            
            return *this;
        }
        
        matrix& operator-=(const T s)
        {
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] -= s;
            }
            
            return *this;
        }
                
        matrix operator*(const matrix<T> &m) const
        {
            /* This must be as high as m is wide */
            assert(m._x == _y);
            
            matrix res(_x, m._y);
            /* For each column of m */
            for (int i = 0; i < m._y; i++)
            {
                /* For each row of this */
                for (int j = 0; j < _x; j++)
                {
                    /* For each element in the row of this or column of m */
                    res._data[(j * m._y) + i] = _data[j * _y] * m._data[i];
                    for (int k = 1; k < _y; k++)
                    {
                        res._data[(j * m._y) + i] += _data[(j * _y) + k] * m._data[(k * m._y) + i];
                    }
                }
            }
            
            return res;
        }
        
        matrix operator*(const T s) const
        {
            matrix res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] * s;
            }
            
            return res;
        }
        
        matrix& operator*=(const matrix<T> &m)
        {
            /* This must be as high as m is wide */
            assert(m._x == _y);
            
            this->multiply_equals(m._data, m._y);
                                    
            return *this;
        }
        
        matrix& operator*=(const T s)
        {
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] *= s;
            }
            
            return *this;
        }
        
        matrix& operator/=(const T s)
        {
            for (int i = 0; i < (_x * _y); i++)
            {
                _data[i] /= s;
            }
            
            return *this;
        }
        
        matrix operator/(const T s) const
        {
            matrix<T> res(_x, _y);
            for (int i = 0; i < (_x * _y); i++)
            {
                res._data[i] = _data[i] / s;
            }
            
            return res;
        }
        
        /* Matrix operations without operators */
        T get_data(const int i, const int j) const
        {
            assert(i < _x);
            assert(j < _y);
            return _data[(i * _y) + j];
        }
        
        matrix<T> extract_column(const int col) const
        {
            assert(col < _x);
            
            matrix<T> res(1, _y);
            for (int i = 0; i < _y; i++)
            {
                res._data[i] = _data[(i * _x) + col];
            }
            
            return res;
        }
        
        T* extract_row(const int row) const
        {
            assert(row < _x);
            
            T* ret = new T [_y];
            memcpy(&ret[0], &_data[row * _y], (_y * sizeof(T)));
            
            return ret;
        }
        
        T trace() const
        {
            /* Trace is only defined to square matrices */
            assert(_x == _y);
            
            T sum = _data[0];
            for (int i = 1; i < _x; i++)
            {
                sum += _data[(i * _x) + i];
            }
            
            return sum;
        }

        
        matrix& transpose()
        {
            T *res = new T[_x * _y];
            for (int i = 0; i < _x; i++)
            {
                for (int j = 0; j < _y; j++)
                {
                    res[(j * _x) + i] = _data[(i * _y) + j];
                }
            }
            
            delete[] _data;
            _data = res;
            std::swap(_x, _y);
            
            return *this;
        }
        
        
        T determinant()
        {
            /* The determinant is not defined for non-square matrices */
            assert(_x == _y);
            return determinant_recursion(_data, _x);
        }
        
        
        matrix& invert()
        {
            /* The inverse is not defined for non-square matrices */
            assert(_x == _y);
            
            /* Determinant of the whole matrix, which must not be 0 */
            T det = determinant();
            assert(det != 0.0);
            T det_inv = 1.0 / det;

            T *tmp = new T[_x * _x];
            T *minor = new T[(_x - 1) * (_x - 1)];
            for(int j = 0; j < _x; j++)
            {
                for(int i = 0; i < _x; i++)
                {
                    /* Get the co-factor matrix of of the sub matrix at i, j */
                    get_minor(_data, minor, j, i, _x);
                    tmp[(i * _x) + j] = (((i + j) & 0x1) ? -1.0 : 1.0) * det_inv * determinant_recursion(minor , _x - 1);
                }
            }

            /* Clean up */
            delete [] _data;
            _data = tmp;
            delete [] minor;
            
            return *this;
        }
        
        matrix& gauss_solve(matrix<T> &eq, T *soln, matrix_optimisation_flags_t opts = None)
        {
            return gauss_solve(eq._data, soln, opts);
        }
        
        matrix& gauss_solve(T *eq, T *soln, matrix_optimisation_flags_t opts = None)
        {
            /* Gaussian elimination is not defined for non-square matrices */
            assert(_x == _y);

            /*  Optimisation enabled */
            const bool tri_opt = ((opts == Tridiagonal) || (opts == TridiagonalPlus));
            
            /* For each column */
            for (int i = 0; i < (_x - 1); i++)
            {
                /* Top row requires an extra row elimination for TridiagonalPlus */
                const int plus_row = ((i == 0) || (i == (_x - 3))) && (opts == TridiagonalPlus);
                
                /* Eliminate all rows unless tridiagonal matrix, then only eliminate 1 */
                const int rows_to_elim = tri_opt ? std::min(_x, (i + 2 + plus_row)) : _x;
                
                /* Eliminate all the rows below the diagonal */
                for (int j = (i + 1); j < rows_to_elim; j++)
                {
                    /* Subtract from this row/column heading right */
                    T factor = _data[(j * _x) + i] / _data[(i * _x) + i];
                    for (int k = i; k < (tri_opt ? std::min((i + 2), _x) : _x); k++)
                    {
                        _data[(j * _x) + k] -= _data[(i * _x) + k] * factor;
                    }
                    eq[j] -= eq[i] * factor;
                }
            }
            
            /* Back substitute */
            /* TODO - Tridiagonal and TridiagonalPlus only, make general */
            soln[_x - 1] = eq[_x - 1] / _data[(_x * _x) - 1];
            for (int i = _x - 2; i > 0; i--)
            {
                soln[i] = (eq[i] - (soln[i + 1] * _data[(i * _x) + i + 1])) / _data[(i * _x) + i];
            }
            
            /* Required for TridiagonalPlus, possibly better to max than if */
            soln[0] = (eq[0] - (soln[1] * _data[1]) - (soln[2] * _data[2])) / _data[0];
            
            return *this;        
        }
        
        matrix& householder_solve(T *eq, T *soln, matrix_optimisation_flags_t opts = None)
        {
            /* Household decomposition is not defined for none-square matrices */
            assert(_x == _y);
            
            /* For each column */
            T *v = new T[_x];
            for (int i = 0; i < _x - 1; i++)
            {
                /* Calculate the normal of the column and an appropriate projection */
                T normal = 0.0;
                for (int j = i; j < _x; j++)
                {
                    normal += _data[(j * _x) + i] * _data[(j * _x) + i];
                }
                normal = std::sqrt(normal);
                
                T d = normal;
                if (_data[(i * _x) + i] > 0.0)
                {
                    d = -normal;
                }
                
                /* Build the reflection vector and eliminate the first column */
                const T w = _data[(i * _x) + i] - d;
                const T f = std::sqrt(-2.0 * d * w);
                v[i] = w / f;
                _data[(i * _x) + i] = d;
                for (int j = i + 1; j < _x; j++)
                {
                    v[j] = _data[(j * _x) + i] / f;
                    _data[(j * _x) + i] = 0.0;
                }
                
                /* Apply the reflection to the remaining columns and the solution */
                for (int j = i + 1; j < _x; j++)
                {
                    T f = 0.0;
                    for (int k = i; k < _x; k++)
                    {
                        f += v[k] * _data[(k * _x) + j];
                    }
                    f *= 2.0;
                    
                    for (int k = i; k < _x; k++)
                    {
                        _data[(k * _x) + j] -= f * v[k];
                    }
                }
                
                T eq_f = 0.0;
                for (int j = i; j < _x; j++)
                {
                    eq_f += v[j] * eq[j];
                }
                eq_f *= 2.0;
                
                for (int j = i; j < _x; j++)
                {
                    eq[j] -= eq_f * v[j];
                }
           }
            
            /* Back substitute */
            soln[_x - 1] = eq[_x - 1] / _data[(_x * _x) - 1];
            for (int i = _x - 2; i >= 0; i--)
            {
                T sub = 0.0;
                for (int j = _x - 1; j > i; j--)
                {
                    sub += _data[(i * _x) + j] * soln[j];
                }
                soln[i] = (eq[i] - sub) / _data[(i * _x) + i];
            }
            
            /* Clean up */
            delete [] v;
            
            return *this;
        }
        
        
        matrix& sor_solve(T *eq, T *soln, T w, T tol, int max_iter, matrix_optimisation_flags_t opts = None)
        {
            /* SOR is not defined for none-square matrices */
            assert(_x == _y);
            
            /* This is optimised for trdiagonal matrices */
            assert(opts == Tridiagonal);

            /* Solve iteratively */
            T res = 0.0;
            int iter = 0;
            do 
            {
                /* Check for none convergence */
                assert(iter < max_iter);

                /* Unroll the loop to avoid range checks */
                T new_s = soln[0] + w * (eq[0] - (_data[1] * soln[1]) - (_data[0] * soln[0])) / _data[0];
                T diff  = new_s - soln[0];
                soln[0] = new_s;
                res     = (diff * diff);

                for (int i = 1; i < _x - 1; i++)
                {
                    new_s = soln[i] + w * (eq[i] - (_data[(i * _x) + i + 1] * soln[i + 1]) 
                        - (_data[(i * _x) + i] * soln[i]) - (_data[(i * _x) + i - 1] * soln[i - 1])) / _data[(i * _x) + i];
                        
                    diff = new_s - soln[i];
                    soln[i] = new_s;
                    res += (diff * diff);
                }
                
                new_s = soln[_x - 1] + w * (eq[_x - 1] - (_data[(_x * _x) - 1] * soln[_x - 1]) - (_data[(_x * _x) - 2] * soln[_x - 2])) / _data[(_x * _x) - 1];
                diff = new_s - soln[_x - 1];
                soln[_x - 1] = new_s;
                res += (diff * diff);

                ++iter;
            } while(res > tol);
            
            return *this;
        }
        
        
        matrix& lu_decomposition(matrix<T> **l)
        {
            /* Decomposition is not defined for non-square matrices */
            assert(_x == _y);
            
            /* For each column */
            T *lower = new T[_x * _x];
            for (int i = 0; i < _x; i++)
            {
                /* Zero the unaffected elements */
                for (int j = 0; j < i; j++)
                {
                    lower[(j * _x) + i] = 0.0;
                }
                lower[(i * _x) + i] = 1.0;
                
                /* Eliminate all the rows below the diagonal */
                for (int j = (i + 1); j < _x; j++)
                {
                    T factor = _data[(j * _x) + i] / _data[(i * _x) + i];
                    lower[(j * _x) + i] = factor;
                    
                    /* Subtract from this row/column heading left */
                    for (int k = i; k < _x; k++)
                    {
                        _data[(j * _x) + k] -= _data[(i * _x) + k] * factor;
                    }
                }
            }
            
            *l = new matrix(lower, _x, _x);
            return *this;
        }

        matrix& cholesky_decomposition()
        {
            /* Decomposition is not defined for non-square matrices */
            assert(_x == _y);
            
            /* For each row */
            T *tmp = new T[_x * _x];
            for (int i = 0; i < _x; i++)
            {
                /* For each column */
                /* Lower matrix */
                //T diagonal_inv = 1.0 / tmp[(i * _x) + i];
                for (int j = 0; j < i; j++)
                {
                    T sum = 0.0;
                    for (int k = 0; k < j; k++)
                    {
                        sum += (tmp[(i * _x) + k] * tmp[(j * _x) + k]);
                    }
                    std::cout << "Sum: " << sum << " data: " << _data[(i * _x) + j] << " divisor: " << tmp[(j * _x) + j] << std::endl;
                    tmp[(i * _x) + j] = (_data[(i * _x) + j] - sum) / tmp[(j * _x) + j];
                }
                
                /* Diagonal element */
                T sum = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sum += tmp[(i * _x) + j] * tmp[(i * _x) + j];
                }
                tmp[(i * _x) + i] = std::sqrt(_data[(i * _x) + i] - sum);
                std::cout << "Diagonal Sum: " << sum << " data: " << _data[(i * _x) + i] << std::endl;
                
                /* Fill in the zero elements */
                for (int j = (i + 1); j < _x; j++)
                {
                    tmp[(i * _x) + j] = 0.0;
                }
            }
            
            delete [] _data;
            _data = tmp;
            return *this;
        }
        
        
        /* ret should be _x by _x */
        matrix<T> single_value_decomposition(T *w)
        {
            /* This algorithm is only valid for over determined systems */
            assert(_x >= _y);
            
            bool flag;
            int i,its,j,jj,k,l = 0,nm = 0;
            T anorm,c,f,g,h,s,scale,x,y,z;

            int m = _x;
            int n = _y;
            T* rv1 = new T [n];
            matrix<T> ret(n, n);
            memset(ret._data, 0, (n * n * sizeof(T)));
            g = scale = anorm = 0.0;
            for (i = 0; i < n; i++)
            {
                l = i + 2;
                rv1[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m)
                {
                    for (k = i; k < m; k++)
                    {
                        scale += fabs(this->_data[(k * _y) + i]);
                    }
                    
                    if (fabs(scale) >= std::numeric_limits<T>::epsilon())
                    {
                        for (k = i; k < m; k++)
                        {
                            this->_data[(k * _y) + i] /= scale;
                            s += this->_data[(k * _y) + i] * this->_data[(k * _y) + i];
                        }
                        
                        f = this->_data[(i * _y) + i];
                        g = (f >= 0.0) ? -std::sqrt(s) : std::sqrt(s);
                        h = f * g - s;
                        this->_data[(i * _y) + i]=f-g;
                        for (j = l - 1; j < n; j++)
                        {
                            for (s = 0.0, k = i; k < m; k++)
                            {
                                s += this->_data[(k * _y) + i]*this->_data[(k * _y) + j];
                            }
                            
                            f = s / h;
                            for (k = i; k < m; k++)
                            {
                                this->_data[(k * _y) + j] += f*this->_data[(k * _y) + i];
                            }
                        }
                        
                        for (k = i; k < m; k++)
                        {
                            this->_data[(k * _y) + i] *= scale;
                        }
                    }
                }
               
                w[i] = scale * g;
                g = s = scale = 0.0;
                if (((i + 1) <= m) && ((i + 1) != n))
                {
                    for (k = l - 1; k < n; k++)
                    {
                        scale += fabs(this->_data[(i * _y) + k]);
                    }
                    
                    if (scale != 0.0)
                    {
                        for (k = l - 1; k < n; k++)
                        {
                            this->_data[(i * _y) + k] /= scale;
                            s += this->_data[(i * _y) + k]*this->_data[(i * _y) + k];
                        }
                        
                        f = this->_data[(i * _y) + l - 1];
                        g = (f >= 0.0) ? -std::sqrt(s) : std::sqrt(s);
                        h = f * g - s;
                        this->_data[(i * _y) + l - 1] = f - g;
                        for (k = l - 1; k < n; k++)
                        {
                            rv1[k]=this->_data[(i * _y) + k]/h;
                        }
                        
                        for (j = l - 1; j < m; j++)
                        {
                            for (s = 0.0, k = l - 1; k < n; k++)
                            {
                                s += this->_data[(j * _y) + k] * this->_data[(i * _y) + k];
                            }
                            
                            for (k = l - 1; k < n; k++)
                            {
                                this->_data[(j * _y) + k] += s * rv1[k];
                            }
                        }
                        
                        for (k = l - 1; k < n; k++)
                        {
                            this->_data[(i * _y) + k] *= scale;
                        }
                    }
                }
                anorm = std::max(anorm,(T)(fabs(w[i]) + fabs(rv1[i])));
            }
           
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n-1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                        {
                            ret._data[(j * _y) + i] = (this->_data[(i * _y) + j] / this->_data[(i * _y) + l]) / g;
                        }
                        
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0, k = l; k < n; k++)
                            {
                                s += this->_data[(i * _y) + k] * ret._data[(k * _y) + j];
                            }
                            
                            for (k = l; k < n; k++)
                            {
                                ret._data[(k * _y) + j] += s * ret._data[(k * _y) + i];
                            }
                        }
                    }
                    
                    for (j = l; j < n; j++)
                    {
                        ret._data[(i * _y) + j] = 0.0;
                        ret._data[(j * _y) + i] = 0.0;
                    }
                }
                
                ret._data[(i * _y) + i] = 1.0;
                g = rv1[i];
                l = i;
            }
      
            for (i = std::min(m, n) - 1; i >= 0; i--)
            {
                l = i + 1;
                g = w[i];
                for (j = l; j < n; j++)
                {
                    this->_data[(i * _y) + j] = 0.0;
                }
                
                if (g != 0.0)
                {
                    g = 1.0 / g;
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                        {
                            s += this->_data[(k * _y) + i] * this->_data[(k * _y) + j];
                        }
                        
                        f = (s / this->_data[(i * _y) + i]) * g;
                        for (k = i; k < m; k++)
                        {
                            this->_data[(k * _y) + j] += f * this->_data[(k * _y) + i];
                        }
                    }
                    for (j = i; j < m; j++)
                    {
                        this->_data[(j * _y) + i] *= g;
                    }
                }
                else
                {
                    for (j = i; j < m; j++)
                    {
                        this->_data[(j * _y) + i] = 0.0;
                    }
                }
                    
                this->_data[(i * _y) + i] = this->_data[(i * _y) + i] + 1;
            }
           
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 0; its < 30; its++)
                {
                    flag = true;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l - 1;
                        if (fabs(rv1[l]) + anorm == anorm)
                        {
                            flag = false;
                            break;
                        }
                        
                        if (fabs(w[nm]) + anorm == anorm)
                        {
                            break;
                        }
                    }
                    
                    if (flag)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i < k + 1; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if (fabs(f) + anorm == anorm)
                            {
                                break;
                            }
                            
                            g = w[i];
                            h = std::sqrt((f * f) + (g * g));
                            w[i] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < m; j++)
                            {
                                y = this->_data[(j * _y) + nm];
                                z = this->_data[(j * _y) + i];
                                this->_data[(j * _y) + nm] = y * c + z * s;
                                this->_data[(j * _y) + i] = z * c - y * s;
                            }
                        }
                    }
                    
                    z = w[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k] = -z;
                            for (j = 0; j < n; j++)
                            {
                                ret._data[(j * _y) + k] = -ret._data[(j * _y) + k];
                            }
                        }
                        break;
                    }
                    
                    /* Failed to converge */
                    if (its == 30)
                    {
                        return ret;
                    }
                    
                    x = w[l];
                    nm = k - 1;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = std::sqrt((f * f) + 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + (f >= 0.0 ? fabs(g) : -fabs(g)))) - h)) / x;
                    c = s = 1.0;
                    
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = std::sqrt((f * f) + (h * h));
                        rv1[j]= z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 0; jj < n; jj++)
                        {
                            x = ret._data[(jj * _y) + j];
                            z = ret._data[(jj * _y) + i];
                            ret._data[(jj * _y) + j] = x * c + z * s;
                            ret._data[(jj * _y) + i] = z * c - x * s;
                        }
                        
                        z = std::sqrt((f * f) + (h * h));
                        w[j] = z;
                        if (z)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 0; jj < m; jj++)
                        {
                            y = this->_data[(jj * _y) + j];
                            z = this->_data[(jj * _y) + i];
                            this->_data[(jj * _y) + j] = y * c + z * s;
                            this->_data[(jj * _y) + i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k] = x;
                }
            }
            
            /* Clean up */
            delete [] rv1;
            
            return ret;
        }
        
        
        matrix<T> moore_penrose_pseudo_inverse(T *w)
        {
            /* Perform single value decomposition */
            matrix<T> v = this->single_value_decomposition(w);
            
            /* Take the recipricol of w to build the pseudo inverse */
            for (int j = 0; j < _y; j++)
            {
                if (w[j] != 0.0)
                {
                    w[j] = 1.0 / w[j];
                }
            }
            matrix<T> svd_w_inv = matrix<T>::diagonal_matrix(w, _y, _y);
            return (*this) * (svd_w_inv * v.transpose());


            /* Take the recipricol of w */
            for (int i = 0; i < _y; i++)
            {
                T scale = 0.0;
                if (w[i])
                {
                    scale = 1.0 / w[i];
                }
                
                /* Multiply w and v */
                for (int j = 0; j < _y; j++)
                {
                    v._data[(j * _y) + i] *= scale;
                }
            }
            
            /* Transpose u and return the pseudo inverse */
            this->transpose();
            return (*this) * v;
        }
        
        
        /* Eigen value and eigen vector functions */        
        /* Uses the Jacobi method to find eigen values */
        T* eigen_values(const T max_error, const matrix_optimisation_flags_t opts)
        {
            /* This algorithm is only valid for symetric matrices */
            assert(opts == Symetric);
            
            T asumsq = this->sum_of_upper_triangular_squares();
            while (asumsq > max_error)
            {
                matrix<T> Anext = this->jacobi_rotate();
                memcpy(this->_data, Anext._data, (_x * _x * sizeof(T)));
                asumsq = this->sum_of_upper_triangular_squares();
            }
            
            T* eigen_values = new T[_x];
            for (int i = 0; i < _x; i++)
            {
                eigen_values[i] = _data[(i * _x) + i];
            }
            
            return eigen_values;
        }

        /*  Uses the jacobi method to find eigen vectors */
        matrix<T>* eigen_vectors(const T max_error, const matrix_optimisation_flags_t opts)
        {
            /* This algorithm is only valid for symetric matrices */
            assert(opts == Symetric);
            
            matrix<T>* ev = matrix<T>::identity_matrix(_x);
            T asumsq = this->sum_of_upper_triangular_squares();
            while(asumsq > max_error)
            {
                matrix<T> Anext = this->jacobi_rotate();
                this->jacobi_eigen_vector_rotate(ev);
                memcpy(this->_data, Anext._data, (_x * _x * sizeof(T)));
                asumsq = this->sum_of_upper_triangular_squares();
            }
            
//            for (int i = 0; i < _x; i++)
//            {
//                for (int j = _x - 1; j >= 0; j--)
//                {
//                    ev->_data[(i * _x) + j] /= ev->_data[(i * _x)];
//                }
//            }
    
            return ev;
        }
        
        void eigen_system(T d[], matrix_optimisation_flags_t opts, const bool with_vectors = true)
        {
            /* This algorithm is only valid for symetric matrices */
            assert(opts == Symetric);
            
            T *e = new T[_x];
            for (int i = _x - 1; i >= 1; i--)
            {
                int l = i - 1;
                T h = 0.0;
                T scale = 0.0;
                if (l > 0)
                {
                    for (int k = 0; k <= l; k++)
                    {
                        scale += fabs(_data[(i * _x) + k]);
                    }
                
                    if (scale == 0.0)
                    {
                        e[i] = _data[(i * _x) + l];
                    }
                    else
                    {
                        for (int k = 0; k <= l; k++)
                        {
                            _data[(i * _x) + k] /= scale;
                            h += _data[(i * _x) + k] * _data[(i * _x) + k];
                        }
            
                        T f = _data[(i * _x) + l];
                        T g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                        e[i] = scale * g;
                        h -= f * g;
                        _data[(i * _x) + l] = f - g;
                        f = 0.0;
                        for (int j = 0; j <= l; j++)
                        {
                            /* Only required for eigen vectors */
                            if (with_vectors)
                            {
                                _data[(j * _x) + i] = _data[(i * _x) + j] / h;
                            }

                            g = 0.0;
                            for (int k = 0; k <= j; k++)
                            {
                                g += _data[(j * _x) + k] * _data[(i * _x) + k];
                            }
                            
                            for (int k = j + 1; k <= l; k++)
                            {
                                g += _data[(k * _x) + j] * _data[(i * _x) + k];
                            }
                            e[j] = g / h;
            
                            f += e[j] * _data[(i * _x) + j];
                        }
            
                        T hh = f / (h + h);
                        for (int j = 0; j <= l; j++)
                        {
                            f = _data[(i * _x) + j];
                            e[j] = g = e[j] - hh * f;
                            for (int k = 0; k <= j; k++)
                            {
                                _data[(j * _x) + k] -= (f * e[k] + g * _data[(i * _x) + k]);
                            }
                        }
                    }
                }
                else
                {
                    e[i] = _data[(i * _x) + l];
                }
                d[i] = h;
            }
        
            /* Only required for eigen vectors */
            d[0] = 0.0;
            e[0] = 0.0;
            for (int i = 0; i < _x; i++)
            {
                int l = i - 1;
                if (with_vectors)
                {
                    if (d[i])
                    {
                        for (int j = 0; j <= l; j++)
                        {
                            T g = 0.0;
                            for (int k = 0; k <= l; k++)
                            {
                                g += _data[(i * _x) + k] * _data[(k * _x) + j];
                            }
                            for (int k = 0; k <= l; k++)
                            {
                                _data[(k * _x) + j] -= g * _data[(k * _x) + i];
                            }
                        }
                    }
                }
        
                d[i] = _data[(i * _x) + i];
                
                if (with_vectors)
                {
                    _data[(i * _x) + i] = 1.0;
                    for (int j = 0; j <= l; j++)
                    {
                        _data[(j * _x) + i] = _data[(i * _x) + j] = 0.0;
                    }
                }
            }

            for (int i = 1; i < _x; i++)
            {
                e[i - 1] = e[i];
            }
            
            e[_x - 1] = 0.0;
            for (int l = 0; l < _x; l++)
            {
                int iter = 0;
                int m;
                do 
                {
                    for (m = l; m < _x - 1; m++) 
                    {
                        T tmp = fabs(d[m]) + fabs(d[m + 1]);
                        if (static_cast<T>(fabs(e[m]) + tmp) == tmp)
                        {
                            break;
                        }
                    }
        
                    if (m != l)
                    {
                        if (iter++ == 30)
                        {
                            std::cout << "Eigen system failed to converge" << std::endl;
                            return;
                        }
                        
                        T g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                        T r = sqrt((g * g) + 1.0);
                        g = d[m] - d[l] + e[l] / (g + ((g < 0.0) ? -fabs(r) : fabs(r)));
                        T s = 1.0;
                        T c = 1.0;
                        T p = 0.0;
        
                        int i;
                        for (i = m - 1; i >= l; i--)
                        {
                            T f = s * e[i];
                            T b = c * e[i];
                            r = std::sqrt((f * f) + (g * g));
                            e[i + 1] = r;
                            if (r == 0.0)
                            {
                                d[i + 1] -= p;
                                e[m] = 0.0;
                                break;
                            }
        
                            s = f / r;
                            c = g / r;
                            g = d[i + 1] - p;
                            r = (d[i] - g) * s + 2.0 * c * b;
                            d[i + 1] = g + (p = s * r);
                            g = c * r - b;
        
                            /* Only required for eigen vectors */
                            if (with_vectors)
                            {
                                for (int k = 0; k < _x; k++)
                                {
                                    f = _data[(k * _x) + i + 1];
                                    _data[(k * _x) + i + 1] = s * _data[(k * _x) + i] + c * f;
                                    _data[(k * _x) + i] = c * _data[(k * _x) + i] - s * f;
                                }
                            }
                        }
        
                        if ((r == 0.0) && (i >= l))
                        {
                            continue;
                        }
                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0;
                    }
                } while (m != l);
            }
            
            /* Clean up */
            delete [] e;
            
            return;
        }

        /* Debug functions */
        matrix& dump(const std::string &s)
        {
            if (!s.empty())
            {
                std::ofstream file(s.c_str());
                assert(file.is_open());
                dump(file);
                file.close();
            }
            
            return *this;
        }
        
        
        matrix& dump(std::ostream &s)
        {
            for (int i = 0; i < _x; i++)
            {
                for (int j = 0; j < _y; j++)
                {
                    s << _data[(i * _y) + j] << ", ";
                }
                
                s << std::endl;
            }
            return *this;
        }
    
    private :
        /* Private CTOR for when we will fill the data ourselves */
        matrix(const int x, const int y)
        {
            _x = x;
            _y = y;
            _data = new T[_x * _y];
        }
        
        void multiply_equals(const T *const m_data, const int m_y)
        {
            /* For each column of m */
            T *res = new T[_x * m_y * sizeof(T)];
            for (int i = 0; i < m_y; i++)
            {
                /* For each row of this */
                for (int j = 0; j < _x; j++)
                {
                    /* For each element in the row of this or column of m */
                    res[(j * m_y) + i] = _data[j * _y] * m_data[i];
                    for (int k = 1; k < _y; k++)
                    {
                        res[(j * m_y) + i] += _data[(j * _y) + k] * m_data[(k * m_y) + i];
                    }
                }
            }
            
            
            delete [] _data;
            _data = res;
            _y = m_y;
            
            return;
        }
        
        void get_minor(const T *const src, T *const dst, const int r, const int c, const int o) const
        {
            int col_num = 0;
            int row_num = 0;
        
            /* For each row that isnt the current row */
            for(int i = 0; i < o; i++ )
            {
                if(i != r)
                {
                    /* For each column that isnt the current column */
                    col_num = 0;
                    for(int j = 0; j < o; j++ )
                    {
                        if(j != c)
                        {
                            dst[(row_num * (o - 1)) + col_num] = src[(i * o) + j];
                            col_num++;
                        }
                    }
                    row_num++;
                }
            }

            return;
        }
        
        T determinant_recursion(const T *const m, const int o) const
        {
            /* Base case for recursion */
            if (o == 1)
            {
                return m[0];
            }
            
            /* Get the minor and recurse */
            T det = 0.0;
            T *minor = new T[(o - 1) * (o - 1)];
            for(int i = 0; i < o; i++ )
            {
                if (m[i] != 0.0)
                {
                    get_minor(m, minor, 0, i, o);
                    det += ((i & 0x1) ? -1.0 : 1.0) * m[i] * determinant_recursion(minor, o - 1);
                }
            }

            /* Clean up */
            delete [] minor;

            return det;
        }
        
        T sum_of_upper_triangular_squares() const
        {
            /* Triangular matrices must be square */
            assert(_x == _y);
            
            /* For each row */
            T sum = 0.0;
            for (int i = 0; i < _x; i++)
            {
                /* For each column */
                for (int j = i + 1; j < _x; j++)
                {
                    sum += (_data[(i * _x) + j] * _data[(i * _x) + j]);
                }
            }
            
            return sum;
        }
        
        matrix<T>* jacobi_rotation_matrix() const
        {
            /* Find the largest value */
            T maxval = -1.0;
            int m_i;
            int m_j;
            for (int i = 0; i < _x; i++)
            {
                for (int j = i + 1; j < _x; j++)
                {
                    T tmp = fabs(_data[(i * _x) + j]);
                    if (tmp > maxval)
                    {
                        maxval = tmp;
                        m_i = i;
                        m_j = j;
                    }
                }
            }
        
            /* Calculate the rotation needed to remove this value */
            T angle;
            if (_data[(m_i * _x) + m_i] == _data[(m_j * _x) + m_j])
            {
                angle = 0.25 * PI * ((_data[(m_i * _x) + m_j] > 0.0) ? 1.0 : -1.0);
            }
            else
            {
                angle = 0.5 * std::atan(2.0 * _data[(m_i * _x) + m_j] / (_data[(m_i * _x) + m_i] - _data[(m_j * _x) + m_j]));
            }
            
            /* Build the rotation matrix */
            matrix<T> *rotation = matrix<T>::identity_matrix(_x);
            rotation->_data[(m_i * _x) + m_i] = std::cos(angle);
            rotation->_data[(m_j * _x) + m_i] = std::sin(angle);
            rotation->_data[(m_i * _x) + m_j] = -std::sin(angle);
            rotation->_data[(m_j * _x) + m_j] = std::cos(angle);
            return rotation;
        }
         
        /* Rotate the largest element of this matrix */
        matrix<T> jacobi_rotate() const
        {
            matrix<T> *rotation = this->jacobi_rotation_matrix();
            matrix<T> ret = rotation->transpose() * ((*this) * (*rotation));
            
            /* Clean up */
            delete rotation;
            
            return ret;
        }

        /* Track the effects of rotation this matrix on another eq/ Eigen vectors */
        void jacobi_eigen_vector_rotate(matrix<T> *const ev) const
        {
            matrix<T> *rotation = this->jacobi_rotation_matrix();
            (*ev) *= (*rotation);
            
            /* Clean up */
            delete rotation;
        }

        
        T*      _data;          /* Data in the matrix                       */ 
        int     _x;             /* Height of the matrix                     */
        int     _y;             /* Width of the matrix                      */
        bool    _share_data;    /* Share the matrix data with third party   */
};





//void tred2(float **a, int n, float d[], float e[]);
//void tqli(float d[], float e[], int n, float **z);
float magn(float a, float b);
//void nrerror(char error_text[]);
//void exit(int status);
//void free(void *block);
//
//int main( void )
//{
//        double **data, mat1[ 2 ], mat2[ 2 ];
//
//        if ( ( data = dmatrix( 2, 2 ) ) == NULL )
//            exit(1);
//
//        data[ 0 ][ 0 ] =  0.35;
//        data[ 0 ][ 1 ] = -0.89;
//        data[ 1 ][ 0 ] = -0.89;
//        data[ 1 ][ 1 ] = -0.2175;
//
//        tred2( data, 2, mat1, mat2 );
//    tqli( mat1 ,mat2 , 2, data );
//    printf("%f == %f", mat1[ 1 ], mat1[ 2 ] );
//
//        free( data );
//        return 0;
//
//}
//
///* Function returns a NULL pointer on failure! */
//
//double **dmatrix( long rows, long columns )
//{
//    double **m;
//    long i;
//
//    /* Check for pathological cases (there is no matrix with less then    */
//    /* one or negative elements)                                          */
//
//    if ( rows < 1  || columns < 1 )
//        return NULL;
//
//    /* Allocate memory for the matrix as well as for an array of pointers */
//    /* to its rows - return NULL pointer on failure                       */
//
//    m = malloc( ( size_t ) rows * ( ( size_t ) columns *
//                                  sizeof( double ) + sizeof( double * ) ) );
//
//    if ( m == NULL )
//        return NULL;                                 /* Allocation failed */
//
//    /* Initialize the array of pointers to the rows of the matrix         */
//
//    m[ 0 ] = ( double * ) ( m + rows );
//    for ( i = 1; i < rows; i++ )
//        m[ i ] = m[ i - 1 ] + columns;
//
//    /* Return pointer to the pointer array to the rows of the matrix      */
//
//    return m;
//}
//
//float sign(float a, float b)
//{
//    return (b < 0.0) ? -fabs(a) : fabs(b);
//}
//
//float magn(float x, float y)
//{ 
//    return sqrt((x * x) + (y * y));
//}


//void nrerror(char error_text[])
///* Numerical Recipes standard error handler */
//{
//printf("Numerical Recipes run-time error...\n");
//printf("%s\n",error_text);
//printf("...now exiting to system...\n");
//exit(1); 

 
/*Function PDcheck(evec)
'   Checks definiteness of symmetric matrices using their eigenvalues
'   Returns 1 (+ve def), 0.5 (+ve semi-def), -0.5 (-ve semi-def), -1 (-ve def)
    Dim pd, sa, smin, smax
    Dim i As Integer, n As Integer, p As Integer
    Dim svec() As Variant
    n = Application.Count(evec)
    ReDim svec(n)
    For i = 1 To n
        svec(i) = Sgn(evec(i))
    Next i
    sa = Application.sum(svec)
    smin = Application.Min(svec)
    smax = Application.Max(svec)
    If sa = n Then
        pd = 1
    ElseIf sa = -n Then
        pd = -1
    ElseIf sa >= 0 And smin >= 0 Then
        pd = 0.5
    ElseIf sa <= 0 And smax <= 0 Then
        pd = -0.5
   Else
        p = 0
    End If
    PDcheck = pd
End Function*/

#endif
