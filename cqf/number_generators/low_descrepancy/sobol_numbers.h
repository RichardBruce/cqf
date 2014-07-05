#ifndef __SOBOL_NUMBERS_H__
#define __SOBOL_NUMBERS_H__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <sstream>

#include <assert.h>

#include "atomic.h"

#include "number_generators_common.h"


template <class T>
class sobol_numbers
{
    public :
        /* CTOR */
        sobol_numbers(const char *const dir_file, const unsigned int n, const unsigned int d)
         : ref_cnt(new tbb::atomic<int>()), last_x(new unsigned int [d]), l((unsigned int)ceil(std::log(static_cast<T>(n))/log(2.0))), 
           cur_n(0), n(n), d(d)
        {
            memset(last_x, 0, d * sizeof(unsigned int));
            cache_params(dir_file);
        }
        
        /* Copy CTOR */
        sobol_numbers(const sobol_numbers &s)
         :  ref_cnt(s.ref_cnt), last_x(new unsigned int [s.d]), v(s.v), l(s.l), 
           cur_n(s.cur_n), n(s.n), d(s.d)
        {
            ref_cnt->fetch_and_add(1);
            memcpy(&last_x[0], &s.last_x[0], d * sizeof(unsigned int));
        }
        
        /* DTOR */
        ~sobol_numbers()
        {
            /* Release a reference to shared cache data c and v */
            const int refs = ref_cnt->fetch_and_add(-1);

            /* All references released */
            if (refs == 0)
            {
                delete [] v;
                delete ref_cnt;
            }            
            
            /* Delete non shared data */
            delete [] last_x;
        }
        
        
        /* Get the next number in all dimensions */
        void get_next(T *const numbers)
        {
            /* Check for overrun */
            assert(cur_n < n);
            
            /* First vector is all 0 */
            if (cur_n == 0)
            {
                memset(numbers, 0, d * sizeof(T));
                ++cur_n;
                return;
            }
            
            /* Only expect each number to be used once, find lowest 0 bit */
            unsigned int c = 0;
            unsigned int cur = cur_n - 1;
            while (cur & (0x1 << c))
            {
                ++c;
            }
            
            /* Create result and scale to -1.0 to 1.0 */
            unsigned int index = c * d;
            for (unsigned int i = 0; i < d; i++)
            {
                last_x[i]  ^= v[index++];
                numbers[i]  = static_cast<T>(last_x[i]) * pow_2_32_inv;
            }
            
            ++cur_n;
            return;
        }
        
        
        /* use the directly evaluated form for Gray coded Sobol numbers to move to 
         * any number in the sequence up to n */
        void skip(T *numbers, const unsigned int s)
        {
            /* Check for over run */
            assert(s < n);

            /* Skip to s */
            skip(s);
            
            /* Scale the output */
            for (unsigned int i = 0; i < d; i++)
            {
                numbers[i] = last_x[i] * pow_2_32_inv;
            }
        }
        
        void skip(const unsigned int s)
        {
            /* Check for over run */
            assert(s < n);

            const unsigned int skip_l   = static_cast<unsigned int>(std::log(static_cast<T>(s)) / std::log(2.0)) + 1;
            const unsigned int gray     = s ^ (s >> 1);
            for (unsigned int i = 0; i < d; i++)
            {
                last_x[i] = 0;
                for (unsigned int j = 0; j < skip_l; j++)
                {
                    if ((gray >> j) & 0x1)
                    {
                        last_x[i] ^= v[i + (j * d)];
                    }
                }
            }
            cur_n = s + 1;
        }
        

    private :
        /* Assignment operator */
        sobol_numbers& operator=(const sobol_numbers &s)
        {
        
        }
        
        void cache_params(const char *const dir_file)
        {
            /* Prepare the direction numbers */
            std::ifstream infile(dir_file, std::ios::in);
            assert(infile);
            
            /* First dimension */
            /* This may seem like a strange way to input the data, but it helps with reading */
            v = new unsigned int [l * d];
            for (unsigned int i = 0; i < l; i++)
            {
                v[i * d] = 1 << (31 - i);
            }
            
            /* Higher dimensions */
            for (unsigned int i = 1; i < d; i++)
            {
                unsigned int fake, s, a, m;
                infile >> fake >> s >> a;
            
                /* Direction numbers */
                if (l <= s)
                {
                    for (unsigned int j = 0; j < l; j++)
                    {
                        infile >> m;
                        v[i + (j * d)] = m << (31 - j);
                    }
                    
                    /* Read past unused entries */
                    for (unsigned int j = l; j < s; j++)
                    {
                        infile >> m;
                    }
                }
                else
                {
                    for (unsigned int j = 0; j < s; j++)
                    {
                        infile >> m;
                        v[i + (j * d)] = m << (31 - j);
                    }
         
                    for (unsigned int j = s; j < l; j++)
                    {
                        v[i + (j * d)] = v[i + ((j - s) * d)] ^ (v[i + ((j - s) * d)] >> s);
                        for (unsigned int k = 1; k < s; k++)
                        {
                            v[i + (j * d)] ^= (((a >> (s - 1 - k)) & 1) * v[i + ((j - k) * d)]);
                        }
                    }
                }
            }
            
            infile.close();
        }
        
        tbb::atomic<int>    *const  ref_cnt;       /* Initialised 0 */
        unsigned int        *const  last_x;
        unsigned int        *       v;
        unsigned int                l;
        unsigned int                cur_n;
        const unsigned int          n;
        const unsigned int          d;
};


#endif
