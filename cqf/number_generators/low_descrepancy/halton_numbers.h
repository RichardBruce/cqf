#ifndef __HALTON_NUMBERS_H__
#define __HATLON_NUMBERS_H__

#include <cmath>
#include <iostream>
#include <assert.h>


const float primes[10] = { 2.0f, 3.0f, 5.0f, 7.0f, 11.0f, 13.0f, 17.0f, 19.0f, 23.0f, 29.0f };

class halton_numbers
{
    public :
        /* CTOR */
        halton_numbers(const unsigned int n, const unsigned int d)
         : n(n), d(d), numbers(new float[n * d])
        {
            /* Only allow dimensions up to 10 or there is too much correlation */
            assert(d < 11);
            
            prepare_numbers();
        }
        
        /* Copy CTOR */
        halton_numbers(const halton_numbers &h) : n(h.n), d(h.d), numbers(h.numbers)
        {
        
        }
        
        /* DTOR */
        ~halton_numbers()
        {
            delete [] numbers;
        }
        
        float get_number(const unsigned int i) const
        {
            return numbers[i];
        }
        
        float get_number(const unsigned int i, const unsigned int j) const
        {
            return get_number((i * d) + j);
        }
        
    private :
        /* Assignment operator */
        halton_numbers& operator=(const halton_numbers &)
        {
            return *this;
        }
        
        void prepare_numbers()
        {
            for (unsigned int i = 0; i < n; i++)
            {
                for (unsigned int j = 0; j < d; j++)
                {
                    numbers[(i * d) + j] = get_halton_number(i, primes[j]);
                }
            }
        }
        
        float get_halton_number(float index, const float base)
        {
            float halton_number = 0.0f;
            float f = 1.0f / base;
            
            while (index > 0.0f)
            {
                halton_number += f * fmod(index, base);
                index = floor(index / base);
                f /= base;
            }
            
            return halton_number;
        }
        
        const unsigned int  n;
        const unsigned int  d;
        float *const        numbers;
};

#endif
