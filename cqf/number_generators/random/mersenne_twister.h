#ifndef __MERSENNE_TWISTER_H__
#define __MERSENNE_TWISTER_H__

#include <cstring>

#include <assert.h>

#include "number_generators_common.h"


/* Mersenne twister random number generator */
class mersenne_twister
{
    public :
        /* CTOR */
        mersenne_twister(const unsigned int seed = 846234)
        {
            _mt = new unsigned int[624];
            _mt[0] = seed;
            for (int i = 1; i < 624; i++)
            {
	            _mt[i] = (0x6c078965 * (_mt[i-1] ^ (_mt[i-1] >> 30))) + i;
            }
        
            _index = 624;
        }
    
        /* Copy CTOR */
        mersenne_twister(mersenne_twister &m)
        {
            _index = m._index;
            _mt = new unsigned int[624];
            memcpy(_mt, m._mt, 624 * sizeof(unsigned int));
        }
    
        /* DTOR */
        ~mersenne_twister()
        {
            delete [] _mt;
        }
        
        /* Get a random number [-1,1] */
        int get_next()
        {
            /* Refill the array */
            if (_index == 624) 
            {
                static unsigned long call = 0;
                ++call;
                
                _index = 0;
                static const unsigned int eor_mask[2] = { 0x00000000, 0x9908b0df };
                for (int i = 0; i < 227; i++)
                {
                    unsigned int y = (_mt[i] & 0x80000000) | (_mt[i + 1] & 0x7fffffff);
                    _mt[i] = (_mt[i + 397] ^ (y >> 1)) ^ eor_mask[y & 0x1];
                }
     
                for (int i = 227; i < 623; i++)
                {
                    unsigned int y = (_mt[i] & 0x80000000) | ((_mt[i + 1]) & 0x7fffffff);
                    _mt[i] = (_mt[i - 227] ^ (y >> 1)) ^ eor_mask[y & 0x1];
                }

                unsigned int y = (_mt[623] & 0x80000000) | ((_mt[0]) & 0x7fffffff);
                _mt[623] = (_mt[623 - 227] ^ (y >> 1)) ^ eor_mask[y & 0x1];
            }
     
            /* Temper and return number */
            unsigned int y = _mt[_index++];
            y = y ^ y >> 11;
            y = y ^ ((y << 7) & 0x9d2c5680); 
            y = y ^ ((y << 15) & 0xefc60000); 
            y = y ^ (y >> 18);
     
            return y;
        }
        
        double get_next_prob()
        {
            unsigned int rand;
            do 
            {
                rand = this->get_next();
            } while (rand == 0);
            
            double ret = static_cast<double>(rand) * pow_2_32_inv;
            assert(ret < 1.0);
            assert(ret > 0.0);
            return ret;
        }

    private :
        unsigned int  _index;
        unsigned int *_mt;
};

#endif
