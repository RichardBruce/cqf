#ifndef __BOX_MULLER_H__
#define __BOX_MULLER_H__

#include <cmath>
#include <limits>
#include <iostream>

#include "constants.h"


class box_muller
{
    public :
        /* CTOR */
        box_muller(mersenne_twister *rand)
        {
            _rand = rand;
            _even_call = false;
        }
        
        /* DTOR */
        ~box_muller()
        {
            delete _rand;
        }
        
        float get_next()
        {
            if (!_even_call)
            {
                _even_call = true;
                
                /* Redraw until not 0 */
                float u0 = 0.0f;
                while (u0 == 0.0f)
                {
                    u0 = static_cast<float>(_rand->get_next()) / static_cast<float>(std::numeric_limits<int>::max());
                }
                const float u1 = static_cast<float>(_rand->get_next()) / static_cast<float>(std::numeric_limits<int>::max());
                
                _z1 = sqrt(-2.0f * std::log(u0)) * std::sin(2.0f * PI * u1);
                return sqrt(-2.0f * std::log(u0)) * std::cos(2.0f * PI * u1);
            }
            else
            {
                _even_call = false;
                return _z1;
            }            
        }
    
    private :
        mersenne_twister *_rand;
        float _z1;
        bool _even_call;
    
};

#endif
