#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <vector>


template<class T>
void pad_data(std::vector<T> *padded_data, const std::vector<T> &required_data, const T max_step)
{
    /* For each step */
    for (unsigned int i = 0; i < required_data.size() - 1; i++)
    {
        /* Check it is small enough */
        const T interval = required_data[i + 1] - required_data[i];
        if (interval > max_step)
        {
            /* Pad with the minimum number of steps */
            const T nr_steps = ceil(interval / max_step);
            const T dt = interval / nr_steps;
            for (float j = 0.0; j < nr_steps; j++)
            {
                padded_data->push_back(required_data[i] + (j * dt));
            }
        }
        else
        {
            /* Step size is small enough */
            padded_data->push_back(required_data[i]);
        }
    }
    
    /* Add the last element */
    padded_data->push_back(required_data.back());
    
    return;
}


template<class T>
class container_cleaner
{
    public :
        explicit container_cleaner() {  };
        
        /* DTOR to clean up the contents of the container */
        ~container_cleaner()
        {
            for (typename T::iterator i = _container.begin(); i != _container.end(); i++)
            {
                delete (*i);
            }
        }
        
        /* Retreive the conatiner. */
        /* The purpose is not to encapsulate the container, but to clean its contents */
        T& get() { return _container; }
        
        const T& get() const { return _container; }

    private :
        T _container;
};

#endif
