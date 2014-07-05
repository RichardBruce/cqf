#ifndef __EXTRAPOLATION_H__
#define __EXTRAPOLATION_H__

/* Look in 'in' for 'a' and then extrapolate from 'from' */
template<class T>
T linear_extrapolation(const T *const in, const T *const from, const int size, const T a)
{
    /* Find one past a */
    const T *const iter = std::lower_bound(in, &in[size], a);
    
    /* The index of one before a */
    const int v_idx = std::distance(in, iter) - 1;
    
    const T d = a - in[v_idx];
    const T g = (from[v_idx] - from[v_idx - 1]) / (in[v_idx] - in[v_idx - 1]);
    return from[v_idx] + (g * d);
}

#endif


