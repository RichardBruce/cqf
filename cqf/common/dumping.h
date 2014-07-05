#ifndef __DUMPING_H__
#define __DUMPING_H__

#include <iostream>
#include <fstream>
#include <string>

#include <assert.h>


template<class T>
void dump_1d_data(std::ostream &s, const T *data, const int x)
{
    for (int i = 0; i < x; i++)
    {
        s << (*data) << ",";
        ++data;
    }
}


template<class T>
void dump_1d_data_to_file(const std::string &s, const T *data, const int x)
{
    if (!s.empty())
    {
        std::ofstream file(s.c_str());
        assert(file.is_open());
        dump_1d_data(file, data, x);
        file.close();
    }
}


template<class T>
void dump_2d_data(std::ostream &s, const T *data, const int x, const int y)
{
    for (int i = 0; i < y; i++)
    {
        for (int j = 0; j < x; j++)
        {
            s << (*data) << ",";
            ++data;
        }
        s << std::endl;
    }
}


template<class T>
void dump_2d_data_to_file(const std::string &s, const T *data, const int x, const int y)
{
    if (!s.empty())
    {
        std::ofstream file(s.c_str());
        assert(file.is_open());
        dump_2d_data(file, data, x, y);
        file.close();
    }
}


template<class T, class S>
void dump_row_to_gnuplot(std::ostream &s, const T *data, const int x, const int y, const S x_inc, const S y_inc)
{
    const S x_val = static_cast<S>(x) * x_inc;
    for (int i = 0; i < y; i++)
    {
        s << x_val << " " << (static_cast<S>(i) * y_inc) << " " << data[i] << std::endl;   
    }
    s << std::endl;
}


template<class T, class S>
void dump_2d_data_to_gnuplot(std::ostream &s, const T *data, const int x, const int y, const S x_inc, const S y_inc)
{
    for (int i = 0; i < x; i++)
    {
        dump_row_to_gnuplot(s, data, i, y, x_inc, y_inc);
    }
}


template<class T, class S>
void dump_2d_data_to_gnuplot(const std::string &s, const T *data, const int x, const int y, const S x_inc, const S y_inc)
{
    if (!s.empty())
    {
        std::ofstream file(s.c_str());
        assert(file.is_open());
        dump_2d_data_to_gnuplot(file, data, x, y, x_inc, y_inc);
        file.close();
    }
}


template<class T, class S>
void dump_row_to_gnuplot(std::ostream &s, const T *const data, const T *const ys, const T x_val, const S size)
{
    for (unsigned int i = 0; i < size; i++)
    {
        s << x_val << " " << ys[i] << " " << data[i] << std::endl;   
    }
    s << std::endl;
}


template<class T, class S>
void dump_row_to_gnuplot(std::ostream &s, const T *const data, const std::vector<T> &ys, const S x_val)
{
    dump_row_to_gnuplot(s, data, ys.data(), x_val, ys.size());
}


template<class T, class S>
void dump_row_to_gnuplot(std::ostream &s, const T *data, const typename std::vector<T>::const_iterator firsty, const typename std::vector<T>::const_iterator lasty, const S x_val)
{
    int data_idx = 0;
    for (typename std::vector<T>::const_iterator i = firsty; i != lasty; i++)
    {
        s << x_val << " " << (*i) << " " << data[data_idx++] << std::endl;   
    }
    s << std::endl;
}



#endif
