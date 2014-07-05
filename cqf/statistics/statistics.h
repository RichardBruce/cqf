#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <vector>

#include "matrix.h"


template<class T>
matrix<T> build_covariance_matrix(const T *const diff, const int x, const int y)
{
    T *covar = new T[x * x];
    T x_inv = 1.0 / y;
    T x_m1_inv = 1.0 / (y - 1.0);
    
    /* For each column */
    for (int i = 0; i < x; i++)
    {
        /* Sum the current column */
        T sum_x = 0.0;
        for (int k = 0; k < y; k++)
        {
            sum_x += diff[(k * x) + i];
        }
        
        /* For each column including the outer loops */
        for (int j = 0; j < x; j++)
        {
            /* Take advantage of the problem's symetry */
            if (j < i)
            {
                covar[(i * x) + j] = covar[(j * x) + i];
            }
            else
            {
                /* Calculate the co-variance of these pair */
                T sum_y = 0.0;
                T sum_xy = 0.0;
                for (int k = 0; k < y; k++)
                {
                    sum_y += diff[(k * x) + j];
                    sum_xy += (diff[(k * x) + j] * diff[(k * x) + i]);
                }

                covar[(i * x) + j] = (sum_xy - (sum_x * sum_y * x_inv)) * x_m1_inv * (252.0 / 10000.0);
            }
        }
    }
    
    matrix<T> ret(covar, x, x);
    
    /* Clean up */
    delete [] covar;
    
    return ret;
}


template<class T>
T* build_difference_matrix(const std::vector<T> &data, const int x, const int y)
{
    T *diff = new T[x * (y - 1)];
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y - 1; j++)
        {
            diff[(j * x) + i] = data[((j + 1) * x) + i] - data[(j * x) + i];
        }
    }
    
    return diff;
}


template<class T>
T calculate_correlation(const T *const a, const T *const b, const int n)
{
    const T n_inv = 1.0 / static_cast<T>(n);
    const T n_m1_inv = 1.0 / static_cast<T>(n - 1);
    
    /* Sample averages */
    T a_avg = 0.0;
    T b_avg = 0.0;
    for (int i = 0; i < n; i++)
    {
        a_avg += a[i];
        b_avg += b[i];
    }
    a_avg *= n_inv;
    b_avg *= n_inv;
    
    /* Sample standard deviations */
    T a_sd = 0.0;
    T b_sd = 0.0;
    for (int i = 0; i < n; i++)
    {
        const T a_dif = a[i] - a_avg;
        const T b_dif = b[i] - b_avg;
        a_sd += a_dif * a_dif;
        b_sd += b_dif * b_dif;
    }
    a_sd = sqrt(a_sd * n_m1_inv);
    b_sd = sqrt(b_sd * n_m1_inv);
    
    /* Correlation */
    T covar = 0.0;
    for (int i = 0; i < n; i++)
    {
        covar += (a[i] - a_avg) * (b[i] - b_avg);
    }
    
    return (covar * n_m1_inv) / (a_sd * b_sd);
}


template<class T>
T calculate_rss(const T *const s, const T *const p, const int n)
{
    /*  Standard deviation of the sample vs. population error */
    T rss = 0.0;
    for (int i = 0; i < n; i++)
    {
        rss += (s[i] - p[i]) * (s[i] - p[i]);
    }
    
    return rss;
}


template<class T>
T calculate_standard_error(const T *const s, const T *const p, const int n)
{
    /* Residual sum of squares */
    const T rss = calculate_rss(s, p, n);

    /* Scaled to standard error */
    return sqrt(rss / static_cast<T>(n - 1));
}


template<class T>
T calculate_t_statistic(const T *const s, const T *const p, const T param, const T est, const int n)
{
    /* Standard error */
    const T se = calculate_standard_error(s, p, n);

    /* Test statistic */
    return (est - param) / se;
}


template<class T>
T calculate_f_statistic(const T *const s0, const T *const s1, const T *const p, const int p0, const int p1, const int n)
{
    /* Residual sum of squares */
    const T rss0 = calculate_rss(s0, p, n);
    const T rss1 = calculate_rss(s1, p, n);
    
    /* Test statistic */
    return ((rss0 - rss1) / (p1 - p0)) / (rss1 / (n - p1));
}

template<class T>
T calculate_auto_correlation(const T *const s, const int n, const int l)
{
    /* Whole sample average */
    T avg = 0.0;
    for (int i = 0; i < n; i++)
    {
        avg += s[i];
    }
    avg /= static_cast<T>(n);
    
    /* Whole sample var */
    T var = 0.0;
    for (int i = 0; i < n; i++)
    {
        const T diff = s[i] - avg;
        var += diff * diff;
    }
    var /= static_cast<T>(n);
    
    /* Sample - l covariance */
    T covar = 0.0;
    for (int i = l; i < n; i++)
    {
        covar += (s[i] - avg) * (s[i - l] - avg);
    }
    covar /= static_cast<T>(n - l);
    
    /* Auto correlation */
    return covar / var;
}


template<class T>
T calculate_auto_correlation_non_stationary(const T *const s, const int n, const int l)
{
    /* Sample - l average  */
    T avg_s = 0.0;
    T avg_l = 0.0;
    for (int i = l; i < n; i++)
    {
        avg_s += s[i - l];
        avg_l += s[i];
    }
    avg_s /= static_cast<T>(n - l);
    avg_l /= static_cast<T>(n - l);
    
    /* Sample - l var and covar */
    T var_s = 0.0;
    T var_l = 0.0;
    T covar = 0.0;
    for (int i = l; i < n; i++)
    {
        const T diff_s = s[i - l] - avg_s;
        const T diff_l = s[i] - avg_l;
        var_s += diff_s * diff_s;
        var_l += diff_l * diff_l;
        covar += (s[i] - avg_l) * (s[i - l] - avg_s);
    }
    var_s /= static_cast<T>(n - l);
    var_l /= static_cast<T>(n - l);
    covar /= static_cast<T>(n - l);
    
    /* Auto correlation */
    return covar / (sqrt(var_s) * sqrt(var_l));
}


#endif
