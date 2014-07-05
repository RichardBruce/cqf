#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>

#include <assert.h>


#include "matrix.h"


void solve()
{    
    /* Economics */
    const double k       = 100.0;
    const double r       = 0.05;
    const double t       = 1.0;
    const double sigma   = 0.2;
    
    /* Grid set up */
    const int ngn = nas * nts;
    
    const double s_max = 300.0;
    const double s_inc = s_max / static_cast<double>(nas);
    
    double *fd_grid = new double[ngn];
    
    /* Pre-computes */
    const double s_inc_sq           = s_inc * s_inc;
    const double half_sigma_sq      = 0.5 * sigma * sigma;
    const double t_inc              = t / static_cast<double>(nts);
    const double t_inc_s_inc_sq     = t_inc / s_inc_sq;
    const double t_inc_s_inc        = t_inc / s_inc;
    const double dbl_s_inc_inv      = 1.0 / (2.0 * s_inc);
    const double s_inc_sq_inv       = 1.0 / s_inc_sq;
    
    
    /* Set the end condition (ie payoff) */
    int payoff_offset = ngn - nas;
    for (int i = 0; i < nas; i++)
    {
        /* Call payoff */
//        fd_grid[payoff_offset++] = std::max((s_inc * static_cast<double>(i)) - k, 0.0f);

        /* Cash or nothing binary call */
        if ((s_inc * static_cast<double>(i)) > k)
        {
            fd_grid[payoff_offset++] = 1.0;
        }
        else
        {
            fd_grid[payoff_offset++] = 0.0;
        }
    }
    
    /* Work backwards in time */
    double theta      = 0.5 + (s_inc_sq / (12.0 * t_inc));
    double m_theta    = 1.0 - theta;
    std::cout << "Theta is " << theta << std::endl;
    assert(theta > 0.0);
    assert(theta < 1.0);
    for (int i = nts - 2; i >= 0; i--)
    {
        double *matrix_data_k    = new double [nas * nas];
        double *matrix_data_kp1  = new double [nas];
        memset(matrix_data_k, 0, (nas * nas * sizeof(double)));
        
        if ((i == nts - 2) || (i == nts - 3))
        {
            theta = 1.0;
            m_theta = 0.0;
        }
        else
        {
            theta      = 0.5 + (s_inc_sq / (12.0 * t_inc));
            m_theta    = 1.0 - theta;
        }
    
        /* Populate the matrix */
        for (int j = 1; j < nas - 1; j++)
        {
            const double s      = static_cast<double>(j) * s_inc;
            const double s_sq   = s * s;
            
            const double g = half_sigma_sq * s_sq;
            const double r_s = r * s;
            
            const double a = (t_inc_s_inc_sq * g) - (0.5 * t_inc_s_inc * r_s);
            const double b = (t_inc * -r) - (2.0 * t_inc_s_inc_sq * g);
            const double c = (t_inc_s_inc_sq * g) + (0.5 * t_inc_s_inc * r_s);
            
            matrix_data_k[(j * nas) + j - 1] = -a * theta;
            matrix_data_k[(j * nas) + j    ] = 1.0 - b * theta;
            matrix_data_k[(j * nas) + j + 1] = -c * theta;
            
            matrix_data_kp1[j]  = a * m_theta * fd_grid[(i * nas) + nas + j - 1];
            matrix_data_kp1[j] += (1.0 + b * m_theta) * fd_grid[(i * nas) + nas + j];
            matrix_data_kp1[j] += c * m_theta * fd_grid[(i * nas) + nas + j + 1];
  
//            const double delta = (fd_grid[(i * nas) + nas + j + 1] - fd_grid[(i * nas) + nas + j - 1]) * dbl_s_inc_inv;
//            std::cout << (j * s_inc) << " " << (i * t_inc) << " " << delta << std::endl;
  
//            const double gamma = (fd_grid[(i * nas) + nas + j + 1] - (2.0 * fd_grid[(i * nas) + nas + j]) + fd_grid[(i * nas) + nas + j - 1]) * s_inc_sq_inv;
//            std::cout << (j * s_inc) << " " << (i * t_inc) << " " << gamma << std::endl;
        }
//        std::cout << std::endl;
        
        matrix_data_k[0]                    = exp(-r * t_inc);
        matrix_data_k[(nas * nas) - 1]      = 2.0;
        matrix_data_k[(nas * nas) - 2]      = -1.0;
        matrix_data_kp1[0]                  = exp(-r * t_inc) * fd_grid[(i * nas) + nas];
        matrix_data_kp1[nas - 1]            = 2.0 * fd_grid[(i * nas) + nas + nas - 1] - fd_grid[(i * nas ) + nas + nas - 2];
        
        matrix<double> matrix_kp1(&matrix_data_kp1[0], nas, 1, true);
        matrix<double> matrix_k(&matrix_data_k[0], nas, nas, true);
        matrix_k.gauss_solve(matrix_kp1, &fd_grid[i * nas], Tridiagonal);
    }
    
    
    /* Dump the grid */
    for (int j = 0; j < nts; j++)
    {
        //std::cout << (i * s_inc) << ", ";
        for (int i = 0; i < nas; i++)
        {
            std::cout << (i * s_inc) << " " << (j * t_inc) << " " << fd_grid[(j * nas) + i] << std::endl;
        }
        std::cout << std::endl;
    }
    
    
    /* Clean up */
    delete [] fd_grid;
    
    return;
}
