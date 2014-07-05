#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <assert.h>


void help_missing_value(const std::string &option)
{
    std::cout << "Warning: No value for option " << option << std::endl;
}


int get_arg_int(const std::string &option, const char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        help_missing_value(option);
        return 0;
    }
    else
    {
        return atoi(argv[i]);
    }
}

int main(int argc, char *argv[])
{
    /* Parse inputs */
    int nas1 = 600;
    int nas2 = 600;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-nas1") == 0)
        {
            nas1 = get_arg_int("-nas1", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nas2") == 0)
        {
            nas2 = get_arg_int("-nas2", argv, argc, ++i);
        }
        else
        {
            std::cout << "Unknown arguement" << std::endl;
            assert(false);
        }
    }
    
    /* Economics */
    const double t                  = 1.0;
    const double k_s1               = 100.0;
    const double r_s1               = 0.05;
    const double sigma_s1           = 0.2;
    const double half_sigma_s1_sq   = 0.5 * sigma_s1 * sigma_s1;
    
    const double k_s2               = 100.0;
    const double r_s2               = 0.05;
    const double sigma_s2           = 0.2;
    const double half_sigma_s2_sq   = 0.5 * sigma_s2 * sigma_s2;
    
    const double corr_s1_s2         = 0.5;
    
    /* Grid set up */
    /* TODO -- What is the limit on the time step */
    double t_inc        = 0.05 / (static_cast<double>(nas1 * nas1) * sigma_s1 * sigma_s1);
    const int nts       = (int)(t / t_inc) + 1; /* Number of time steps */
    t_inc               = t / nts;
    const int nsn       = nas1 * nas2;
    double *fd_grid     = new double[nsn * 2];
    
    const double s1_max          = 300.0;
    const double s1_inc          = s1_max / static_cast<double>(nas1);
    const double s1_inc_inv      = 1.0 / s1_inc;
    const double s1_inc2_inv     = 0.5 * s1_inc_inv;
    const double s1_inc_sq_inv   = 1.0 / (s1_inc * s1_inc);
    
    const double s2_max          = 300.0;
    const double s2_inc          = s2_max / static_cast<double>(nas2);
    const double s2_inc_inv      = 1.0 / s2_inc;
    const double s2_inc2_inv     = 0.5 * s2_inc_inv;
    const double s2_inc_sq_inv   = 1.0 / (s2_inc * s2_inc);
        
    std::cout << "Number of time steps " << nts << std::endl;
    std::cout << "Number of S1 asset steps " << nas1 << std::endl;
    std::cout << "Number of S2 asset steps " << nas2 << std::endl;

    /* Set the end condition (ie payoff) */
    int payoff_offset = nsn;
    for (int i = 0; i < nas1; i++)
    {
        const double s1 = s1_inc * static_cast<double>(i);
        for (int j = 0; j < nas2; j++)
        {
            /* Call payoff */
//            fd_grid[payoff_offset++] = std::max((s1_inc * static_cast<double>(i)) - k_s1, 0.0f);

            /* Cash or nothing binary call */
            if ((s1 > k_s1) && ((s2_inc * static_cast<double>(j)) > k_s2))
            {
                fd_grid[payoff_offset++] = 1.0;
            }
            else
            {
                fd_grid[payoff_offset++] = 0.0;
            }
        }
    }

    
    /* Work backwards in time */
    int t_offset      = 0; 
    int tp1_offset    = nsn;
    for (int i = nts - 2; i >= 0; i--)
    {
        /* Fill in the row */
        for (int j = 1; j < nas1 - 1; j++)
        {
            const int s1_idx        = j * nas2;
            const int s1_idx_p1     = s1_idx + nas2;
            const int s1_idx_m1     = s1_idx - nas2;
            
            /* Fill in the column */
            for (int k = 1; k < nas2 - 1; k++)
            {
                /* S1 greeks */
                const double s1         = static_cast<double>(j) * s1_inc;
                const double delta_s1   = (fd_grid[tp1_offset + s1_idx_p1 + k] - fd_grid[tp1_offset + s1_idx_m1 + k]) * s1_inc2_inv;
                const double gamma_s1   = (fd_grid[tp1_offset + s1_idx_p1 + k] - (2.0 * fd_grid[tp1_offset + s1_idx + k]) + fd_grid[tp1_offset + s1_idx_m1 + k]) * s1_inc_sq_inv;
                const double rv_s1      = r_s1 * fd_grid[tp1_offset + s1_idx + k];

                /* s2 greeks */
                const double s2         = static_cast<double>(k) * s2_inc;
                const double delta_s2   = (fd_grid[tp1_offset + s1_idx + k + 1] - fd_grid[tp1_offset + s1_idx + k - 1]) * s2_inc2_inv;
                const double gamma_s2   = (fd_grid[tp1_offset + s1_idx + k + 1] - (2.0 * fd_grid[tp1_offset + s1_idx + k]) + fd_grid[tp1_offset + s1_idx + k - 1]) * s2_inc_sq_inv;
                const double rv_s2      = r_s2 * fd_grid[tp1_offset + s1_idx + k];
                
                /* Cross partial */
                const double delta_m1   = (fd_grid[tp1_offset + s1_idx_m1 + k + 1] - fd_grid[tp1_offset + s1_idx_m1 + k - 1]) * s2_inc2_inv;
                const double delta_p1   = (fd_grid[tp1_offset + s1_idx_p1 + k + 1] - fd_grid[tp1_offset + s1_idx_p1 + k - 1]) * s2_inc2_inv;
                const double x_gamma    = (delta_p1 - delta_m1) * s1_inc2_inv;

                /* Update */
                const double theta      = (half_sigma_s1_sq * s1 * s1 * gamma_s1) + (r_s1 * s1 * delta_s1) - rv_s1 +
                                          (half_sigma_s2_sq * s2 * s2 * gamma_s2) + (r_s2 * s2 * delta_s2) - rv_s2 +
                                          (sigma_s1 * s1 * sigma_s2 * s2 * corr_s1_s2 * x_gamma);
                fd_grid[t_offset + s1_idx + k]   = (theta * t_inc) + fd_grid[tp1_offset + s1_idx + k];
            }
    
            /* Apply the s2 = 0 boundary condition */
            fd_grid[t_offset + s1_idx] = (2.0 * fd_grid[t_offset + s1_idx + 1]) - fd_grid[t_offset + s1_idx + 2];
        
            /* Apply the s2 = s2_max boundary condition */
            fd_grid[t_offset + s1_idx_p1 - 1] = (2.0 * fd_grid[t_offset + s1_idx_p1 - 2]) - fd_grid[t_offset + s1_idx_p1 - 3];
        }
    
        /* Upper and lower boundary condition in s1 */
        const int s1_idx_upper  = (nas1 - 1) * nas2;
        const int s1_idx_m1     = s1_idx_upper - nas2;
        const int s1_idx_m2     = s1_idx_m1 - nas2;
        
        const int s1_idx_lower  = 0;
        const int s1_idx_p1     = nas2;
        const int s1_idx_p2     = (nas2 << 1);
        for (int j = 1; j < nas2 - 1; j++)
        {
            fd_grid[t_offset + s1_idx_upper + j] = (2.0 * fd_grid[t_offset + s1_idx_m1 + j]) - fd_grid[t_offset + s1_idx_m2 + j];
            fd_grid[t_offset + s1_idx_lower + j] = (2.0 * fd_grid[t_offset + s1_idx_p1 + j]) - fd_grid[t_offset + s1_idx_p2 + j];
        }
        
        /* The corners */
        /* Apply the s1 = 0 , s2 = 0 boundary condition */
        fd_grid[t_offset] = (2.0 * fd_grid[t_offset + 1]) - fd_grid[t_offset + 2];
        
        /* Apply the s1 = 0, s2 = s2_max boundary condition */
        fd_grid[t_offset + nas2 - 1] = (2.0 * fd_grid[t_offset + nas2 - 2]) - fd_grid[t_offset + nas2 - 3];
        
        /* Apply the s1 = s1_max, s2 = 0 boundary condition */
        fd_grid[t_offset + s1_idx_upper] = (2.0 * fd_grid[t_offset + s1_idx_upper + 1]) - fd_grid[t_offset + s1_idx_upper + 2];

        /* Apply the s1 = s1_max, s2 = s2_max boundary condition */
        fd_grid[t_offset + nsn - 1] = (2.0 * fd_grid[t_offset + nsn - 2]) - fd_grid[t_offset + nsn - 3];

        std::swap(tp1_offset, t_offset);
    }
    
    
    /* Dump the grid */
    for (int i = 0; i < nas1; i++)
    {
        for (int j = 0; j < nas2; j++)
        {
            std::cout << (i * s1_inc) << " " << (j * s2_inc) << " " << fd_grid[(i * nas1) + j] << " " << std::endl;   
        }
        std::cout << std::endl;
    }
    std::cout << fd_grid[static_cast<int>(((1.0/3.0) * nas1 * nas2) + ((1.0/3.0) * nas2))] << std::endl;
    
    /* Clean up */
    delete [] fd_grid;
    
    return 0;
}
