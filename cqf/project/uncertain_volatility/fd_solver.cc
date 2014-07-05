#include <algorithm>
#include <numeric>
#include <vector>

#include "fd_solver.h"

#include "option.h"

#include "matrix.h"

#include "utility.h"


double fd_solver::solve(const std::vector<option*> &options, const double s) const
{
//    for (unsigned int i = 0; i < options.size(); i++)
//    {
//        options[i]->dump();
//        std::cout << std::endl;
//    }
    
    /* Get the expiry dates */
    double max_strike = 0.0;
    std::vector<double> expiries;
    expiries.push_back(0.0);
    for (unsigned int i = 0; i < options.size(); i++)
    {
        max_strike = std::max(max_strike, options[i]->get_strike());
        expiries.push_back(options[i]->get_expiry());
    }
//    std::cout << "max strike: " << max_strike << std::endl;

    /* Calculate the time steps */
    std::sort(expiries.begin(), expiries.end());
    std::vector<double>::iterator expiries_iter = std::unique(expiries.begin(), expiries.end());
    expiries.resize(expiries_iter - expiries.begin());

    std::vector<double> padded_expiries;
    pad_data<double>(&padded_expiries, expiries, max_t_inc);
    std::vector<double> dt(padded_expiries.size());
    std::adjacent_difference(padded_expiries.begin() + 1, padded_expiries.end(), dt.begin());
//    for (unsigned int i = 0; i < padded_expiries.size(); i++)
//    {
//        std::cout << padded_expiries[i] << ", ";
//    }
//    std::cout << std::endl;
//    for (unsigned int i = 0; i < dt.size(); i++)
//    {
//        std::cout << dt[i] << ", ";
//    }
//    std::cout << std::endl;

        
    /* Grid set up */
    const int       nts         = dt.size();

    /* Make sure the asset steps land on the asset price */
    const int       nas         = 3 * ceil(s / max_s_inc);
    const double    s_inc       = s / ceil(s / max_s_inc);
//    std::cout << s_max << ", " <<  max_s_inc << std::endl;
//    std::cout << nts << ", " << nas << ", " << ngn << std::endl;
    double          *fd_grid    = new double[2 * nas];
    double          *payoff     = new double[nas];
//    std::cout << s_inc << std::endl;
    
    /* Add any end condition */
    memset(&fd_grid[nas], 0, nas * sizeof(double));
    for (unsigned int i = 0; i < options.size(); i++)
    {
        int payoff_offset = nas;
        if (options[i]->get_expiry() == padded_expiries.back())
        {
//            std::cout << "option " << i << std::endl;
            options[i]->payoff(&payoff[0], s_inc, nas);
            for (int j = 0; j < nas; j++)
            {
                fd_grid[payoff_offset++] += payoff[j];
            }
        }
    }
        
    /* Work backwards in time */
    const double s_inc_sq   = s_inc * s_inc;
    unsigned int grid_wr_idx = 0;
    unsigned int grid_rd_idx = nas;
    for (int i = nts - 2; i >= 0; i--)
    {
        /* Precomputes */
        const double t_inc              = dt[i];
//        std::cout << t_inc << std::endl;
        const double t_inc_s_inc_sq     = t_inc / (s_inc * s_inc);
        const double t_inc_s_inc        = t_inc / s_inc;
        
        const double theta      = 0.5 + (s_inc_sq / (12.0 * t_inc));
        const double m_theta    = 1.0 - theta;
//        std::cout << t_inc << ", " << theta << std::endl;
        assert(theta > 0.0);
        assert(theta < 1.0);

        /* Populate the matrix data */
        double *matrix_data_k    = new double [nas * nas];
        double *matrix_data_kp1  = new double [nas];
        memset(matrix_data_k, 0, (nas * nas * sizeof(double)));
        for (int j = 1; j < nas - 1; j++)
        {
            const double s      = static_cast<double>(j) * s_inc;
            const double s_sq   = s * s;
            
            double sigma        = min_vol;
            const double gamma  = (fd_grid[grid_rd_idx + j + 1] - (2.0 * fd_grid[grid_rd_idx + j]) + fd_grid[grid_rd_idx + j - 1]);
//            std::cout << i << " -> " << gamma << std::endl;
            if (gamma < 0.0)
            {
                sigma = max_vol;
            }
            const double half_sigma_sq  = sigma * sigma * 0.5;
        
            const double g = half_sigma_sq * s_sq;
            const double r_s = r * s;
            
            const double a = (t_inc_s_inc_sq * g) - (0.5 * t_inc_s_inc * r_s);
            const double b = (t_inc * -r) - (2.0 * t_inc_s_inc_sq * g);
            const double c = (t_inc_s_inc_sq * g) + (0.5 * t_inc_s_inc * r_s);
            
            matrix_data_k[(j * nas) + j - 1] = -a * theta;
            matrix_data_k[(j * nas) + j    ] = 1.0 - b * theta;
            matrix_data_k[(j * nas) + j + 1] = -c * theta;
            
            matrix_data_kp1[j]  = a * m_theta * fd_grid[grid_rd_idx + j - 1];
            matrix_data_kp1[j] += (1.0 + b * m_theta) * fd_grid[grid_rd_idx + j];
            matrix_data_kp1[j] += c * m_theta * fd_grid[grid_rd_idx + j + 1];
        }
        
        matrix_data_k[0]                    = exp(-r * t_inc);
        matrix_data_k[(nas * nas) - 1]      = 2.0;
        matrix_data_k[(nas * nas) - 2]      = -1.0;
        matrix_data_kp1[0]                  = exp(-r * t_inc) * fd_grid[grid_rd_idx];
        matrix_data_kp1[nas - 1]            = 2.0 * fd_grid[grid_rd_idx + nas - 1] - fd_grid[grid_rd_idx + nas - 2];
        
        matrix<double> matrix_kp1(&matrix_data_kp1[0], nas, 1, true);
        matrix<double> matrix_k(&matrix_data_k[0], nas, nas, true);
        matrix_k.gauss_solve(matrix_kp1, &fd_grid[grid_wr_idx], Tridiagonal);

        /* Add any payoff */
        for (unsigned int j = 0; j < options.size(); j++)
        {
            if (options[j]->get_expiry() == padded_expiries[i])
            {
                options[j]->payoff(&payoff[0], s_inc, nas);
//                std::cout << "payoff at: " << i << std::endl;
                for (int k = 0; k < nas; k++)
                {
                    fd_grid[grid_wr_idx + k] += payoff[k];
                }
            }
        }
        std::swap(grid_rd_idx, grid_wr_idx);
    }
    
    /* Dump the grid */
//    for (int i = 0; i < nas; i++)
//    {
//        std::cout << (i * s_inc) << ", ";
//        for (int j = 0; j < nts; j++)
//        {
//            std::cout << fd_grid[(j * nas) + i] << ", ";
//        }
//        std::cout << std::endl;
//    }
    
    
    /* Subtract the price of the vanilla (not the exotic) */
    double value = fd_grid[grid_rd_idx + static_cast<int>(s / s_inc)];
    for (unsigned int i = 0; i < options.size(); i++)
    {
        value -= options[i]->get_price();
    }
    
    /* Clean up */
    delete [] fd_grid;
    delete [] payoff;
 
//    std::cout << "Value: " << value << std::endl;       
    return value;
}
