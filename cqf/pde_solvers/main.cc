#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <map>

#include <assert.h>

#include "pde_solver.h"
#include "explicit.h"
#include "implicit.h"
#include "crank_nicolson.h"
#include "douglas.h"
#include "douglas_adi_2d.h"
#include "cs_adi_2d.h"
#include "mcs_adi_2d.h"
#include "hv_adi_2d.h"
#include "ade.h"

#include "cashflow.h"

#include "economics.h"
#include "scenario_point.h"
#include "risk.h"

#include "boundary_condition.h"

#include "grid_mapper.h"

#include "utility.h"


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


double get_arg_float(const std::string &option, const char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        help_missing_value(option);
        return 0;
    }
    else
    {
        return atof(argv[i]);
    }
}


char* get_arg_string(const std::string &option, char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        help_missing_value(option);
        return 0;
    }
    else
    {
        return argv[i];
    }
}


bump_dir_t get_bump_dir(const std::string &option, const char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        help_missing_value(option);
        return 0;
    }
    else
    {
        if (strcmp(argv[i], "DOWN") == 0)
        {
            return DOWN;
        }
        else if (strcmp(argv[i], "UP") == 0)
        {
            return UP;
        }
        else if (strcmp(argv[i], "BI") == 0)
        {
            return BI;
        }
        else
        {
            help_missing_value(option);
        }
    }
}


bump_style_t get_bump_style(const std::string &option, const char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        help_missing_value(option);
        return 0;
    }
    else
    {
        if (strcmp(argv[i], "ADD") == 0)
        {
            return ADD;
        }
        else if (strcmp(argv[i], "MULT") == 0)
        {
            return MULT;
        }
        else if (strcmp(argv[i], "OVERRIDE") == 0)
        {
            return OVERRIDE;
        }
        else
        {
            help_missing_value(option);
        }
    }
}


int main(int argc, char *argv[])
{
    /* Parse inputs */
    /* Solver */
    const char *solver_name = nullptr;
    
    /* Boundary Condition */
    std::auto_ptr<boundary_condition<double> > bc;
    
    /* Grid definitions */
    container_cleaner<std::vector<grid<double>*> > grid_cleaner;
    std::vector<grid<double>*> &grids = grid_cleaner.get();
    double grid_ratio   = 20.0;
    double grid_center  = 100.0;
    
    /* Payoff */
    container_cleaner<std::vector<cashflow<double>*> > cash_cleaner;
    std::vector<cashflow<double>*> &cashflows = cash_cleaner.get();
    
    /* Risks to calculate */
    container_cleaner<std::vector<risk<double>*> > risk_cleaner;
    std::vector<risk<double>*> &risks = risk_cleaner.get();
    
    
    /* PDE Solver Parameters */
    double max_dt = 0.01;
    double s1_max = 300.0;
    double s2_max = 300.0;
    double s1_min = 0.0;
    double s2_min = 0.0;
    int nas1 = 600;
    int nas2 = 600;
    int nts  = 100;
    
    /* Economics */
    container_cleaner<std::vector<equity_economics<double>*> > econ_cleaner;
    std::vector<equity_economics<double>*> &economics = econ_cleaner.get();
    double t    = 1.0;
    double corr = 0.0;
    
    /* Logging */
    std::string log_file;
    bool log_final_only = false;

    /* Parse args */    
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-max_dt") == 0)
        {
            max_dt = get_arg_float("-max_dt", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nas1") == 0)
        {
            nas1 = get_arg_int("-nas1", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nas2") == 0)
        {
            nas2 = get_arg_int("-nas2", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-s1_max") == 0)
        {
            s1_max = get_arg_float("-s1_max", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-s2_max") == 0)
        {
            s2_max = get_arg_float("-s2_max", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-s1_min") == 0)
        {
            s1_min = get_arg_float("-s1_min", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-s2_min") == 0)
        {
            s2_min = get_arg_float("-s2_min", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nts") == 0)
        {
            nts = get_arg_int("-nts", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-corr") == 0)
        {
            corr = get_arg_float("-corr", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-log_file") == 0)
        {
            log_file = get_arg_string("-log_file", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-log_final_only") == 0)
        {
            log_final_only = true;
        }
        else if (strcmp(argv[i], "-solver") == 0)
        {
            solver_name = argv[++i];
        }
        else if (strcmp(argv[i], "-risk") == 0)
        {
            assert(!grids.empty() || !"Error: Grid must be specified before risk.");
            
            double bump_size = 0.01;
            if (argv[i][0] != '-')
            {
                bump_size = get_arg_float("-risk", argv, argc, ++i);
            }
            
            bump_dir_t bump_dir = BI;
            if (argv[i][0] != '-')
            {
                bump_dir = get_bump_dir("-risk", argv, argc, ++i);
            }
            
            
            bump_style_t bump_style = ADD;
            if (argv[i][0] != '-')
            {
                bump_style = get_bump_style("-risk", argv, argc, ++i);
            }
            risks.add(new risk(grids[0], bump_size, bump_dir, bump_style));
        }
        else if (strcmp(argv[i], "-cashflow") == 0)
        {
            ++i;
            const double s1_inc = s1_max / static_cast<double>(nas1);
            const double s2_inc = s2_max / static_cast<double>(nas2);
            if (strcmp(argv[i], "wo_call") == 0)
            {
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                cashflows.push_back(new two_asset_payoff_cashflow<double, wo_call_cashflow_func<double> >(wo_call_cashflow_func<double>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2));
            }
            else if (strcmp(argv[i], "bo_call") == 0)
            {
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                cashflows.push_back(new two_asset_payoff_cashflow<double, bo_call_cashflow_func<double> >(bo_call_cashflow_func<double>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2));
            }
            else if (strcmp(argv[i], "wo_binary_call") == 0)
            {
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                cashflows.push_back(new two_asset_payoff_cashflow<double, wo_binary_call_cashflow_func<double> >(wo_binary_call_cashflow_func<double>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2));
            }
            else if (strcmp(argv[i], "call") == 0)
            {
                const bool   amer = (strcmp(argv[++i], "amer") == 0);
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat  = get_arg_float("-cashflow", argv, argc, ++i);
                if (amer)
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, amer_call_cashflow_func<double> >(amer_call_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
                else
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, euro_call_cashflow_func<double> >(euro_call_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
            }
            else if (strcmp(argv[i], "put") == 0)
            {
                const bool   amer = (strcmp(argv[++i], "amer") == 0);
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat  = get_arg_float("-cashflow", argv, argc, ++i);
                if (amer)
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, amer_put_cashflow_func<double> >(amer_put_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
                else
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, euro_put_cashflow_func<double> >(euro_put_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
            }
            else if (strcmp(argv[i], "binary_call") == 0)
            {
                const bool   amer = (strcmp(argv[++i], "amer") == 0);
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat  = get_arg_float("-cashflow", argv, argc, ++i);
                if (amer)
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, amer_binary_call_cashflow_func<double> >(amer_binary_call_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
                else
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, euro_binary_call_cashflow_func<double> >(euro_binary_call_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
            }
            else if (strcmp(argv[i], "binary_put") == 0)
            {
                const bool   amer = (strcmp(argv[++i], "amer") == 0);
                const double k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat  = get_arg_float("-cashflow", argv, argc, ++i);
                if (amer)
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, amer_binary_put_cashflow_func<double> >(amer_binary_put_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
                else
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, euro_binary_put_cashflow_func<double> >(euro_binary_put_cashflow_func<double>(k_s1), mat, s1_inc, nas1));
                }
            }
            else if (strcmp(argv[i], "dividend") == 0)
            {
                const bool   fix = (strcmp(argv[++i], "fix") == 0);
                const double div = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat = get_arg_float("-cashflow", argv, argc, ++i);
                if (fix)
                {
                    cashflows.push_back(new one_asset_dividend_cashflow<double, fixed_dividend<double> >(fixed_dividend<double>(div), mat));
                }
                else
                {
                    cashflows.push_back(new one_asset_dividend_cashflow<double, percent_dividend<double> >(percent_dividend<double>(div), mat));
                }
            }
            else if (strcmp(argv[i], "payment") == 0)
            {
                const bool   fix = (strcmp(argv[++i], "fix") == 0);
                const double div = get_arg_float("-cashflow", argv, argc, ++i);
                const double mat = get_arg_float("-cashflow", argv, argc, ++i);
                if (fix)
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, fixed_payment<double> >(fixed_payment<double>(div), mat, s1_inc, nas1));
                }
                else
                {
                    cashflows.push_back(new one_asset_payoff_cashflow<double, percent_payment<double> >(percent_payment<double>(div), mat, s1_inc, nas1));
                }
            }
            else
            {
                std::cout << "Error: Unknown cashflow " << argv[++i] << std::endl;
                assert(false);
            }
        }
        else if (strcmp(argv[i], "-economics") == 0)
        {
            const double s      = get_arg_float("-economics", argv, argc, i + 2);
            const double r      = get_arg_float("-economics", argv, argc, i + 3);
            const double sigma  = get_arg_float("-economics", argv, argc, i + 4);
            if (strcmp(argv[i + 1], "s1") == 0)
            {
                economics.push_back(new equity_economics<double>(s, r, sigma));
            }
            else if (strcmp(argv[i + 1], "s2") == 0)
            {
                economics.push_back(new equity_economics<double>(s, r, sigma));
            }
            else
            {
                std::cout << "Error: Unknown asset " << argv[i + 1] << std::endl;
                assert(false);
            }
            i += 4;
        }
        else if (strcmp(argv[i], "-bc") == 0)
        {
            const double max_s = get_arg_float("-bc", argv, argc, i + 2);
            const double min_s = get_arg_float("-bc", argv, argc, i + 3);
            if (strcmp(argv[i + 1], "const_ds") == 0)
            {
                bc.reset(new const_ds_boundary_condition<double>(max_s, min_s));
            }
            else if (strcmp(argv[i + 1], "barrier") == 0)
            {
                bc.reset(new barrier_boundary_condition<double>(max_s, min_s));
            }
            else
            {
                std::cout << "Error: Unknown boundary condition " << argv[i + 1] << std::endl;
                assert(false);
            }
            i += 3;
        }
        else if (strcmp(argv[i], "-grid") == 0)
        {
            const char* grid_name = argv[++i];
            if (strcmp(grid_name, "uniform") == 0)
            {
                grids.push_back(new grid<double>(new uniformed_grid_mapper<double>(s1_max, s1_min, nas1)));
            }
            else if (strcmp(grid_name, "sinh") == 0)
            {
                const T grid_ratio  = get_arg_float("-grid", argv, argc, ++i);
                const T grid_center = get_arg_float("-grid", argv, argc, ++i);
                grids.push_back(new grid<double>(new sinh_grid_mapper<double>(grid_center, grid_ratio, s1_max, s1_min, nas1)));
            }
            else
            {
                std::cout << "Error: Unknown grid type " << grid_names[i] << std::endl;
                assert(false);
            }
        }
        else
        {
            std::cout << "Error: Unknown arguement " << argv[i] << std::endl;
            assert(false);
        }
    }
    
    
    /* Build the solver */
    //assert(bc.get() != nullptr);
    assert(economics.size());
    assert(cashflows.size());
    std::auto_ptr<pde_solver<double> > solver;
    if (strcmp(solver_name, "explicit") == 0)
    {
        solver.reset(new explicit_pde_solver<double>(cashflows, economics, grids, (*bc.get()), log_file, max_dt, log_final_only));
    }
    else if (strcmp(solver_name, "implicit") == 0)
    {
        solver.reset(new implicit_pde_solver<double>(cashflows, economics, grids, (*bc.get()), log_file, max_dt, log_final_only));
    }
    else if (strcmp(solver_name, "cn") == 0)
    {
        solver.reset(new crank_nicolson_pde_solver<double>(cashflows, economics, grids, (*bc.get()), log_file, max_dt, log_final_only));
    }
    else if (strcmp(solver_name, "doug") == 0)
    {
        solver.reset(new douglas_pde_solver<double>(cashflows, economics, grids, (*bc.get()), log_file, max_dt, log_final_only));
    }
    else if (strcmp(solver_name, "doug_adi") == 0)
    {
        solver.reset(new douglas_adi_2d_pde_solver<double>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only));
    }
    else if (strcmp(solver_name, "cs") == 0)
    {
        solver.reset(new cs_adi_2d_pde_solver<double>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only));
    }
    else if (strcmp(solver_name, "mcs") == 0)
    {
        solver.reset(new mcs_adi_2d_pde_solver<double>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only));
    }
    else if (strcmp(solver_name, "hv") == 0)
    {
        solver.reset(new hv_adi_2d_pde_solver<double>(cashflows, economics, log_file, corr, nts, nas1, nas2, t, log_final_only));
    }
    else if (strcmp(solver_name, "ade") == 0)
    {
        solver.reset(new ade_pde_solver<double>(cashflows, economics, grids, (*bc.get()), log_file, max_dt, log_final_only));
    }
    else
    {
        std::cout << "Error: Unknown solver " << solver_name << std::endl;
        assert(false);
    }
    
    std::cout << "Number of time steps " << nts << std::endl;
    std::cout << "Number of S1 asset steps " << nas1 << std::endl;
    std::cout << "Number of S2 asset steps " << nas2 << std::endl;
    std::cout << "Option lifetime " << t << std::endl;
    
    /* Run the solver */
    assert(solver.get() != nullptr);
    const double result = solver->solve();
    std::cout << "Result: " << result << std::endl;
    
    return 0;
}
