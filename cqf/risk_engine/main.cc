#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <map>

#include <assert.h>

#include "pde_solver.h"
#include "crank_nicolson.h"
#include "risk_engine.h"

#include "cashflow_factory.h"

#include "scenario_point.h"
#include "risk_factory.h"

#include "grid_factory.h"

#include "utility.h"

//#include "correctors.h"

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
        return 0.0;
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
        return nullptr;
    }
    else
    {
        return argv[i];
    }
}


int main(int argc, char *argv[])
{
    /* Grid definitions */
    grid_factory<double> grid_parser;
    std::auto_ptr<grid<double>> grids(nullptr);
    
    /* Payoff */
    cashflow_factory<double> cashflow_parser;
    container_cleaner<std::vector<cashflow<double>*>> cash_cleaner;
    auto &cashflows = cash_cleaner.get();
    
    /* Risks to calculate */
    risk_factory<double> risk_parser;
    container_cleaner<std::vector<risk<double>*>> risk_cleaner;
    auto &risks = risk_cleaner.get();
    
    /* PDE solver */
    std::auto_ptr<pde_solver<double>> solver(nullptr);
    
    /* Null scenario */
    std::auto_ptr<scenario_point<double>> null_bump(nullptr);
    
    /* Economics */
    double spot = 100.0;
    double rate = 0.05;
    double vol = 0.2;
    
    /* PDE Solver Parameters */
    double max_dt = 0.01;
    //double s2_max = 300.0;
    //double s2_min = 0.0;
    //int nas2 = 600;


    /* Parse args */    
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-max_dt") == 0)
        {
            max_dt = get_arg_float("-max_dt", argv, argc, ++i);
        }
        /*else if (strcmp(argv[i], "-nas2") == 0)
        {
            nas2 = get_arg_int("-nas2", argv, argc, ++i);
        }*/
        /*else if (strcmp(argv[i], "-s2_max") == 0)
        {
            s2_max = get_arg_float("-s2_max", argv, argc, ++i);
        }*/
        /*else if (strcmp(argv[i], "-s2_min") == 0)
        {
            s2_min = get_arg_float("-s2_min", argv, argc, ++i);
        }*/
        /*else if (strcmp(argv[i], "-corr") == 0)
        {
            corr = get_arg_float("-corr", argv, argc, ++i);
        }*/
        else if (strcmp(argv[i], "-solver") == 0)
        {
            solver.reset(new crank_nicolson_pde_solver<double>());
        }
        else if (strcmp(argv[i], "-risk") == 0)
        {
            /* Find next arg */
            ++i;
            std::vector<std::string> risk_args;
            while ((i < argc) && (argv[i][0] != '-'))
            {
                risk_args.push_back(argv[i++]);
            }
            
            --i;
            risks.push_back(risk_parser.build(risk_args));
        }
        else if (strcmp(argv[i], "-cashflow") == 0)
        {
            assert((null_bump.get() == nullptr) || !"Error: All cashflows must be define before the economics.");
            assert((grids.get() == nullptr) || !"Error: All cashflows must be defined before the grids.");

            /* Find next arg */
            ++i;
            std::vector<std::string> cashflow_args;
            while ((i < argc) && (argv[i][0] != '-'))
            {
                cashflow_args.push_back(argv[i++]);
            }
            
            --i;
            cashflows.push_back(cashflow_parser.build(cashflow_args, spot));
        }
        else if (strcmp(argv[i], "-grid") == 0)
        {
            assert((cashflows.size() > 0) || !"Error: Cashflows must be defined before grids.");
            
            /* Find next arg */
            ++i;
            std::vector<std::string> grid_args;
            while ((i < argc) && (argv[i][0] != '-'))
            {
                grid_args.push_back(argv[i++]);
            }
            
            --i;
            grids.reset(grid_parser.build(grid_args, cashflows, spot));
        }
        else if (strcmp(argv[i], "-economics") == 0)
        {
            assert(cashflows.empty() || !"Error: Economics must be defined before cashflows.");
            assert((grids.get() == nullptr) || !"Error: Economics must be defined before grids.");
            spot = get_arg_float("-economics", argv, argc, ++i);
            rate = get_arg_float("-economics", argv, argc, ++i);
            vol = get_arg_float("-economics", argv, argc, ++i);
        }
        else
        {
            std::cout << "Error: Unknown arguement " << argv[i] << std::endl;
            assert(false);
        }
    }
    
    assert((cashflows.size() > 0) || !"Error: No cashflows specified.");
    assert((grids.get() != nullptr) || !"Error: No grids specified.");
    null_bump.reset(new scenario_point<double>(cashflows, (*grids), max_dt, spot, rate, vol));

    assert((solver.get() != nullptr) || !"Error: No solver specified.");
    assert((null_bump.get() != nullptr) || !"Error: No economics specified");
    assert((risks.size() > 0) || !"Error: No risks specified.");
    

    risk_engine<double> re(risks, *null_bump);
    re.run(*solver);
    
    return 0;
}
