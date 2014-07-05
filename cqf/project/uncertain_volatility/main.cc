#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <assert.h>

#include "fd_solver.h"

#include "downhill_simplex.h"
#include "quasi_newton_solver.h"
#include "metropolis.h"

#include "option.h"
#include "call_option.h"
#include "put_option.h"
#include "binary_call_option.h"
#include "binary_put_option.h"


enum optimiser_algo_t { DOWNHILL_SIMPLEX = 0, QUASI_NEWTON = 1, METROPOLIS = 2 };

void usage()
{
    std::cout << "Usage:"                                                                   << std::endl;
    std::cout << "uncertainty_tool"                                                         << std::endl;
    std::cout << "Economics"                                                                << std::endl;
    std::cout << "  -stock_price d          The current stock price"                        << std::endl;
    std::cout << "  -interest_rate d        The current risk free rate"                     << std::endl;
    std::cout << "  -max_vol d              The maximum bound for volatility"               << std::endl;
    std::cout << "  -min_vol d              The minimum bound for volatility"               << std::endl;
    std::cout << "  -hedge opt k t n b a    An option to hedge with"                        << std::endl;
    std::cout << "  -target opt k t n b a   The target option"                              << std::endl;
    std::cout << "Expected values for opt are call, put, binary_call and binary_put"        << std::endl;
    std::cout << "s is the option strike"                                                   << std::endl;
    std::cout << "t is the option tenor"                                                    << std::endl;
    std::cout << "n is the option notional"                                                 << std::endl;
    std::cout << "b is the option bid"                                                      << std::endl;
    std::cout << "a is the option ask"                                                      << std::endl;
    std::cout << "0 bid and ask may be entered for the target option"                       << std::endl;
    std::cout << "The last specified target option will be the one that is priced"          << std::endl;
    std::cout << std::endl;
    std::cout << "Finite Difference Control"                                                << std::endl;
    std::cout << "  -max_t_inc d            Maximum time increment for finite difference"   << std::endl;
    std::cout << "  -max_s_inc d            Maximum asset increment for finite difference"  << std::endl;
    std::cout << "  -price_only             Only price the portfolio, dont optimise"        << std::endl;
    std::cout << "If not provided suitable defaults will be used"                           << std::endl;
    std::cout << std::endl;
    std::cout << "Optimiser Control"                                                        << std::endl;
    std::cout << "  -optimiser_algo s       The optimiser algorithm to use"                 << std::endl;
    std::cout << "  -opt_tolerance d        The tolerance of the optimiser"                 << std::endl;
    std::cout << "  -opt_max_iter d         The maximum number of optimisation iterations"  << std::endl;
    std::cout << "  -opt_annealing_period d Factor of iterations before quenching"          << std::endl;
    std::cout << "Optimiser algoritms may be one of"                                        << std::endl;
    std::cout << "  downhill_simplex"                                                       << std::endl;
    std::cout << "  quasi_newton"                                                           << std::endl;
    std::cout << "  metropolis"                                                             << std::endl;
    std::cout << "-opt_annealing_period is only used with the Metropolis algorithm"         << std::endl;
    std::cout << "If not provided suitable defaults will be used"                           << std::endl;
}

void help_missing_value(const std::string &option)
{
    std::cout << "Warning: No value for option " << option << std::endl;
    usage();
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
        return NULL;
    }
    else
    {
        return argv[i];
    }
}

option* get_arg_option(const std::string &option, char *const argv[], const int argc, const int i)
{
    if ((i + 5) >= argc)
    {
        help_missing_value(option);
        return NULL;
    }

    /* Parse the option details */
    const double t = get_arg_float("option tenor", argv, argc, i + 1);
    const double k = get_arg_float("option strike", argv, argc, i + 2);
    const double n = -get_arg_float("option notional", argv, argc, i + 3);
    const double b = get_arg_float("option bid", argv, argc, i + 4);
    const double a = get_arg_float("option ask", argv, argc, i + 5);
    if (strcmp(argv[i], "call") == 0)
    {
        return new call_option(t, k, n, b, a);   
    }
    else if (strcmp(argv[i], "put") == 0)
    {
        return new put_option(t, k, n, b, a);
    }
    else if (strcmp(argv[i], "binary_call") == 0)
    {
        return new binary_call_option(t, k, n, b, a);
    }
    else if (strcmp(argv[i], "binary_put") == 0)
    {
        return new binary_put_option(t, k, n, b, a);
    }
    else
    {
        std::cout << "Warning: Unrecognised option " << argv[i] << std::endl;
        std::cout << "Expected one of:"                         << std::endl;
        std::cout << "  call"                                   << std::endl;
        std::cout << "  put"                                    << std::endl;
        std::cout << "  binary_call"                            << std::endl;
        std::cout << "  binary_put"                             << std::endl;
    }
    
    return NULL;
}

int main(int argc, char *argv[])
{
    /* Parse inputs */
    std::vector<option*>    hedge;
    optimiser_algo_t        optimiser_algo = METROPOLIS;
    option* target          = NULL;
    double s                = 100.0;
    double r                = 0.06;
    double max_vol          = 0.23;
    double min_vol          = 0.17;
    double max_s_inc        = 0.1;
    double max_t_inc        = 0.025;
    double tolerance        = 0.0000001;
    double anlg_period      = 0.5;
    int    max_iter         = 500;
    bool   price_only       = false;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-max_s_inc") == 0)
        {
            max_s_inc = get_arg_float("-max_s_inc", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-max_t_inc") == 0)
        {
            max_t_inc = get_arg_float("-max_t_inc", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-stock_price") == 0)
        {
            s = get_arg_float("-stock_price", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-interest_rate") == 0)
        {
            r = get_arg_float("-interest_rate", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-max_vol") == 0)
        {
            max_vol = get_arg_float("-max_vol", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-min_vol") == 0)
        {
            min_vol = get_arg_float("-min_vol", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-price_only") == 0)
        {
            price_only = true;
        }
        else if (strcmp(argv[i], "-hedge") == 0)
        {
            hedge.push_back(get_arg_option("-hedge", argv, argc, ++i));
            i += 5;
        }
        else if (strcmp(argv[i], "-target") == 0)
        {
            target = get_arg_option("-target", argv, argc, ++i);
            i += 5;
        }
        else if (strcmp(argv[i], "-opt_tolerance") == 0)
        {
            tolerance = get_arg_float("-opt_tolerance", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-opt_max_iter") == 0)
        {
            max_iter = get_arg_int("-opt_max_iter", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-opt_annealing_period") == 0)
        {
            anlg_period = get_arg_float("-opt_annealing_period", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-optimiser_algo") == 0)
        {
            ++i;
            if (strcmp(argv[i], "downhill_simplex") == 0)
            {
                optimiser_algo = DOWNHILL_SIMPLEX;
            }
            else if (strcmp(argv[i], "quasi_newton") == 0)
            {
                optimiser_algo = QUASI_NEWTON;
            }
            else if (strcmp(argv[i], "metropolis") == 0)
            {
                optimiser_algo = METROPOLIS;
            }
            else
            {
                std::cout << "Warning: Unknown optimiser algorithm " << argv[i] << std::endl;
                std::cout << "Expected one of:"         << std::endl;
                std::cout << "  downhill_simplex"       << std::endl;
                std::cout << "  quasi_newton"           << std::endl;
                std::cout << "  metropolis"             << std::endl;
                usage();
            }
        }
        else if (strcmp(argv[i], "-help") == 0)
        {
            usage();
        }
        else
        {
            std::cout << "Warning: Unknown option " << argv[i] << std::endl;
            usage();
        }
    }
    

    /* Initialise the finite difference solver */
    fd_solver fd_sol(max_vol, min_vol, r, max_t_inc, max_s_inc);
    
    /* Just pricing */
    if (price_only)
    {
        if (target != NULL)
        {
            hedge.push_back(target);
        }
        assert(hedge.size() > 0);
        
        const double val = fd_sol.solve(hedge, s);
        std::cout << "Value: " << val << std::endl;
    }
    /* Optimising */
    else
    {
        assert(target != NULL);
        assert(hedge.size() > 0);

        /* Pick the selected optimiser */
        optimiser *opt;
        switch (optimiser_algo)
        {
            case DOWNHILL_SIMPLEX   :
                opt = new downhill_simplex(fd_sol, tolerance, max_iter);
                break;
            case QUASI_NEWTON       :
                opt = new quasi_newton_solver(fd_sol, tolerance, max_iter);
                break;
            case METROPOLIS         :
                opt = new metropolis(fd_sol, tolerance, anlg_period, max_iter);
                break;
            default :
                std::cout << "Error: This isnt possible, unknown optimiser selected" << std::endl;
                return 1;
        }
        
        /* Find the best price for the target */
        const double max_val = opt->maximise(hedge, *target, s);
        
        /* Clean up */
        delete opt;
        delete target;
        std::cout << "The maximum value is " << max_val << std::endl;
    }

    /* Clean up */
    for (unsigned int i = 0; i < hedge.size(); i++)
    {
        delete hedge[i];
    }
    
    return 0;
}
