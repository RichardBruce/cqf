#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <xmmintrin.h>

#include <assert.h>

#include "blocked_range.h"
#include "parallel_reduce.h"
#include "hjm_simulator.h"

#include "caplet_payoff.h"
#include "floorlet_payoff.h"
#include "ko_caplet_payoff.h"
#include "ko_floorlet_payoff.h"
#include "zcb_payoff.h"
#include "option_payoff.h"

#include "matrix.h"

#include "utility.h"

#include "statistics.h"

#include "dumping.h"

#include "mersenne_twister.h"
#include "box_muller.h"
#include "sobol_numbers.h"
#include "standard_normal.h"
#include "brownian_bridge.h"


using std::vector;


int parse_line(vector<double> *data, const std::string &line, bool ignore_border)
{
    int x = 1;
    size_t comma = line.find(',');
    data->push_back(atof(&line.c_str()[ignore_border ? comma + 1 : 0]));
    if (ignore_border)
    {
        comma = line.find(',', comma + 1);
    }
    
    while((comma != std::string::npos) &&(comma != line.length() - 1))
    {
        ++x;
        data->push_back(atof(&line.c_str()[comma + 1]));
        comma = line.find(',', comma + 1);
    }
    
    return x;
}


vector<double>* parse_csv(const std::string &file_name, vector<double> *maturities, int *x_sz, int *y_sz, bool ignore_border = true)
{
    std::ifstream csv_file(file_name.c_str());
    if (!csv_file)
    {
        std::cout << "Error: Couldn't open input file: " << file_name << std::endl;
        return NULL;
    }
    
    /* Read the maturities (top line of the file) */
    std::string line;
    getline(csv_file, line);
    int x = parse_line(maturities, line, ignore_border);
    int x_old = x;
    
    /* Read all the lines in the file */
    int y = 0;
    vector<double> *data = new vector<double>();
    while (!csv_file.eof())
    {
        std::string line;
        getline(csv_file, line);
        if (line.length() < 2)
        {
            break;
        }
        ++y;
        
        /* Split the line by commas and parse the data */
        x = parse_line(data, line, ignore_border);
        
        /* Check enough data was read */
        if (x != x_old)
        {
            std::cout << "Error: Inconsistant data width" << std::endl;
            return NULL;
        }
    }
    
    /* Clean up */
    csv_file.close();
    
    (*x_sz) = x;
    (*y_sz) = y;
    return data;
}


/* Calculate the volatility at rate x given the polynomial coefficients and 
   order of the polynomial. */
double evaluate_vol(const double *const vol_coeffs, const double order, const double x)
{
    /* Calculate the fitted y value at this x */
    double rolling_pow = x;
    double fitted = vol_coeffs[0];
    for (int i = 1; i < order; i++)
    {
        fitted += rolling_pow * vol_coeffs[i];
        rolling_pow *= x;
    }
    
    return fitted;
}


/* Use linear regression to fit a curve to the raw volatilities */
double* fit_vol_matrix(matrix<double> **vol, const vector<double> &maturities, const std::string dump, const int factors, const int max_order, const int min_order)
{
    /* Debug */
    std::ofstream debug_file;
    if (!dump.empty())
    {
        debug_file.open(dump.c_str());
        assert(debug_file.is_open());
    }
    
    /* Perform linear regression to find a best fit curve */
    /* Assumes the number of factors is less than the number of maturities */
    double *svd_w       = new double [max_order];
    double *data        = new double [maturities.size() * max_order];
    double *fitted_data = new double [maturities.size() * (max_order + 1)];
    double *vol_coeffs  = new double [factors * max_order];
    
    /* For each sequence of raw vol data */
    for (int i = 0; i < factors; i++)
    {
        const int polynomial_order = std::max(min_order, std::min(max_order, i + 1));

        /* Build a matrix of [const, T^1, T^2, T^3, ..., T^i] */
        for (unsigned int j = 0; j < maturities.size(); j++)
        {
            data[j * polynomial_order] = 1.0;
            double mat_pow = maturities[j]; /* Keep a rolling copy of maturities[j]^k to save calling std::pow */
            for (int k = 1; k < polynomial_order; k++)
            {
                data[(j * polynomial_order) + k] = mat_pow;
                mat_pow *= maturities[j];
            }
        }
        
        /* Build the pseudo inverse of the data */
        matrix<double> reg_values(data, maturities.size(), polynomial_order);
        matrix<double> pseudo_inv = reg_values.moore_penrose_pseudo_inverse(&svd_w[0]);
        
        /* Multiply inverse by raw data to get the fitting coefficients */
        matrix<double> reg_coeff = (*(vol[i])) * pseudo_inv;

        /* Build the fitted data */       
        /* One row for each factor and one row for the combine factors */
        for (unsigned int j = 0; j < maturities.size(); j++)
        {
            double fitted = 0.0;
            for (int k = 0; k < polynomial_order; k++)
            {
                fitted_data[(k * maturities.size()) + j] = data[(j * polynomial_order) + k] * reg_coeff.get_data(0, k);
                fitted += fitted_data[(k * maturities.size()) + j];
            }
            fitted_data[(polynomial_order * maturities.size()) + j] = fitted;
        }

        /* Debug */
        if (!dump.empty())
        {
            for (unsigned int j = 0; j < maturities.size(); j++)
            {
                debug_file << maturities[j] << "," << vol[i]->get_data(0, j) << "," << fitted_data[(polynomial_order * maturities.size()) + j] << std::endl;
            }
        }

        /* Calculate the correlation and t-stats of real and fitted data */
        if (polynomial_order > 1)
        {
            const double *const actual_data = vol[i]->extract_row(0);
            const double corr = calculate_correlation<double>(&fitted_data[polynomial_order * maturities.size()], actual_data, maturities.size());
            std::cout << "Info: Componant " << i << " has a fitted correlation of " << corr << std::endl;
            
            for (int j = 0; j < polynomial_order; j++)
            {
                const double tstat = calculate_t_statistic<double>(&fitted_data[j * maturities.size()], actual_data, 0.0, 1.0, maturities.size());
                std::cout << "Info: Order " << j << " has a T statistic of " << (reg_coeff.get_data(0, j) / tstat) << std::endl;
            }
            std::cout << std::endl;
            delete [] actual_data;
        }
        else
        {
            std::cout << "Info: Componant " << i << " is an average fitting" << std::endl << std::endl;;
        }
        
        /* Copy out the result */
        for (int j = 0; j < polynomial_order; j++)
        {
            vol_coeffs[(i * max_order) + j] = reg_coeff.get_data(0, j);
        }
        
        /* 0 Unused co efficients */
        for (int j = polynomial_order; j < max_order; j++)
        {
            vol_coeffs[(i * max_order) + j] = 0.0;
        }
    }
    
    /* Clean up */
    delete [] svd_w;
    delete [] data;
    delete [] fitted_data;
        
    debug_file.close();
    
    return vol_coeffs;
}


/* Calculate the volatilities for Monte Carlo simulation using the fitted curve */
double* build_fitted_vol_matrix(const double *const vol_coeffs, const vector<double> &rates, const int factors, const int order)
{
    /* For each required rate and factor fit a vol */
    double *vols = new double [rates.size() * factors];
    for (unsigned int i = 0; i < rates.size(); i++)
    {
        for (int j = 0; j < factors; j++)
        {
            vols[(j * rates.size()) + i] = evaluate_vol(&vol_coeffs[j * order], order, rates[i]);
        }
    }
    
    return vols;
}


/* Evaluate the HJM equation to calculate drifts from volatilities */
double* build_fitted_drift_matrix(const double *const vol_coeffs, const vector<double> &rates, const int factors, const int order)
{
    /* Build the polynomials */
    const double max_dt = 0.01;

    /* For each required rate interpolate a drift */
    double *drifts = new double[rates.size()];
    for (unsigned int i = 0; i < rates.size(); i++)
    {
        /* At time 0 the drift is 0 */
        drifts[i] = 0.0;
        if (rates[i] == 0.0)
        {
            continue;
        }
        
        /* Calculate the dt for this rate */
        /* Uses numerical integration of the vol up to point i for each factor */
        int strips = ceil(rates[i] / max_dt);
        strips += strips & 0x1;                 /* Must be even for Simpsons rule */
        const double dt = rates[i] / static_cast<double>(strips);
        for (int j = 0; j < factors; j++)
        {
            double integ = 0.0;
            for (int k = 0; k < strips; k += 2)
            {
                const double flt_j = static_cast<double>(j);
                const double s0 = evaluate_vol(&vol_coeffs[j * order], order, dt * flt_j);
                const double s1 = evaluate_vol(&vol_coeffs[j * order], order, dt * (flt_j + 1.0));
                const double s2 = evaluate_vol(&vol_coeffs[j * order], order, dt * (flt_j + 2.0));
                integ += (s0 + (4.0 * s1) + s2) * 2.0 * dt * (1.0 / 6.0);
            }
            
            const double val = evaluate_vol(&vol_coeffs[j * order], order, rates[i]);
            drifts[i] += integ * val;
        }
    }
    
    return drifts;
}


/* Use interpolation or extrapolation from the calibration data to calculate
   the initial forward curve for Monte Carlo simulation */
double* build_initial_forward_curve(const vector<double> &rates, const vector<double> &raw_rates, const double *const data, const int nr_dates)
{
    /* Must have at least 2 rates */
    assert(raw_rates.size() > 1);
    
    /* For each required rate interpolate the initial point */
    /* Assumes the raw rates are sorted */
    const double percent_to_factor = 1.0 / 100.0;
    unsigned int raw_index = 0;
    double *fwd_curve = new double [rates.size() * nr_dates];
    for (unsigned int i = 0; i < rates.size(); i++)
    {
        /* Iterate to or pasted the required rate */
        while ((raw_index < raw_rates.size()) && (raw_rates[raw_index] < rates[i]))
        {
            ++raw_index;
        }
        
        /* The raw data doesnt extend far enough for this required rate */
        assert(raw_index < raw_rates.size());
        
        /* Exact match, copy it over */
        if (rates[i] == raw_rates[raw_index])
        {
            fwd_curve[i] = data[raw_index] * percent_to_factor;
        }
        /* Extrapolate the spot rate */
        else if (raw_index == 0)
        {
            const double d = (data[1] - data[0]) / (raw_rates[1] - raw_rates[0]);
            fwd_curve[i] = (data[0] - (d * raw_rates[0])) * percent_to_factor;
        }
        /* Interpolate between raw_index - 1 and raw_index */
        else
        {
            const double d = raw_rates[raw_index] - raw_rates[raw_index - 1];
            const double w0 = (raw_rates[raw_index] - rates[i]) / d;
            const double w1 = (rates[i] - raw_rates[raw_index - 1]) / d;
            fwd_curve[i] = ((data[raw_index - 1] * w0) + (data[raw_index] * w1)) * percent_to_factor;
        }
    }
    
    return fwd_curve;
}


#ifdef USE_SIMD
void build_monte_carlo_points_simd(double *rands, double *fwd_curve, const vector<double> &dr_inv, const vector<double> &dt, const vector<double> &sqrt_dt, double *vol, const double *drifts, const int factors)
{
    /* Number of rates must be a multiple of SIMD width */
    assert((dr_inv.size() & 0x1) == 0);
    
    unsigned int t_idx = dr_inv.size();
    unsigned int last_t_idx = 0;
    for (unsigned int i = 1; i < dt.size(); i++)
    {
        /* Build the forward curve a point at a time */
        for (unsigned int j = 0; j < dr_inv.size() - SIMD_WIDTH; j += SIMD_WIDTH)
        {
            const __m128d last_t = _mm_loadu_pd(&fwd_curve[last_t_idx]);
            last_t_idx += SIMD_WIDTH;
            
            const __m128d drift = _mm_loadu_pd(&drifts[j]);
            
            unsigned int rand_idx = i - 1;
            unsigned int vol_idx = j;
            __m128d diff = _mm_mul_pd(_mm_load_pd1(&rands[rand_idx]), _mm_loadu_pd(&vol[vol_idx]));
            for (int k = 1; k < factors; k++)
            {
                rand_idx += (dt.size() - 1);
                vol_idx += dr_inv.size();
                __m128d mul = _mm_mul_pd(_mm_load_pd1(&rands[rand_idx]), _mm_loadu_pd(&vol[vol_idx]));
                diff = _mm_add_pd(diff, mul);
            }
            diff = _mm_mul_pd(diff, _mm_load_pd1(&sqrt_dt[i]));
            
            const __m128d musiela = _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(&fwd_curve[last_t_idx - 1]), last_t), _mm_loadu_pd(&dr_inv[j + 1]));
            const __m128d res = _mm_add_pd(_mm_add_pd(last_t, diff), _mm_mul_pd(_mm_add_pd(drift, musiela), _mm_load1_pd(&dt[i])));
            _mm_store_pd(&fwd_curve[t_idx], res);
            t_idx += SIMD_WIDTH;
        }
        
        /* Last point is built slightly differently */
        const __m128d last_t = _mm_loadu_pd(&fwd_curve[last_t_idx]);
        last_t_idx += SIMD_WIDTH;

        const __m128d drift = _mm_loadu_pd(&drifts[dr_inv.size() - SIMD_WIDTH]);
        
        unsigned int rand_idx = i - 1;
        unsigned int vol_idx = dr_inv.size() - SIMD_WIDTH;
        __m128d diff = _mm_mul_pd(_mm_load_pd1(&rands[rand_idx]), _mm_loadu_pd(&vol[vol_idx]));
        for (int k = 1; k < factors; k++)
        {
            rand_idx += (dt.size() - 1);
            vol_idx += dr_inv.size();
            __m128d mul = _mm_mul_pd(_mm_load_pd1(&rands[rand_idx]), _mm_loadu_pd(&vol[vol_idx]));
            diff = _mm_add_pd(diff, mul);
        }
        diff = _mm_mul_pd(diff, _mm_load_pd1(&sqrt_dt[i]));
        
        const __m128d musiela = _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(&fwd_curve[last_t_idx - 1]), last_t), _mm_loadu_pd(&dr_inv[dr_inv.size() - 2]));
        const __m128d res = _mm_add_pd(_mm_add_pd(last_t, diff), _mm_mul_pd(_mm_add_pd(drift, musiela), _mm_load1_pd(&dt[i])));
        _mm_store_pd(&fwd_curve[t_idx], res);
        t_idx += SIMD_WIDTH;
    }
}
#endif


void build_monte_carlo_points(double *rands, double *fwd_curve, const vector<double> &dr_inv, const vector<double> &dt, const vector<double> &sqrt_dt, double *vol, const double *drifts, const int factors)
{
    unsigned int t_idx = dr_inv.size();
    unsigned int last_t_idx = 0;
    for (unsigned int i = 1; i < dt.size(); i++)
    {
        /* Build the forward curve a point at a time */
        double musiela = 0.0;
        for (unsigned int j = 0; j < dr_inv.size() - 1; j++)
        {
            const double last_t = fwd_curve[last_t_idx++];
            
            double diff = rands[i - 1] * vol[j];
            for (int k = 1; k < factors; k++)
            {
                diff += rands[i - 1 + (k * (dt.size() - 1))] * vol[(k * dr_inv.size()) + j];
            }
            diff *= sqrt_dt[i];
            
            musiela = ((fwd_curve[last_t_idx] - last_t) * dr_inv[j + 1]);
            fwd_curve[t_idx++] = last_t + diff + (drifts[j] + musiela) * dt[i];
        }
        
        /* Last point is built slightly differently */
        const double last_t = fwd_curve[last_t_idx++];
        
        double diff = rands[i - 1] * vol[dr_inv.size() - 1];
        for (int j = 1; j < factors; j++)
        {
            diff += rands[i - 1 + (j * (dt.size() - 1))] * vol[(j * dr_inv.size()) + dr_inv.size() - 1];
        }
        diff *= sqrt_dt[i];
        
        fwd_curve[t_idx++] = last_t + diff + (drifts[dr_inv.size() - 1] + musiela) * dt[i];
    }
}


void build_monte_carlo_points(double *rands, double *fwd_curve, const vector<double> &dr_inv, brownian_bridge<double> **bridges, const vector<double> &dt, const vector<double> &sqrt_dt, double *vol, const double *drifts, const int factors)
{
    double *mapped_rands = new double [factors * dr_inv.size() * dt.size()];
    for (unsigned int i = 0; i < dr_inv.size(); i++)
    {
        bridges[i                      ]->map_randoms(&rands[0                  ], &mapped_rands[(dt.size() - 1) * i]);
        bridges[i +      dr_inv.size() ]->map_randoms(&rands[(dt.size() - 1)    ], &mapped_rands[(dt.size() - 1) * (i + 1)]);
        bridges[i + (2 * dr_inv.size())]->map_randoms(&rands[(dt.size() - 1) * 2], &mapped_rands[(dt.size() - 1) * (i + 2)]);
    }    
    
    unsigned int t_idx = dr_inv.size();
    unsigned int last_t_idx = 0;
    for (unsigned int i = 1; i < dt.size(); i++)
    {
        /* Build the forward curve a point at a time */
        double musiela = 0.0;
        for (unsigned int j = 0; j < dr_inv.size() - 1; j++)
        {
            const double last_t = fwd_curve[last_t_idx++];
            
            double diff = rands[i - 1] * vol[j];
            for (int k = 1; k < factors; k++)
            {
                diff += mapped_rands[(i - 1) + ((dt.size() - 1) * ((j * factors) + k))];//rands[i - 1 + (k * (dt.size() - 1))] * vol[(k * dr_inv.size()) + j];
            }
            diff *= sqrt_dt[i];
            
            musiela = ((fwd_curve[last_t_idx] - last_t) * dr_inv[j + 1]);
            fwd_curve[t_idx++] = last_t + diff + (drifts[j] + musiela) * dt[i];
        }
        
        /* Last point is built slightly differently */
        const double last_t = fwd_curve[last_t_idx++];
        
        double diff = rands[i - 1] * vol[dr_inv.size() - 1];
        for (int j = 1; j < factors; j++)
        {
            diff += mapped_rands[(i - 1) + ((dt.size() - 1) * (((dr_inv.size() - 1) * factors) + j))];//rands[i - 1 + (j * (dt.size() - 1))] * vol[(j * dr_inv.size()) + dr_inv.size() - 1];
        }
        diff *= sqrt_dt[i];
        
        fwd_curve[t_idx++] = last_t + diff + (drifts[dr_inv.size() - 1] + musiela) * dt[i];
    }
    
    delete [] mapped_rands;
}


void usage()
{
    std::cout << "Usage:" << std::endl;
    std::cout << "hjm -historic_data <str> -zcb 1 10 1000000"                                       << std::endl;
    std::cout << "Calibration Control"                                                              << std::endl;
    std::cout << "  -historic_data f            Input csv file of historic rate data from f."       << std::endl;
    std::cout << "  -factors f                  Use f factors for PCA, range 0-9."                  << std::endl;
    std::cout << "  -pca_explain p              Explain p percent of the curve by PCA, range 0-1."  << std::endl;
    std::cout << "  -min_curve_order p          Minimum order of vol curve to fit, range 0-20."     << std::endl;
    std::cout << "  -max_curve_order p          Maximum order of vol curve to fit, range 0-20."     << std::endl;
    std::cout << "  -calibrate_only             Only perform PCA and line fitting. Dont run MC."    << std::endl;
    std::cout << std::endl;
    std::cout << "If -factors and -pca_explain are specified then both criteria"                    << std::endl;
    std::cout << "will be met. -factors and -pca_explain are both minimum criteria."                << std::endl;
    std::cout << std::endl;
    std::cout << "MC Simulation Control"                                                            << std::endl;
    std::cout << "  -nr_of_paths d              The number of path to use in MC simulation."        << std::endl;
    std::cout << "  -max_dt d                   Longest time step to allow in MC simulation."       << std::endl;
    std::cout << "  -max_dr d                   Longest rate step to allow in MC simulation."       << std::endl;
    std::cout << "  -min_long_rate r            Simulate rates up to atleast r."                    << std::endl;
    std::cout << "  -sobol f                    Use Sobol numbers with the direction numbers in f." << std::endl;
    std::cout << "  -nbrownian_bridge           Dont use Brownian Bridge path construction."        << std::endl;
    std::cout << "  -threaded                   Use threaded Monte Carlo simulation."               << std::endl;
    std::cout << std::endl;
    std::cout << "Brownian bridge path contruction will by default be used with sobol numbers."     << std::endl;
    std::cout << std::endl;
    std::cout << "Payoffs"                                                                          << std::endl;
    std::cout << "  -compound_freq d            Set the compounding frequency."                     << std::endl;
    std::cout << "  -option d e p s c           Price an option on the following payoffs."          << std::endl;
    std::cout << "  -option_end                 Stop adding payoffs to the current option."         << std::endl;
    std::cout << "  -caplet d e p t n s r       Price a caplet."                                    << std::endl;
    std::cout << "  -floorlet d e p t n s r     Price a floorlet."                                  << std::endl;
    std::cout << "  -ko_caplet d e p t n s r b  Price a knock out caplet."                          << std::endl;
    std::cout << "  -ko_floorlet d e p t n s r bPrice a knock out floorlet."                        << std::endl;
    std::cout << "  -zcb s t n                  Price a ZCB."                                       << std::endl;
    std::cout << "                              notional n."                                        << std::endl;
    std::cout << std::endl;
    std::cout << "For all payoffs: "                                                                << std::endl;
    std::cout << "  d is the start date."                                                           << std::endl;
    std::cout << "  e is the evaluation date."                                                      << std::endl;
    std::cout << "  p is the payment date."                                                         << std::endl;
    std::cout << "  t is the tenor."                                                                << std::endl;
    std::cout << "  n is the notional."                                                             << std::endl;
    std::cout << "  s is the strike."                                                               << std::endl;
    std::cout << "  r is the underlying rate."                                                      << std::endl;
    std::cout << "  b is the barrier."                                                              << std::endl;
    std::cout << "  c is 0 for a put option else a call option."                                    << std::endl;
    std::cout << "Atleast one payoff must be specified."                                            << std::endl;
    std::cout << "All times are in years so you can enjoy playing with the day count"               << std::endl;
    std::cout << "conventions you like."                                                            << std::endl;
    std::cout << "The compounding frequency should be set before the payoffs it affect."            << std::endl;
    std::cout << "The default compounding frequency is 12 times a year."                            << std::endl;
    std::cout << "All payoffs requested between -option and -option_end are added to that option."  << std::endl;
    std::cout << "Multiple options may be defined in one pricing."                                  << std::endl;
    std::cout << std::endl;
    std::cout << "Debug"                                                                            << std::endl;
    std::cout << "-dump_sim_rates f             Dump the rates that will be simulated to f."        << std::endl;
    std::cout << "-dump_sim_dates f             Dump the dates that will be simulated to f."        << std::endl;
    std::cout << "-dump_raw_data f              Dump the raw data used for calibration to f."       << std::endl;
    std::cout << "-dump_dif_data f              Dump the differenced data used for calibration."    << std::endl;
    std::cout << "-dump_covar_data f            Dump the co variance data used for calibration."    << std::endl;
    std::cout << "-dump_pca_size f              Dump the eigen values from pca to f."               << std::endl;
    std::cout << "-dump_pca_vecs f              Dump the eigen vectors from pca to f."              << std::endl;
    std::cout << "-dump_pca_main f              Dump the most significant eigen vectors to f."      << std::endl;
    std::cout << "-dump_fitted_vol f            Dump the volatiliy curves after fitting to f."      << std::endl;
    std::cout << "-dump_vol_coeffs f            Dump the coefficients from linear regression to f." << std::endl;
    std::cout << "-dump_sim_vols f              Dump the volatiliy curves after fitting to f."      << std::endl;
    std::cout << "-dump_sim_drifts f            Dump the calculated drifts to f."                   << std::endl;
    std::cout << "-dump_initial_data f          Dump the initial rates used for simulation to f."   << std::endl;
    std::cout << "-dump_conv_data f             Dump convergence graph to f."                       << std::endl;
}


void help_missing_value(const std::string &option)
{
    std::cout << "Error: No value for option " << option << std::endl;
    usage();
    assert(false);
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

void add_payoff(std::vector<payoff*> *alone, std::vector<payoff*> *option, payoff *p, const double option_e, const double payoff_e)
{
    if (option == NULL)
    {
        alone->push_back(p);
    }
    /* Add the payoff to an option */
    else
    {
        /* Check the payoff begins on the option evaluation date */
        assert(payoff_e == option_e);
        option->push_back(p);
    }
}

int main(int argc, char *argv[])
{
    /* Default parameters */
    bool calibrate_only         = false;
    bool bridge_en              = true;
    bool threaded_mc            = false;
    int factors                 = 3;
    int min_curve_order         = 0;
    int max_curve_order         = 9;
    int nr_of_paths             = 100000;
    double pca_explain          = 0.0;
    double max_dt               = 0.1;
    double max_dr               = 1.0;
    double min_long_rate        = 0.0;
    double comp_freq            = 12.0;
    std::string historic_data("bond_data.csv");
    std::string direction_numbers;
    
    /* Payoff construction */
    double option_d             = 0.0;
    double option_e             = 0.0;
    double option_p             = 0.0;
    double option_s             = 0.0;
    bool option_c               = 0.0;
    vector<payoff*> *option_on  = NULL;
    vector<payoff *> payoffs;
    
    /* File names to dump to */
    std::string dump_sim_rates;
    std::string dump_sim_dates;
    std::string dump_raw_data;
    std::string dump_dif_data;
    std::string dump_covar_data;
    std::string dump_pca_size;
    std::string dump_pca_vecs;
    std::string dump_pca_main;
    std::string dump_fitted_vol;
    std::string dump_vol_coeffs;
    std::string dump_sim_vols;
    std::string dump_sim_drifts;
    std::string dump_initial_data;
    std::string dump_conv_data;
    
    /* Parse inputs */
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-historic_data") == 0)
        {
            historic_data = get_arg_string("-historic_data", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-help") == 0)
        {
            usage();
        }
        else if (strcmp(argv[i], "-caplet") == 0)
        {
            if ((i + 7) >= argc)
            {
                help_missing_value("-caplet");
            }
            else
            {
                /* Read caplet data */
                const double d = atof(argv[i + 1]);
                const double e = atof(argv[i + 2]);
                const double p = atof(argv[i + 3]);
                const double t = atof(argv[i + 4]);
                const double n = atof(argv[i + 5]);
                const double s = atof(argv[i + 6]);
                const double r = atof(argv[i + 7]);
                payoff* caplet = new caplet_payoff(d, e, p, t, n, s, r, comp_freq);
                add_payoff(&payoffs, option_on, caplet, option_e, d);
                i += 7;
            }
        }
        else if (strcmp(argv[i], "-floorlet") == 0)
        {
            if ((i + 7) >= argc)
            {
                help_missing_value("-floorlet");
            }
            else
            {
                /* Read floorlet data */
                const double d = atof(argv[i + 1]);
                const double e = atof(argv[i + 2]);
                const double p = atof(argv[i + 3]);
                const double t = atof(argv[i + 4]);
                const double n = atof(argv[i + 5]);
                const double s = atof(argv[i + 6]);
                const double r = atof(argv[i + 7]);
                payoff *floorlet = new floorlet_payoff(d, e, p, t, n, s, r, comp_freq);
                add_payoff(&payoffs, option_on, floorlet, option_e, d);
                i += 7;
            }
        }
        else if (strcmp(argv[i], "-ko_caplet") == 0)
        {
            if ((i + 8) >= argc)
            {
                help_missing_value("-ko_caplet");
            }
            else
            {
                /* Read knock out caplet data */
                const double d = atof(argv[i + 1]);
                const double e = atof(argv[i + 2]);
                const double p = atof(argv[i + 3]);
                const double t = atof(argv[i + 4]);
                const double n = atof(argv[i + 5]);
                const double s = atof(argv[i + 6]);
                const double r = atof(argv[i + 7]);
                const double b = atof(argv[i + 8]);
                payoff* ko_caplet = new ko_caplet_payoff(d, e, p, t, n, s, r, b, comp_freq);
                add_payoff(&payoffs, option_on, ko_caplet, option_e, d);
                i += 8;
            }
        }
        else if (strcmp(argv[i], "-ko_floorlet") == 0)
        {
            if ((i + 8) >= argc)
            {
                help_missing_value("-ko_floorlet");
            }
            else
            {
                /* Read knock out floorlet data */
                const double d = atof(argv[i + 1]);
                const double e = atof(argv[i + 2]);
                const double p = atof(argv[i + 3]);
                const double t = atof(argv[i + 4]);
                const double n = atof(argv[i + 5]);
                const double s = atof(argv[i + 6]);
                const double r = atof(argv[i + 7]);
                const double b = atof(argv[i + 8]);
                payoff *ko_floorlet = new ko_floorlet_payoff(d, e, p, t, n, s, r, b, comp_freq);
                add_payoff(&payoffs, option_on, ko_floorlet, option_e, d);
                i += 8;
            }
        }
        else if (strcmp(argv[i], "-zcb") == 0)
        {
            if ((i + 3) >= argc)
            {
                help_missing_value("-zcb");
            }
            else
            {
                /* Read ZCB data */
                const double d = atof(argv[i + 1]);
                const double t = atof(argv[i + 2]);
                const double n = atof(argv[i + 3]);
                
                /* Add the payoff standalone */
                payoff * zcb = new zcb_payoff(d, t, n);
                add_payoff(&payoffs, option_on, zcb, option_e, d);
                i += 3;
            }
        }
        else if (strcmp(argv[i], "-option") == 0)
        {
            if ((i + 5) >= argc)
            {
                help_missing_value("-option");
            }
            else
            {
                /* Read option data */
                option_d = atof(argv[i + 1]);
                option_e = atof(argv[i + 2]);
                option_p = atof(argv[i + 3]);
                option_s = atof(argv[i + 4]);
                option_c = (atoi(argv[i + 5]) != 0);
                option_on = new std::vector<payoff*>();
                i += 5;
            }
        }
        else if (strcmp(argv[i], "-option_end") == 0)
        {
            /* The option must be on something */
            assert((option_on != NULL) && (option_on->size() > 0));
            payoffs.push_back(new option_payoff(option_on, option_d, option_p, option_s, option_c));
            option_on = NULL;
        }
        else if (strcmp(argv[i], "-compound_freq") == 0)
        {
            comp_freq = get_arg_float("-compound_freq", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nr_of_paths") == 0)
        {
            nr_of_paths = get_arg_int("-nr_of_paths", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-max_dt") == 0)
        {
            max_dt = get_arg_float("-max_dt", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-max_dr") == 0)
        {
            max_dr = get_arg_float("-max_dr", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-min_long_rate") == 0)
        {
            min_long_rate = get_arg_float("-min_long_rate", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-sobol") == 0)
        {
            direction_numbers = get_arg_string("-sobol", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-nbrownian_bridge") == 0)
        {
            bridge_en = false;
        }
        else if (strcmp(argv[i], "-threaded") == 0)
        {
            threaded_mc = true;
        }
        else if (strcmp(argv[i], "-factors") == 0)
        {
            factors = std::max(1, get_arg_int("-factors", argv, argc, ++i));
        }
        else if (strcmp(argv[i], "-min_curve_order") == 0)
        {
            min_curve_order = std::min(20, std::max(0, get_arg_int("-min_curve_order", argv, argc, ++i))) + 1;
        }
        else if (strcmp(argv[i], "-max_curve_order") == 0)
        {
            max_curve_order = std::min(20, std::max(0, get_arg_int("-max_curve_order", argv, argc, ++i))) + 1;
        }
        else if (strcmp(argv[i], "-pca_explain") == 0)
        {
            pca_explain = std::max(0.0, std::min(1.0, get_arg_float("-pca_explain", argv, argc, ++i)));
        }
        else if (strcmp(argv[i], "-calibrate_only") == 0)
        {
            calibrate_only = true;
        }
        else if (strcmp(argv[i], "-dump_sim_rates") == 0)
        {
            dump_sim_rates = get_arg_string("-dump_sim_rates", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_sim_dates") == 0)
        {
            dump_sim_dates = get_arg_string("-dump_sim_dates", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_raw_data") == 0)
        {
            dump_raw_data = get_arg_string("-dump_raw_data", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_dif_data") == 0)
        {
            dump_dif_data = get_arg_string("-dump_dif_data", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_covar_data") == 0)
        {
            dump_covar_data = get_arg_string("-dump_covar_data", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_pca_size") == 0)
        {
            dump_pca_size = get_arg_string("-dump_pca_size", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_pca_vecs") == 0)
        {
            dump_pca_vecs = get_arg_string("-dump_pca_vecs", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_pca_main") == 0)
        {
            dump_pca_main = get_arg_string("-dump_pca_main", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_fitted_vol") == 0)
        {
            dump_fitted_vol = get_arg_string("-dump_fitted_vol", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_vol_coeffs") == 0)
        {
            dump_vol_coeffs = get_arg_string("-dump_vol_coeffs", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_sim_vols") == 0)
        {
            dump_sim_vols = get_arg_string("-dump_sim_vols", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_sim_drifts") == 0)
        {
            dump_sim_drifts = get_arg_string("-dump_sim_drifts", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_initial_data") == 0)
        {
            dump_initial_data = get_arg_string("-dump_initial_data", argv, argc, ++i);
        }
        else if (strcmp(argv[i], "-dump_conv_data") == 0)
        {
            dump_conv_data = get_arg_string("-dump_conv_data", argv, argc, ++i);
        }
        else
        {
            std::cout << "Error: Unknown input " << argv[i] << std::endl;
            usage();
            assert(false);
        }
    }
    
    /* There must be some payoffs for us to simulate */
    assert(calibrate_only || (payoffs.size() != 0));
    
    /* Max fitting curve order must be greater than min fitting curve order */
    assert(max_curve_order >= min_curve_order);
    
    /* Extract the rates and dates that need simulating */
    vector<double> required_rates;
    vector<double> required_dates;
    
    /* Always simulate the spot rate for discounting and include today in the data */
    required_rates.push_back(0.0);
    required_dates.push_back(0.0);
    std::back_insert_iterator<std::vector<double> > rates_inserter(required_rates);
    std::back_insert_iterator<std::vector<double> > dates_inserter(required_dates);
    for (unsigned int i = 0; i < payoffs.size(); i++)
    {
        payoffs[i]->required_rate(rates_inserter);
        payoffs[i]->required_date(dates_inserter);
    }
    
    /* Clean up any duplicates */
    std::sort(required_rates.begin(), required_rates.end());
    std::sort(required_dates.begin(), required_dates.end());
    vector<double>::iterator rates_iter = std::unique(required_rates.begin(), required_rates.end());
    vector<double>::iterator dates_iter = std::unique(required_dates.begin(), required_dates.end());
    required_rates.resize(rates_iter - required_rates.begin());
    required_dates.resize(dates_iter - required_dates.begin());
    if (min_long_rate > required_rates.back())
    {
        required_rates.push_back(min_long_rate);
    }
    
    /* Pad the dates and rates for the minimum step size */
    vector<double> padded_rates;
    vector<double> padded_dates;
    pad_data<double>(&padded_rates, required_rates, max_dr);
    pad_data<double>(&padded_dates, required_dates, max_dt);
    while (padded_rates.size() < 2)
    {
        padded_rates.push_back(padded_rates.back() + max_dr);
    }
#ifdef USE_SIMD
    if ((padded_rates.size() % SIMD_WIDTH) != 0)
    {
        padded_rates.push_back(padded_rates.back() + max_dr);
    }
#endif
    
    dump_1d_data_to_file(dump_sim_dates, &padded_dates[0], padded_dates.size());
    dump_1d_data_to_file(dump_sim_rates, &padded_rates[0], padded_rates.size());
 
    
    /* Parse the csv file */
    int x, y;
    vector<double> maturities;
    vector<double> *data = parse_csv(historic_data, &maturities, &x, &y);
    if (data == NULL)
    {
        std::cout << "Fatal: Can't recover from previous error(s) in parse_csv" << std::endl;
        std::cout << "Fatal: Exitting" << std::endl;
        return 1;
    }
    std::cout << "Info: Parsed " << y << " lines, with " << x << " points per line" << std::endl;;
    dump_2d_data_to_file(dump_raw_data, &((*data)[0]), x, y);

    /* Difference the time series data */
    double* diff = build_difference_matrix(*data, x, y);
    --y;
    dump_2d_data_to_file(dump_dif_data, &diff[0], x, y);
    
    /* Build the co variance matrix from the differences */
    matrix<double> covar = build_covariance_matrix(diff, x, y);
    covar.dump(dump_covar_data);

    /* Principle componant analysis */
    double *eigen_values = new double[x];
    covar.eigen_system(eigen_values, Symetric);
    
    dump_1d_data_to_file(dump_pca_size, &eigen_values[0], x);
    covar.dump(dump_pca_vecs);

    /* Pick the main componants */
    double eigen_sum = 0.0;
    vector<std::pair<double, int> > ev_indexes;
    for (int i = 0; i < x; i++)
    {
        eigen_sum += eigen_values[i];
        ev_indexes.push_back(std::make_pair(eigen_values[i], i));
    }
    
    /* The inequality operator is defined for std::pair to compare .first */
    std::sort(ev_indexes.begin(), ev_indexes.end());

    /* Build the raw volatility data */
    /* This maybe over sized, but covers the case of -pca_explain = 1 */
    matrix<double> **vol = new matrix<double>* [x];
    const double eigen_sum_inv = 1.0 / eigen_sum;
    double factor_sum = 0.0;
    int cur_factor = 0;
    while ((cur_factor < factors) || ((factor_sum * eigen_sum_inv) < pca_explain))
    {
        const std::pair<double, int> cur_pair = ev_indexes[x - cur_factor - 1];
        vol[cur_factor] = new matrix<double>(covar.extract_column(cur_pair.second) * std::sqrt(cur_pair.first));
        factor_sum += cur_pair.first;
        ++cur_factor;
    }
    factors = cur_factor;
    std::cout << "Info: Using " << factors << " factors" << std::endl;
    std::cout << "Info: Factors explain " << (factor_sum * eigen_sum_inv) * 100.0 << "% of the variance" << std::endl;
    if (!dump_pca_main.empty())
    {
        std::ofstream file(dump_pca_main.c_str());
        assert(file.is_open());
        for (int i = 0; i < factors; i++)
        {
            file << ev_indexes[x - i - 1].first << " -> ";
            vol[i]->dump(file);
        }
        file.close();
    }    
    
    /* Fit the volatility matrix with linear regression */
    double *vol_coeffs = fit_vol_matrix(vol, maturities, dump_fitted_vol, factors, max_curve_order, min_curve_order);
    dump_2d_data_to_file(dump_vol_coeffs, vol_coeffs, max_curve_order, factors);
    
    /* Build the fitted volatility matrix */
    double *fitted_vols = build_fitted_vol_matrix(vol_coeffs, padded_rates, factors, max_curve_order);
    dump_2d_data_to_file(dump_sim_vols, fitted_vols, padded_rates.size(), factors);
    
    /* Build the drift matrix */
    double *fitted_drifts = build_fitted_drift_matrix(vol_coeffs, padded_rates, factors, max_curve_order);
    dump_1d_data_to_file(dump_sim_drifts, fitted_drifts, padded_rates.size());
    
    /* Build the initial forward curve */
    double *fwd_curve = build_initial_forward_curve(padded_rates, maturities, &(*data)[x * y], padded_dates.size());
    dump_1d_data_to_file(dump_initial_data, fwd_curve, padded_rates.size());
    
    /* Check if complete */
    double average_value = 0.0;
    if (calibrate_only)
    {
        std::cout << "Calibration complete exitting" << std::endl;    
    }
    else if (threaded_mc)
    {
        /* Run the Monte Carlo simualtion */
        hjm_simulator hjm_sim(direction_numbers, payoffs, fwd_curve, padded_rates, padded_dates, fitted_vols, fitted_drifts, nr_of_paths, factors);
        tbb::parallel_reduce(tbb::blocked_range<size_t>(0, nr_of_paths), hjm_sim);
        
        average_value = hjm_sim.value();
    }
    else
    {
        /* Dumping of convergence diagram */
        std::ofstream conv_file;
        if (!dump_conv_data.empty())
        {
            conv_file.open(dump_conv_data.c_str());
            assert(conv_file);
            conv_file << "Iteration,Value" << std::endl;
        }

        const double nr_of_paths_inv = 1.0 / static_cast<double>(nr_of_paths);

        /* Calculate the time steps */
        vector<double> dt(padded_dates.size());
        std::adjacent_difference(padded_dates.begin(), padded_dates.end(), dt.begin());
        std::vector<double> sqrt_dt(dt);
        for (unsigned int i = 0; i < sqrt_dt.size(); i++)
        {
            sqrt_dt[i] = std::sqrt(sqrt_dt[i]);
        }

        /* Calculate the rate steps */
        vector<double> dr_inv(padded_rates.size());
        std::adjacent_difference(padded_rates.begin(), padded_rates.end(), dr_inv.begin());
        for (unsigned int i = 0; i < dr_inv.size(); i++)
        {
            dr_inv[i] = 1.0 / dr_inv[i];
        }


        std::string mc_points_dump = "mc_points_";
        double *gen_nums    = new double [(padded_dates.size() - 1) * factors];
        double *normal_nums = NULL;
        mersenne_twister        *mt_rand    = NULL;
        sobol_numbers<double>   *sobol_rand = NULL;

        /* Brownian Bridge to help reduce variance with sobol numbers */
        if (!direction_numbers.empty())
        {
            /* Build low descrepancy number generators */
            sobol_rand = new sobol_numbers<double>(direction_numbers.c_str(), nr_of_paths + 1, ((padded_dates.size() - 1) * factors));

            /* Skip 0 because it doesnt go well in an inverse normal */
            sobol_rand->get_next(&gen_nums[0]);
            
            /* Extra storage required for brownian bridge */
            normal_nums = new double [(padded_dates.size() - 1) * factors];
        }
        else
        {
            /* Build mersenne twister normally distributed random number generators */
            mt_rand = new mersenne_twister(5489);
        }

    /* Foreach time step */
//    brownian_bridge<double> **bridges = new brownian_bridge<double>*[factors * padded_rates.size()];
    brownian_bridge<double> *bridge = new brownian_bridge<double>(padded_dates.size() - 1);
//    for (unsigned int i = 0; i < (factors * padded_rates.size()); i++)
//    {
//        bridges[i] = new brownian_bridge<double>(padded_dates, fitted_vols[i], padded_dates.size() - 1);
//    }

        /*  Run lots of sims */
        for (int i = 0; i < nr_of_paths; i++)
        {
            if (sobol_rand != NULL)
            {
                /* Get the low descrepency numbers for all dimensions */
                sobol_rand->get_next(&gen_nums[0]);
                for (unsigned int j = 0; j < ((padded_dates.size() - 1) * factors); j++)
                {
                    normal_nums[j] = inverse_standard_normal_cdf(gen_nums[j]);
                }
                
                /* Apply the brownian bridge or not */
                if (bridge_en)
                {
                    for (int j = 0; j < factors; j++)
                    {
                        bridge->map_randoms(&normal_nums[j * (padded_dates.size() - 1)], &gen_nums[j * (padded_dates.size() - 1)]);
                    }
                }
                else
                {
                    memcpy(&gen_nums[0], &normal_nums[0], (factors * (padded_dates.size() - 1) * sizeof(double)));
                }
            }
            else
            {
                /* MT random numbers for this sim */
                for (unsigned int j = 0; j < ((padded_dates.size() - 1) * factors); j++)
                {
                    gen_nums[j] = inverse_standard_normal_cdf(mt_rand->get_next_prob());
                }
            }
            
            /* Run the sim */
#ifdef USE_SIMD
            build_monte_carlo_points_simd(&gen_nums[0], fwd_curve, dr_inv, dt, sqrt_dt, fitted_vols, fitted_drifts, factors);
#else
            build_monte_carlo_points(&gen_nums[0], fwd_curve, dr_inv, dt, sqrt_dt, fitted_vols, fitted_drifts, factors);
#endif
//            if (i == 10000)
//            {
//                std::stringstream out;
//                out << i;
//                dump_2d_data_to_file(mc_points_dump + out.str(), fwd_curve, padded_rates.size(), padded_dates.size());
//            }
                        
            /* Calculate the total payoff of this sim */
            double v = 0.0;
            for (unsigned int j = 0; j < payoffs.size(); j++)
            {
                const double payoff_comp = payoffs[j]->evaluate(padded_rates, padded_dates, fwd_curve);
                v += payoff_comp;
            }
            
            average_value += (v * nr_of_paths_inv);
            if (!dump_conv_data.empty() && ((i % 10) == 0))
            {
                conv_file << i << "," << (average_value * (nr_of_paths / static_cast<double>(i + 1))) << std::endl;
            }
        }

        /* Clean up */
        delete [] gen_nums;
        
        if (sobol_rand != NULL)
        {
            delete sobol_rand;
//            delete [] normal_nums;
        }
        else
        {
            delete mt_rand;
        }
    }
        
    /* The result */
    std::cout << "Value: " << average_value << std::endl;

    /* Clean up */
    delete data;
    delete [] diff;
    delete [] eigen_values;
    delete [] vol_coeffs;
    delete [] fitted_vols;
    delete [] fitted_drifts;
    delete [] fwd_curve;
    
    for (int i = 0; i < factors; i++)
    {
        delete vol[i];
    }
    delete [] vol;
    
    for (unsigned int i = 0; i < payoffs.size(); i++)
    {
        delete payoffs[i];
    }
    
    return 0;
}
