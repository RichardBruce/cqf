#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <assert.h>

#include "matrix.h"

#include "statistics.h"


/* STL typedefs */
typedef std::vector<double> rawdata_container;


char* get_arg_string(const std::string &option, char *const argv[], const int argc, const int i)
{
    if (i == argc)
    {
        /* TODO add help message */
        //help_missing_value(option);
        return NULL;
    }
    else
    {
        return argv[i];
    }
}


int parse_line(rawdata_container *data, const std::string &line, bool ignore_border)
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


rawdata_container* parse_csv(const std::string &file_name, int *x_sz, int *y_sz)
{
    std::ifstream csv_file(file_name.c_str());
    if (!csv_file)
    {
        std::cout << "Error: Couldn't open input file: " << file_name << std::endl;
        return nullptr;
    }
    
    /* Read all the lines in the file */
    rawdata_container *data = new rawdata_container();
    
    std::string line;
    getline(csv_file, line);
    int x       = parse_line(data, line, false);
    int x_old   = x;
    int y       = 1;
    while (!csv_file.eof())
    {
        std::string line;
        getline(csv_file, line);
        ++y;
        
        /* Split the line by commas and parse the data */
        x = parse_line(data, line, false);
        
        /* Check enough data was read */
        if (x != x_old)
        {
            std::cout << "Error: Inconsistant data width" << std::endl;
            return nullptr;
        }
    }
    
    /* Clean up */
    csv_file.close();
    
    (*x_sz) = x;
    (*y_sz) = y;
    return data;
}


/* data : The data to test for stationarity                         */
/* debug: Debug file to dump debug info to                          */
/* c    : A constant for the regression                             */
/* dt   : A trend for the regression                                */
/* order: The number of delay points to include in the regression   */
template<class T>
bool dicky_fuller_test(rawdata_container &data, const std::string &debug, 
    const T c = 0.0, const T dt = 0.0, const int order = 0)
{
    /* Deltas for delay taps */
    rawdata_container deltas(data.size());
    std::adjacent_difference(data.begin(), data.end(), deltas.begin());
    
    /* Pick the number of delay points */
    /* TODO - This needs to be more intelligent */
    const int poly_order = order;
    std::cout << "ADF using order: " << order << std::endl;
    
    const int reg_size = data.size() - poly_order - 1;
    
    
    /* Perform linear regression to find a best fit curve */
    /* Build a matrix of: */
    /* [const, dt * (t - 3), Y(t - 3), dY(t - 3)...Y(t - 3 - poly_order)] */
    /* [const, dt * (t - 2), Y(t - 2), dY(t - 2)...Y(t - 2 - poly_order)] */
    /* [const, dt * (t - 1), Y(t - 1), dY(t - 1)...Y(t - 1 - poly_order)] */
    std::auto_ptr<T> reg_data(new T [reg_size * (poly_order + 3)]);
    T time          = 0.0;
    const T t_inc   = 1.0;
    for (int i = 0; i < reg_size; i++)
    {
        const int time_adr = (i * (poly_order + 3));
        const int data_adr = poly_order + 1;
        reg_data.get()[time_adr    ]  = c;
        reg_data.get()[time_adr + 1]  = dt * time;
        reg_data.get()[time_adr + 2]  = data[data_adr + i];
        for (int j = 0; j < poly_order; j++)
        {
            reg_data.get()[time_adr + j] = deltas[data_adr + i - j];
        }
        
        time += t_inc;
    }
    
    
    /* Build the pseudo inverse of the data */
    std::auto_ptr<T> svd_w(new T [poly_order + 3]);
    matrix<T> reg_values(reg_data.get(), reg_size, poly_order + 3);
    matrix<T> pseudo_inv = reg_values.moore_penrose_pseudo_inverse(svd_w.get());
    
    /* Multiply inverse by raw data to get the fitting coefficients */
    matrix<T> reg_coeff = matrix<T>(deltas.data(), 1, reg_size) * pseudo_inv;


    /* Build the fitted data */       
    std::auto_ptr<T> fitted_data(new T [reg_size * (poly_order + 4)]);
    for (int i = 0; i < reg_size; i++)
    {
        T fitted = 0.0;
        for (int j = 0; j < poly_order + 3; j++)
        {
            /* One per coefficient */
            fitted_data.get()[(j * reg_size) + i] = reg_data.get()[(i * (poly_order + 3)) + j] * reg_coeff.get_data(0, j);
            fitted += fitted_data.get()[(j * reg_size) + i];
        }
        
        /* Top level data */
        fitted_data.get()[((poly_order + 3) * reg_size) + i] = fitted;
    }
    

    /* Debug */
    if (!debug.empty())
    {
        std::ofstream debug_file(debug.c_str());
        assert(debug_file || !"ERROR: Failed to open debug file.");
        
        debug_file << "Regression Values:" << std::endl;
        reg_values.dump(debug_file);
        debug_file << std::endl; 
        
        debug_file << "Regressed Coefficients:" << std::endl;
        reg_coeff.dump(debug_file);
        debug_file << std::endl;
        
        debug_file << "Fitted Data:" << std::endl;
        for (int i = 0; i < reg_size; i++)
        {
            debug_file << deltas[i] << "," << fitted_data.get()[((poly_order + 3) * reg_size) + i] << std::endl;
        }

        debug_file.close();
    }
    

    /* Calculate t-stats of real and fitted data */
    T tstat = calculate_t_statistic(&fitted_data.get()[0], deltas.data(), 0.0, reg_coeff.get_data(0, 0), reg_size);
    std::cout << "Info: Constant has a T statistic of " << tstat << std::endl;

    tstat = calculate_t_statistic(&fitted_data.get()[reg_size], deltas.data(), 0.0, reg_coeff.get_data(0, 1), reg_size);
    std::cout << "Info: Trend has a T statistic of " << tstat << std::endl;

    tstat = calculate_t_statistic(&fitted_data.get()[(reg_size << 1)], deltas.data(), 0.0, reg_coeff.get_data(0, 2), reg_size);
    std::cout << "Info: Y(t - 1) has a T statistic of " << tstat << std::endl;
    for (int i = 0; i < poly_order; i++)
    {
        tstat = calculate_t_statistic(&fitted_data.get()[i * reg_size], deltas.data(), 0.0, reg_coeff.get_data(0, i + 3), reg_size);
        std::cout << "Info: Order " << i << " has a T statistic of " << tstat << std::endl;
    }
    
    const T se = calculate_standard_error(&fitted_data.get()[reg_size  * (poly_order + 3)], deltas.data(), reg_size);
    std::cout << "Info: Overall has a standard error of " << se << std::endl;
    std::cout << std::endl;

    
    return false;
}


int main(int argc, char *argv[])
{
    std::string historic_data;
    
    /* Parse inputs */
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-historic_data") == 0)
        {
            historic_data = get_arg_string("-historic_data", argv, argc, ++i);
        }
        else
        {
            assert(!"Unknown arguement: " + argv[i]);
        }
    }
    
    /* Read data from file */
    int raw_data_width, raw_data_length;
    rawdata_container *raw_data = parse_csv(historic_data, &raw_data_width, &raw_data_length);
    assert(raw_data != nullptr);
 
    /* Auto correlation */
    for (int i = 0; i < 10; i++)
    {
        const double ac = calculate_auto_correlation(raw_data->data(), raw_data->size(), i);
        std::cout << i << " " << ac << std::endl;
    }
    
    /* Check for stationarity */
    std::string df_dump("df_dump");
    dicky_fuller_test<double>(*raw_data, df_dump/* const, dt, order */);

    return 0;
}
