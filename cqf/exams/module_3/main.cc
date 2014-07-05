#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "standard_normal.h"
#include "implied_volatility.h"


int main(int argc, char *argv[])
{
    /* Read arguments */
    if (argc != 4)
    {
        std::cout << "Usage: " << std::endl;
        std::cout << "./module_3 <n start> <n end> <n step>" << std::endl;
        return 1;
    }
    
    const int start_n   = atoi(argv[1]);
    const int end_n     = atoi(argv[2]);
    const int n_inc     = atoi(argv[3]);
    if (fmod((double)(end_n - start_n), (double)n_inc) != 0.0)
    {
        std::cout << "Number of increments must be an integer multiple of range." << std::endl;
    }

    std::ofstream error_rates("error_rates.csv");
    error_rates << "Samples,Jn Estimate,Jn Error Rate,Kn Estimate,Kn Error Rate" << std::endl;
    for (int n = start_n; n < end_n; n += n_inc)
    {    
        const double n_dbl = (double)n;
        std::cout << "Samples: " << n << std::endl;
    
        /* Iterate */
        double jn = 0.0;
        double kn = 0.0;
        for (int i = 0; i < n; i++)
        {
            const double std_norm = standard_normal_cdf(i);
            jn += (std_norm * std_norm);
            kn += (std_norm * std_norm * 3.0) / (std_norm * std_norm * std_norm * std_norm);
        }
        jn /= n_dbl;
        kn /= n_dbl;
    
        /* Output estimate */
        //std::cout << "Jn: " << jn << std:: endl;
        std::cout << "Jn Error: " << 1.0 - jn << std::endl;
        //std::cout << "Kn Error: " << 3.0 - kn << std::endl; /* Uncomment to see kurtosis estimate */
        error_rates << n << "," << jn << "," << (1.0 - jn) << "," << kn << "," << (3.0 - kn) << std::endl;
    }
    error_rates.close();
        
    /* Back out implied vol */
    implied_volatility *imp_vol = new implied_volatility(1e-5, 10);
    double vol = imp_vol->back_out_from_call(101.5, 100.0, 0.08, 4.0 / 12.0, 6.51);
    std::cout << "Implied Volatility: " << vol * 100.0 << "%" << std::endl;

    return 0;
}
