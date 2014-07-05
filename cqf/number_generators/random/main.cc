#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#include "mersenne_twister.h"
#include "standard_normal.h"
#include "bivariate_standard_normal.h"
#include "box_muller.h"
#include "european_call.h"


using namespace std;

int main(int argc, char *argv[])
{
    const double xs[28]     = {  0.0,  0.0,  0.0,
                                 0.0,  0.0,  0.0,
                                 0.0,  0.0,  0.0,
                                -0.5, -0.5, -0.5,
                                -0.5, -0.5, -0.5,
                                -0.5, -0.5, -0.5,
                                 0.5,  0.5,  0.5,
                                 0.5,  0.5,  0.5,
                                 0.5,  0.5,
                                 0.0,  0.00001
                                };
                                
    const double ys[28]     = {  0.0,  0.0,  0.0,
                                -0.5, -0.5, -0.5,
                                 0.5,  0.5,  0.5,
                                 0.0,  0.0,  0.0,
                                -0.5, -0.5, -0.5,
                                 0.5,  0.5,  0.5,
                                 0.0,  0.0,  0.0,
                                -0.5, -0.5, -0.5,
                                 0.5,  0.5,
                                -0.9999999, -0.9999999 };

    const double rhos[28]   = { 0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5, 0.5,
                                0.0, -0.5,
                                -0.9999999, -0.9999999 };
                                
    for (int i = 0; i < 28; i++)
    {
        const double bv_cnd = bivariate_standard_normal_cdf(xs[i], ys[i], rhos[i]);
//        std::cout << i << ": " << bv_cnd << std::endl;
    }
    
    /* Draw Bi-variate CDF */
//    const double min_x = -6.0f;
//    const double max_x =  6.0f;
//    const double min_y = -6.0f;
//    const double max_y =  6.0f;
//    const unsigned long steps = 100;
//    const double x_range = max_x - min_x;
//    const double y_range = max_y - min_y;
//    const double x_sample_width = x_range / static_cast<double>(steps);
//    const double y_sample_width = y_range / static_cast<double>(steps);
//
//    clock_t bv_cdf_init_time = clock();
//    for (double i = min_x + x_sample_width; i < max_x; i += x_sample_width)
//    {
//        for (double j = min_y + y_sample_width; j < max_y; j += y_sample_width)
//        {
//            double cdf_sample = bivariate_standard_normal_cdf(i, j, -0.99999);
//            std::cout << cdf_sample << " ";
//        }
//        std::cout << std::endl;
//    }
//    clock_t bv_cdf_end_time = clock();
//    cout << "Bi-variate CDF Start Time: " << bv_cdf_init_time << endl;
//    cout << "Bi-variate CDF Finish Time: " << bv_cdf_end_time << endl;
//    cout << "Bi-variate CDF Duration: " << bv_cdf_end_time - bv_cdf_init_time << endl;

    /* Draw Bi-variate PDF */
    const double min_x = -6.0f;
    const double max_x =  6.0f;
    const double min_y = -6.0f;
    const double max_y =  6.0f;
    const unsigned long steps = 100;
    const double x_range = max_x - min_x;
    const double y_range = max_y - min_y;
    const double x_sample_width = x_range / static_cast<double>(steps);
    const double y_sample_width = y_range / static_cast<double>(steps);

    clock_t bv_pdf_init_time = clock();
    for (double i = min_x + x_sample_width; i < max_x; i += x_sample_width)
    {
        for (double j = min_y + y_sample_width; j < max_y; j += y_sample_width)
        {
            double cdf_sample = bivariate_standard_normal_pdf(i, j, 0.0);
            std::cout << cdf_sample << " ";
        }
        std::cout << std::endl;
    }
//    clock_t bv_pdf_end_time = clock();
//    cout << "Bi-variate PDF Start Time: " << bv_pdf_init_time << endl;
//    cout << "Bi-variate PDF Finish Time: " << bv_pdf_end_time << endl;
//    cout << "Bi-variate PDF Duration: " << bv_pdf_end_time - bv_pdf_init_time << endl;
    
    
//    const int number_of_values = atoi(argv[1]);
//    
//    /* Generate mersenne twister numbers and time how long it takes */
//    int *mt_output_values = new int[number_of_values];
//    clock_t mt_init_time = clock();
//    mersenne_twister mt(0);
//    for (unsigned long i = 0; i < 62400000000000; i++)
//    {
//        mt_output_values[i % number_of_values] = mt.get_next();
//    }
//    clock_t mt_end_time = clock();
//    cout << "MT Start Time: " << mt_init_time << endl;
//    cout << "MT Finish Time: " << mt_end_time << endl;
//    cout << "MT Duration: " << mt_end_time - mt_init_time << endl;
//    
//    return 0;
//
//    /* Generate SIMD mersenne twister numbers and time how long it takes */
//    int *sfmt_output_values = new int[number_of_values];
//    clock_t sfmt_init_time = clock();
//    mersenne_twister sfmt(0);
//    for (int i = 0; i < number_of_values; i++)
//    {
//        sfmt_output_values[i] = sfmt.get_next();
//    }
//    clock_t sfmt_end_time = clock();
//    cout << "SFMT Start Time: " << sfmt_init_time << endl;
//    cout << "SFMT Finish Time: " << sfmt_end_time << endl;
//    cout << "SFMT Duration: " << sfmt_end_time - sfmt_init_time << endl;
//
//    
//    /* Check Black-Scholes price using inverse cdf algorithm */
//    const float ir = 0.05f;
//    const float vol = 0.2f;
//    const float s0 = 100.0f;
//    const float k = 100.0f;
//    const float t = 1.0f;
//    cout << "BS value: " << european_call_value(s0, k, vol, ir, t) << endl;
//
//
//    /* Plot a standard normal pdf and cdf */
//    const float min_x = 0.0f;
//    const float max_x = 1.0f;
//    const int steps = 100000;
//    const float range = max_x - min_x;
//    const float sample_width = range / static_cast<float>(steps);
//
//    float *cdf_values = new float [steps];
//    int cdf_sample = 0;
//    clock_t cdf_init_time = clock();
//    for (float i = min_x + sample_width; i < max_x; i+= sample_width)
//    {
//        cdf_values[cdf_sample++] = inverse_standard_normal_cdf(i, false);
//    }
//    clock_t cdf_end_time = clock();
//    cout << "CDF Start Time: " << cdf_init_time << endl;
//    cout << "CDF Finish Time: " << cdf_end_time << endl;
//    cout << "CDF Duration: " << cdf_end_time - cdf_init_time << endl;
//    
//    float *pdf_values = new float [steps];
//    int pdf_sample = 0;
//    clock_t pdf_init_time = clock();
//    for (float i = min_x; i < max_x; i+= sample_width)
//    {
//        pdf_values[pdf_sample++] = standard_normal_pdf(i);
//    }
//    clock_t pdf_end_time = clock();
//    cout << "PDF Start Time: " << pdf_init_time << endl;
//    cout << "PDF Finish Time: " << pdf_end_time << endl;
//    cout << "PDF Duration: " << pdf_end_time - pdf_init_time << endl;
//    
//    ofstream cdf_csv("inverse_cdf.csv");
//    ofstream pdf_csv("inverse_pdf.csv");
//    int dump_sample = 0;
//    for (float i = min_x + sample_width; i < max_x; i+= sample_width)
//    {
//        cdf_csv << i << ", " << cdf_values[dump_sample  ] << endl;
//        pdf_csv << i << ", " << pdf_values[dump_sample++] << endl;
//    }
//    cdf_csv.close();
//    pdf_csv.close();
//    cout << "Standard normals dumped" << endl;
//    
//
//    /* Test inverse cumulative normal */
//    const int icn_samples = atoi(argv[1]);
//    float *icn_values = new float [icn_samples];
//    mersenne_twister icn_mt(0);
//    for (int i = 0; i < icn_samples; i++)
//    {
//        icn_values[i] = inverse_standard_normal_cdf(icn_mt.get_next(), true);
//    }
//    
//    /* Test box muller with mersenne twister for generating normal numbers */
//    const int bm_samples = 10000000;
//    float* bm_values = new float[bm_samples];
//    box_muller *norm_rand = new box_muller(new mersenne_twister(0));
//    clock_t bm_init_time = clock();
//    for (int i = 0; i < bm_samples; i++)
//    {
//        bm_values[i] = norm_rand->get_next();
//    }
//    clock_t bm_end_time = clock();
//    cout << "BM Start Time: " << bm_init_time << endl;
//    cout << "BM Finish Time: " << bm_end_time << endl;
//    cout << "BM Duration: " << bm_end_time - bm_init_time << endl;
//    
//    /* Truncate and bin values */
//    const float min_bm = -6.0f;
//    const float max_bm =  6.0f;
//    const float bin_width_bm = 0.01f;
//    const int nr_bm_bins = (max_bm - min_bm) / bin_width_bm;
//    int *bm_bins = new int[nr_bm_bins];
//    memset(bm_bins, 0, nr_bm_bins * sizeof(int));
//    for (int i = 0; i < bm_samples; i++)
//    {
//        const float value = bm_values[i];
//        if ((value < min_bm) || (value > max_bm))
//        {
//            continue;
//        }
//        
//        bm_bins[(int)((value - min_bm) / bin_width_bm)]++;
//    }
//    
//    ofstream bm_csv("bm_values.csv");
//    for (int i = 0; i < nr_bm_bins; i++)
//    {
//        bm_csv << (min_bm + (bin_width_bm * 0.5f) + (bin_width_bm * (float)i)) << "," << (((float)bm_bins[i] / (float)bm_samples) / bin_width_bm)  << endl;
//    }
//    bm_csv.close();
//    
//    
//    /* Dump the random numbers to file */
//    if (false)
//    {
//        ofstream mt_numbers("mt_numbers.txt");
//        for (int i = 0; i < number_of_values; i++)
//        {
//            mt_numbers << i << ", " << mt_output_values[i] << endl;
//        }    
//        mt_numbers.close();
//        cout << "Reference values dumped" << endl;
//    }
//    /*  Check the results */
//    else if (true)
//    {
//        bool passed = true;
//        ifstream mt_numbers("mt_numbers.txt");
//        string expected;
//        for (int i = 0; i < number_of_values; i++)
//        {
//            stringstream actual;
//            actual << i << ", " << mt_output_values[i];
//            getline(mt_numbers, expected);
//            if (expected.compare(actual.str()) != 0)
//            {
//                passed &= false;
//                cout << "Mismatch on line: " << i << endl;
//                cout << "Expected: " << expected << " Actual: " << actual.str() << endl;
//            }
//        }    
//        mt_numbers.close();
//        
//        if (passed)
//        {
//            cout << "Test Passed" << endl;
//        }
//        else
//        {
//            cout << "Test Failed" << endl;
//        }
//    } 
//
//    delete [] mt_output_values;
//    delete [] sfmt_output_values;
//    delete [] cdf_values;
//    delete [] pdf_values;
//    
    return 1; 
}
