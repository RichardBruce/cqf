#include "halton_numbers.h"
#include "sobol_numbers.h"


double *sobol_points(unsigned n, unsigned D, char *dir_file)
{
    std::ifstream infile(dir_file, std::ios::in);
    if (!infile)
    {
        std::cout << "Input file containing direction numbers cannot be found!\n";
        exit(1);
    }
 
    // L = max number of bits needed
    int L = (int)ceil(log((double)n)/log(2.0));
 
    // C[i] = index from the right of the first zero bit of i
    double *numbers = new double [n * D];
    unsigned *C = new unsigned [n];
    unsigned *V = new unsigned [L + 1];
    unsigned *X = new unsigned [n];
    C[0] = 1;
    for (unsigned int i = 1; i < n; i++)
    {
        C[i] = 1;
        int value = i;
        while (value & 0x1)
        {
            value >>= 1;
            C[i]++;
        }
    }
 
    // POINTS[i][j] = the jth component of the ith point
    // with i indexed from 0 to n-1 and j indexed from 0 to D-1
    // ----- Compute the first dimension -----
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    for (int i = 1; i <= L; i++)
    {
        V[i] = 1 << (32 - i);
    }
 
    // all m's = 1
    // Evalulate X[0] to X[n-1], scaled by pow(2,32)
    X[0] = 0;
    for (unsigned int i = 1; i < n; i++)
    {
        X[i] = X[i-1] ^ V[C[i-1]];
        numbers[i] = (double)X[i]/pow(2.0,32);
    }
 
    // ----- Compute the remaining dimensions -----
    for (unsigned int i = 1; i < D; i++)
    {
        // Read in parameters from file
        int d, s, a;
        infile >> d >> s >> a;
        int *m = new int [s + 1];
        for (int j = 1; j <= s; j++)
        {
            infile >> m[j];
        }
 
        // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
        if (L <= s)
        {
            for (int j = 1; j <= L; j++)
            {
                V[j] = m[j] << (32 - j);
            }
        }
        else
        {
            for (int j = 1; j <= s; j++)
            {
                V[j] = m[j] << (32 - j);
            }
 
            for (int j = s + 1; j <= L; j++)
            {
                V[j] = V[j - s] ^ (V[j - s] >> s);
                for (int k = 1; k < s; k++)
                {
                    V[j] ^= (((a >> (s - 1 - k)) & 1) * V[j - k]);
                }
            }
        }
 
        // Evalulate X[0] to X[n-1], scaled by pow(2,32)
        X[0] = 0;
        for (unsigned int j = 1; j < n; j++)
        {
            X[j] = X[j - 1] ^ V[C[j - 1]];
            numbers[j + (i * n)] = (double)X[j] / pow(2.0,32);
        }
 
        // Clean up
        delete [] m;
    }
 
    /* Clean up */
    delete [] C;
    delete [] V;
    delete [] X;
    
    return numbers;
}
 
int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << std::endl << "input format: sobol n D FILENAME" << std::endl << std::endl;
        std::cout << "The program prints the first n sobol points in D dimensions." << std::endl;
        std::cout << "The points are generated in graycode order." << std::endl;
        std::cout << "The primitive polynomials and initial direction numbers are" << std::endl << "given by the input file FILENAME." << std::endl << std::endl;
        return 0;
    }
 
    const int n = atoi(argv[1]);
    const int d = atoi(argv[2]);

#if 0
    halton_numbers hn(n, std::min(d, 10));
    
    /* Dump the low descrepancy numbers to file */
    if (false)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < std::min(d, 10); j++)
            {
                std::cout << hn.get_number(i, j) << ", ";
            }
            std::cout << std::endl;
        }
    }
    /* Check the results */
    else if (true)
    {
        bool passed = true;
        std::ifstream halton_numbers("halton_numbers.txt");
        std::string expected;
        for (int i = 0; i < n; i++)
        {
            std::stringstream actual;
            for (int j = 0; j < std::min(d, 10); j++)
            {
                actual << hn.get_number(i, j) << ", ";
            }
            
            getline(halton_numbers, expected);
            if (expected.compare(actual.str()) != 0)
            {
                passed &= false;
                std::cout << "Halton mismatch on line: " << i << std::endl;
                std::cout << "Expected: " << expected << " Actual: " << actual.str() << std::endl;
            }
        }
        halton_numbers.close();
        
        if (passed)
        {
            std::cout << "Test Passed" << std::endl;
        }
        else
        {
            std::cout << "Test Failed" << std::endl;
        }
    }
    
#else
    // display points std::cout << setprecision(20);
    //std::cout << setiosflags(ios::scientific) << setprecision(10);
    /* Dump the low descrepancy numbers to file */
    float *nums = new float [d];
    const int skip = 12439;
    sobol_numbers<float> sn0(argv[3], n, d);
    sobol_numbers<float> sn1(sn0);
    if (false)
    {
        for (int i = 0; i< n; i++)
        {
            sn0.get_next(nums);
            for (int j = 0; j < d; j++)
            {
                std::cout << nums[j] << ", " ;
            }
            std::cout << std::endl;
        }
    }
    /* Check the results */
    else if (true)
    {    
        bool passed = true;
        std::ifstream sobol_numbers("sobol_numbers.txt");
        std::string expected;
        for (int i = 0; i < n; i++)
        {
            std::stringstream actual;
            if (i < skip)
            {
                sn0.get_next(nums);
            }
            else if (i == skip)
            {
                sn1.skip(nums, skip);
            }
            else
            {
                sn1.get_next(nums);
            }
            
            for (int j = 0; j < d; j++)
            {
                actual << nums[j] << ", ";
            }
        
            getline(sobol_numbers, expected);
            if (expected.compare(actual.str()) != 0)
            {
                passed &= false;
                std::cout << "Sobol mismatch on line: " << i << std::endl;
                std::cout << "Expected: " << expected << std::endl;
                std::cout << "Actual:   " << actual.str() << std::endl;
            }
        }    
        sobol_numbers.close();
        
        if (passed)
        {
            std::cout << "Test Passed" << std::endl;
        }
        else
        {
            std::cout << "Test Failed" << std::endl;
        }
    }
    delete [] nums;
#endif
}
