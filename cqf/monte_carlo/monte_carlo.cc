#include <algorithm>
#include <assert.h>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <string>

#include "monte_carlo.h"
#include "mersenne_twister.h"
#include "box_muller.h"


using std::cout;
using std::endl;


/* Use built in rand() to get a normally distributed random number */
float get_normal_rand()
{
    float sum = 0.0f;
    for (int i = 0; i < 12; i++)
    {
        sum += (float)rand() / (float)RAND_MAX;
    }
    
    return sum - 6.0f;
}


float get_normal_number(mersenne_twister *rand)
{
    float sum = 0.0;
    for (int i = 0; i < 12; i++)
    {
        sum += static_cast<float>(static_cast<unsigned>(rand->get_next())) / static_cast<float>(std::numeric_limits<int>::max());
    }
    
    return sum - 6.0f;
}


int main(int argc, char *argv[])
{
    /* Parse args */
    int mc_samples = 1000000;
    for (int i = 1; i < argc; i++) 
    {
        if (std::string(argv[i]) == "-mc_samples")
        {
            if ((i + 1) == argc)
            {
                cout << "Missing parameter for '-mc_samples'. Using default." << endl;
            }
            else
            {
                mc_samples = atoi(argv[i + 1]);
                ++i;
            }
        }
        else
        {
            cout << "Invalid arguement " << argv[i] << endl;
        }
    }


    /* Simulation parameters */
    const float mc_samples_inv = 1.0f / static_cast<float>(mc_samples);
    
    /* Economics */
    const float s0  = 100.0f;
    const float vol = 0.2f;
    const float ir  = 0.05f;
    const float t   = 1.0f;
    const float k   = 100.0f;
    
    /* Random number generator */
    box_muller *normal_dist = new box_muller(new mersenne_twister(5489));
    
    /* Simulations */
    int *sim_results = new int[mc_samples];
    for (int i = 0; i < mc_samples; i++)
    {
        /* Simulate to the end in one big hop */
        sim_results[i] = s0 * (1.0f + (ir * t) + (sqrt(t) * vol * normal_dist->get_next()));
    }
    
    
    /* Apply the payoff function */
    float * payoffs = new float[mc_samples];
    float average_payoff = 0.0f;
    for (int i = 0; i < mc_samples; i++)
    {
        /* Call option */
        payoffs[i] = std::max(sim_results[i] - k, 0.0f);
        average_payoff += payoffs[i] * mc_samples_inv;
    }
    
    
    /* Discount */
    float v = average_payoff * exp(-ir * t);
    
    
    /* Output results */
    //for (int i =0; i < mc_samples; i++)
    //{
        //cout << sim_results[i] << " " << payoffs[i] << " " << normal_dist->get_next() << endl;
    //}
    
    //cout << "Average payoff: " << average_payoff << endl;
    cout << "MC Samples: " << mc_samples << endl;
    cout << "V: " << v << endl;
    
    /* Clean up */
    delete [] sim_results;
    delete [] payoffs;
    
    return 1;
}
