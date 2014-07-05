#include <iostream>

#include "implied_volatility.h"

using namespace std;


int main()
{
    implied_volatility *imp_vol = new implied_volatility(0.0001, 10);
    const float vol = imp_vol->back_out_from_call(100.0, 100.0, 0.05, 1, 10.40058);
    
    cout << "Implied vol: " << vol << endl;
    
    /* Clean up */
    delete imp_vol;
    
    return 0;
}
