#include <iostream>

#include "european_call.h"
#include "european_put.h"
#include "european_fwd_start_call.h"
#include "european_fwd_start_put.h"
#include "european_cash_binary_call.h"
#include "european_cash_binary_put.h"
#include "european_two_asset_cash_binary.h"
#include "european_min_of_two_call.h"


using namespace std;

int main()
{
    double v;
//    v = european_cash_binary_call_value(100.0f, 100.0f, 0.2f, 0.05f, 1.0f);
//    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, -0.5, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, -0.5, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, -0.5, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, -0.5, 10.0, false, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.0, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.0, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.0, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.0, 10.0, false, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.5, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.5, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.5, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 0.5, 0.5, 10.0, false, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, -0.5, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, -0.5, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, -0.5, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, -0.5, 10.0, false, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.0, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.0, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.0, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.0, 10.0, false, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.5, 10.0, true, true);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.5, 10.0, false, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.5, 10.0, true, false);
    cout << v << endl;

    v = european_two_asset_cash_binary_value(100.0, 110.0, 0.2, 0.05, 100.0, 90.0, 0.25, 0.06,
        0.1, 1.0, 0.5, 10.0, false, true);
    cout << v << endl;
    
    v = european_two_asset_cash_binary_value(
        100.0, 100.0, 0.2, 0.05, 
        100.0, 100.0, 0.2, 0.05,
        0.05, 1.0, 0.5, 1.0, true, true);
    cout << v << endl;
    
    v = european_min_of_two_call_value(
    100.0, 0.2,  0.0,
    100.0, 0.2,  0.0, 
    0.05,  1.0,  0.5, 100.0);
    cout << v << endl;
    
    v = european_call_value(100.0, 100.0, 0.2, 0.05, 1.0);
    std::cout << "Call: " << v << std::endl;

    v = european_put_value(100.0, 100.0, 0.2, 0.05, 1.0);
    std::cout << "Put: " << v << std::endl;
    
    v = european_cash_binary_call_value(100.0, 100.0, 0.2, 0.05, 1.0);
    std::cout << "Binary Call: " << v << std::endl;

    v = european_cash_binary_put_value(100.0, 100.0, 0.2, 0.05, 1.0);
    std::cout << "Binary Put: " << v << std::endl;
    
    v = european_fwd_start_call_value(60.0, 1.1, 0.3, 0.08, (0.08 - 0.04), 0.25, 1.0);
    std::cout << "Forward Start Call: " << v << std::endl;
    
    v = european_fwd_start_put_value(60.0, 1.1, 0.3, 0.08, (0.08 - 0.04), 0.25, 1.0);
    std::cout << "Forward Start Put: " << v << std::endl;
}
