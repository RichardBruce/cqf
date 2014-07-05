#ifndef __CASHFLOW_FACTORY_H__
#define __CASHFLOW_FACTORY_H__


#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "cashflow.h"
#include "barrier_cashflow.h"
#include "single_asset_cashflow.h"


template<class T>
class cashflow_factory
{
    public :
        
        cashflow<T>* build(const std::vector<std::string> &def, const T spot) const
        {
            /* const T s2_inc = s2_max / static_cast<T>(nas2);
            if (strcmp(argv[i], "wo_call") == 0)
            {
                const T k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const T k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                return new two_asset_payoff_cashflow<T, wo_call_cashflow_func<T> >(wo_call_cashflow_func<T>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2);
            }
            else if (strcmp(argv[i], "bo_call") == 0)
            {
                const T k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const T k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                return new two_asset_payoff_cashflow<T, bo_call_cashflow_func<T> >(bo_call_cashflow_func<T>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2);
            }
            else if (strcmp(argv[i], "wo_binary_call") == 0)
            {
                const T k_s1 = get_arg_float("-cashflow", argv, argc, ++i);
                const T k_s2 = get_arg_float("-cashflow", argv, argc, ++i);
                return new two_asset_payoff_cashflow<T, wo_binary_call_cashflow_func<T> >(wo_binary_call_cashflow_func<T>(), k_s1, k_s2, s1_inc, s2_inc, nas1, nas2);
            }
            else */if (def[0].compare("call") == 0)
            {
                const T k_s1 = boost::lexical_cast<T>(def[2]);
                const barrier_cashflow<T> barr = parse_barrier(&def[4]);
                if (def[1].compare(0, 5, "asian") == 0)
                {
                    const bool asian_in = (def[1].compare("asian_in") == 0);
                    std::vector<T> mat = parse_asian_maturities(def[3]);
                    return new one_asset_asian_payoff_cashflow<T, euro_call_cashflow_func<T>>(euro_call_cashflow_func<T>(barr, k_s1), mat, spot, asian_in);
                }
                else
                {
                    const T mat = boost::lexical_cast<T>(def[3]);
                    if (def[1].compare("amer") == 0)
                    {
                        return new one_asset_payoff_cashflow<T, amer_call_cashflow_func<T>>(amer_call_cashflow_func<T>(barr, k_s1), mat, spot);
                    }
                    else
                    {
                        return new one_asset_payoff_cashflow<T, euro_call_cashflow_func<T>>(euro_call_cashflow_func<T>(barr, k_s1), mat, spot);
                    }
                }
            }
            else if (def[0].compare("put") == 0)
            {
                const T k_s1 = boost::lexical_cast<T>(def[2]);
                const barrier_cashflow<T> barr = parse_barrier(&def[4]);
                if (def[1].compare(0, 5, "asian") == 0)
                {
                    const bool asian_in = (def[1].compare("asian_in") == 0);
                    std::vector<T> mat = parse_asian_maturities(def[3]);
                    return new one_asset_asian_payoff_cashflow<T, euro_put_cashflow_func<T>>(euro_put_cashflow_func<T>(barr, k_s1), mat, spot, asian_in);
                }
                else
                {
                    const T mat = boost::lexical_cast<T>(def[3]);
                    if (def[1].compare("amer") == 0)
                    {
                        return new one_asset_payoff_cashflow<T, amer_put_cashflow_func<T> >(amer_put_cashflow_func<T>(barr, k_s1), mat, spot);
                    }
                    else
                    {
                        return new one_asset_payoff_cashflow<T, euro_put_cashflow_func<T> >(euro_put_cashflow_func<T>(barr, k_s1), mat, spot);
                    }
                }
            }
            else if (def[0].compare("binary_call") == 0)
            {
                const bool  amer = (def[1].compare("amer") == 0);
                const T     k_s1 = boost::lexical_cast<T>(def[2]);
                const T     mat  = boost::lexical_cast<T>(def[3]);
                const barrier_cashflow<T> barr = parse_barrier(&def[4]);
                if (amer)
                {
                    return new one_asset_payoff_cashflow<T, amer_binary_call_cashflow_func<T> >(amer_binary_call_cashflow_func<T>(barr, k_s1), mat, spot);
                }
                else
                {
                    return new one_asset_payoff_cashflow<T, euro_binary_call_cashflow_func<T> >(euro_binary_call_cashflow_func<T>(barr, k_s1), mat, spot);
                }
            }
            else if (def[0].compare("binary_put") == 0)
            {
                const bool  amer = (def[1].compare("amer") == 0);
                const T     k_s1 = boost::lexical_cast<T>(def[2]);
                const T     mat  = boost::lexical_cast<T>(def[3]);
                const barrier_cashflow<T> barr = parse_barrier(&def[4]);
                if (amer)
                {
                    return new one_asset_payoff_cashflow<T, amer_binary_put_cashflow_func<T> >(amer_binary_put_cashflow_func<T>(barr, k_s1), mat, spot);
                }
                else
                {
                    return new one_asset_payoff_cashflow<T, euro_binary_put_cashflow_func<T> >(euro_binary_put_cashflow_func<T>(barr, k_s1), mat, spot);
                }
            }
            else if (def[0].compare("dividend") == 0)
            {
                const bool  fix = (def[1].compare("fix") == 0);
                const T     div = boost::lexical_cast<T>(def[2]);
                const T     mat = boost::lexical_cast<T>(def[3]);
                if (fix)
                {
                    return new one_asset_dividend_cashflow<T, fixed_dividend<T> >(fixed_dividend<T>(div), mat, spot);
                }
                else
                {
                    return new one_asset_dividend_cashflow<T, percent_dividend<T> >(percent_dividend<T>(div), mat, spot);
                }
            }
            else if (def[0].compare("payment") == 0)
            {
                const bool  fix = (def[1].compare("fix") == 0);
                const T     pay = boost::lexical_cast<T>(def[2]);
                const T     mat = boost::lexical_cast<T>(def[3]);
                const barrier_cashflow<T> barr = parse_barrier(&def[4]);
                if (fix)
                {
                    return new one_asset_payoff_cashflow<T, fixed_payment<T> >(fixed_payment<T>(barr, pay), mat, spot);
                }
                else
                {
                    return new one_asset_payoff_cashflow<T, percent_payment<T> >(percent_payment<T>(barr, pay), mat, spot);
                }
            }
            else
            {
                std::cout << "Error: Unknown cashflow " << def[0] << std::endl;
                assert(false);
            }
        }
    
    private :
        std::vector<T> parse_asian_maturities(const std::string &mats) const
        {
            /* Clean and split out the maturities */
            std::vector<std::string> mat_str;
            boost::split(mat_str, mats, boost::is_any_of(","));
           
            /* Convert */
            std::vector<T> mat;
            for (unsigned int i = 0; i < mat_str.size(); i++)
            {
                mat.push_back(boost::lexical_cast<T>(mat_str[i]));
            }
           
            std::sort(mat.begin(), mat.end());
            typename std::vector<T>::iterator uni_iter = std::unique(mat.begin(), mat.end());
            mat.erase(uni_iter, mat.end());
            
            return mat;
        }
        
        barrier_dir_t get_barrier_dir(const std::string &option) const
        {
            if (option.compare("DOWN") == 0)
            {
                return DOWN_BARRIER;
            }
            else if (option.compare("UP") == 0)
            {
                return UP_BARRIER;
            }
            else
            {
                std::cout << "Error: Unknown barrier direction " << option << std::endl;
                assert(false);
            }
        }
        
        barrier_type_t get_barrier_type(const std::string &option) const
        {
            if (option.compare("IN") == 0)
            {
                return IN_BARRIER;
            }
            else if (option.compare("OUT") == 0)
            {
                return OUT_BARRIER;
            }
            else
            {
                std::cout << "Error: Unknown barrier type " << option << std::endl;
                assert(false);
            }
        }
        
        barrier_cashflow<T> parse_barrier(const std::string *const def) const
        {
            /* Parse params, accept default if . */
            bool amer = false;
            if (def[0].compare("amer") == 0)
            {
                amer = true;
            }
            
            T b = -1.0;
            if (def[1].compare(".") != 0)
            {
                b = boost::lexical_cast<T>(def[1]);
            }
         
            barrier_dir_t dir = NO_BARRIER_DIR;
            if (def[2].compare(".") != 0)
            {
                dir = get_barrier_dir(def[2]);
            }            
         
            barrier_type_t type = NO_BARRIER_TYPE;
            if (def[3].compare(".") != 0)
            {
                type = get_barrier_type(def[3]);
            }
            
            return barrier_cashflow<T>(b, dir, type, amer);
        }
};

#endif
