#ifndef __RISK_FACTORY_H__
#define __RISK_FACTORY_H__

#include <initializer_list>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "risk.h"
#include "first_order_risk.h"
#include "second_order_risk.h"
#include "third_order_risk.h"


class compare_except_last_character
{
    public :
        bool operator()(const std::string &lhs, const std::string &rhs)const
        {
            return (lhs.compare(0, lhs.size() - 1, rhs) < 0);
        }
} except_last_character_comparator;



class equal_except_last_character
{
    public :
        bool operator()(const std::string &lhs, const std::string &rhs)const
        {
            return (lhs.compare(0, lhs.size() - 1, rhs, 0, rhs.size() - 1) == 0);
        }
} except_last_character_equal;


class compare_except_last_character_within
{
    public :
        compare_except_last_character_within() 
            : _valid_risk({ "rate", "spot", "time", "valu", "vol" }) {  };
            
        bool operator()(const std::string &lhs)const
        {
            for (unsigned int i = 0; i < _valid_risk.size(); i++)
            {
                if (lhs.compare(0, lhs.size() - 1, _valid_risk[i]) == 0)
                {
                    return false;
                }
            }
            
            return true;
        }
        
    private :
        const std::vector<std::string> _valid_risk;
} except_last_character_within_comparator;


class compare_empty_string
{
    public :
        bool operator()(const std::string &lhs) const
        {
            return lhs.empty();
        }
} empty_string_comparator;


/* Class from parsing risk definitions as strings into risk classes */
template<class T>
class risk_factory
{
    public :
        /* Function to build risks from strings */
        risk<T>* build(const std::vector<std::string> &def) const
        {
            typedef std::vector<std::string>::iterator split_def_iter;
            std::vector<std::string> split_def;
            boost::split(split_def, def[0], boost::is_any_of("d"));

            split_def_iter remove_from = std::remove_if(split_def.begin(), split_def.end(), empty_string_comparator);
            const int split_def_size = std::distance(split_def.begin(), remove_from);
            assert((def.size() == static_cast<unsigned int>((split_def_size * 3) + 4)) || !"Error: Incorrect number of inputs for risk definition.");
            
            /* Check for invalid risk */
            split_def_iter invalid = std::find_if(split_def.begin(), remove_from, except_last_character_within_comparator);
            if (invalid != remove_from)
            {
                std::cout << "Error: Found invalid risk - " << (*invalid) << std::endl;
                assert(false);
            }
            
            /* Check for duplicate risks ignoring order */
            std::sort(split_def.begin(), remove_from, except_last_character_comparator);
            split_def_iter notunique = std::unique(split_def.begin(), remove_from, except_last_character_equal);
            if (notunique != remove_from)
            {
                for ( ; notunique != remove_from; notunique++)
                {
                    std::cout << "Error: Found duplicate risk " << (*notunique) << ". Please merge common risks." << std::endl;
                }
                assert(false);
            }
            
            /* Get log file or risk name */
            std::string log_file;
            if (def[1].compare(".") != 0)
            {
                log_file = def[1];
            }

            /* Get weather to dump grid or single risk values */
            bool grid_risk = false;
            if (def[2].compare("grid") == 0)
            {
                grid_risk = true;
            }
            
            T disc_grid_offset = 0.0;
            if (def[3].compare(".") != 0)
            {
                disc_grid_offset = boost::lexical_cast<T>(def[3]);
            }
            
            /* Build the, always present, leaf node */
            std::string empty;
            const bool root_node = (def.size() == 7);
            const std::string &this_log = root_node ? log_file : empty;
            risk<T> *r = new risk_node<T>(create_scenario_bumper(&def[def.size() - 3], (*(--remove_from))), this_log, disc_grid_offset, root_node, grid_risk);
                
            /* Build the composite risks, if any */
            for (int i = split_def_size - 3; i >= 0; i--)
            {
                const bool root_node = (i == 0);
                const std::string &this_log = root_node ? log_file : empty;
                r = new composite_risk<T>(create_scenario_bumper(&def[4 + (i * 3)], (*(--remove_from))), this_log, r, disc_grid_offset, root_node, grid_risk);
            }
            
            return r;
        }
        
    private :
        bump_dir_t get_bump_dir(const std::string &option) const
        {
            if (option.compare("DOWN") == 0)
            {
                return DOWN;
            }
            else if (option.compare("UP") == 0)
            {
                return UP;
            }
            else if (option.compare("BI") == 0)
            {
                return BI;
            }
            else
            {
                std::cout << "Error: Unknown bump direction " << option << std::endl;
                assert(false);
            }
        }


        bump_style_t get_bump_style(const std::string &option) const
        {
            if (option.compare("ADD") == 0)
            {
                return ADD;
            }
            else if (option.compare("MULT") == 0)
            {
                return MULT;
            }
            else if (option.compare("OVERRIDE") == 0)
            {
                return OVERRIDE;
            }
            else
            {
                std::cout << "Error: Unknown bump style " << option << std::endl;
                assert(false);
            }
        }


        scenario_bumper<T>* create_rate_bumper(const std::string &risk_def, const T bump_size, const bump_dir_t bump_dir, const bump_style_t bump_style) const
        {
            switch (risk_def.data()[risk_def.length() - 1])
            {
                case '1' : 
                    return new drate1_scenario_bumper<T>(1, bump_size, bump_dir, bump_style);

                case '2' : 
                    return new drate2_scenario_bumper<T>(2, bump_size, bump_dir, bump_style);

                case '3' : 
                    return new drate3_scenario_bumper<T>(3, bump_size, bump_dir, bump_style);
                    
                default :
                    std::cout << "Error: Rate risk of unknown order " << risk_def.data()[risk_def.length()] << std::endl;
                    assert(false);
            }
        }
        
        scenario_bumper<T>* create_spot_bumper(const std::string &risk_def, const T bump_size, const bump_dir_t bump_dir, const bump_style_t bump_style) const
        {
            switch (risk_def.data()[risk_def.length() - 1])
            {
                case '1' : 
                    return new dspot1_scenario_bumper<T>(1, bump_size, bump_dir, bump_style);

                case '2' : 
                    return new dspot2_scenario_bumper<T>(2, bump_size, bump_dir, bump_style);

                case '3' : 
                    return new dspot3_scenario_bumper<T>(3, bump_size, bump_dir, bump_style);
                    
                default :
                    std::cout << "Error: Spot risk of unknown order " << risk_def.data()[risk_def.length()] << std::endl;
                    assert(false);
            }
        }
        
        scenario_bumper<T>* create_time_bumper(const std::string &risk_def, const T bump_size, const bump_dir_t bump_dir, const bump_style_t bump_style) const
        {
            switch (risk_def.data()[risk_def.length() - 1])
            {
                case '1' : 
                    return new dtime1_scenario_bumper<T>(1, bump_size, bump_dir, bump_style);

                case '2' : 
                    return new dtime2_scenario_bumper<T>(2, bump_size, bump_dir, bump_style);

                case '3' : 
                    return new dtime3_scenario_bumper<T>(3, bump_size, bump_dir, bump_style);
                    
                default :
                    std::cout << "Error: Time risk of unknown order " << risk_def.data()[risk_def.length()] << std::endl;
                    assert(false);
            }
        }
        
        scenario_bumper<T>* create_vol_bumper(const std::string &risk_def, const T bump_size, const bump_dir_t bump_dir, const bump_style_t bump_style) const
        {
            switch (risk_def.data()[risk_def.length() - 1])
            {
                case '1' : 
                    return new dvol1_scenario_bumper<T>(1, bump_size, bump_dir, bump_style);

                case '2' : 
                    return new dvol2_scenario_bumper<T>(2, bump_size, bump_dir, bump_style);

                case '3' : 
                    return new dvol3_scenario_bumper<T>(3, bump_size, bump_dir, bump_style);
                    
                default :
                    std::cout << "Error: Vol risk of unknown order " << risk_def.data()[risk_def.length()] << std::endl;
                    assert(false);
            }
        }
        
        scenario_bumper<T>* create_scenario_bumper(const std::string *const def, const std::string &risk) const
        {
            /* Defaults */
            T bump_size             = 0.01; 
            bump_dir_t bump_dir     = BI;
            bump_style_t bump_style = ADD;
         
            /* Parse params, accept default if . */
            if (def[0].compare(".") != 0)
            {
                bump_size = boost::lexical_cast<T>(def[0]);
            }
         
            if (def[1].compare(".") != 0)
            {
                bump_dir = get_bump_dir(def[1]);
            }            
         
            if (def[2].compare(".") != 0)
            {
                bump_style = get_bump_style(def[2]);
            }
         
            /* Construct and return the bumper */
            if (except_last_character_equal(risk, "valux"))
            {
                return new value_scenario_bumper<T>();
            }
            else if (except_last_character_equal(risk, "volx"))
            {
                return create_vol_bumper(risk, bump_size, bump_dir, bump_style);
            }
            else if (except_last_character_equal(risk, "ratex"))
            {
                return create_rate_bumper(risk, bump_size, bump_dir, bump_style);
            }
            else if (except_last_character_equal(risk, "timex"))
            {
                return create_time_bumper(risk, bump_size, bump_dir, bump_style);
            }
            else if (except_last_character_equal(risk, "spotx"))
            {
                return create_spot_bumper(risk, bump_size, bump_dir, bump_style);
            }
            else
            {
                std::cout << "Error: Unknown risk " << risk << std::endl;
                assert(false);
            }
        }
};

#endif
