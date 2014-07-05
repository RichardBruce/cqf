#ifndef __GRID_FACTORY_H__
#define __GRID_FACTORY_H__

#include "grid_mapper.h"
#include "grid.h"


template<class T>
class grid_factory
{
    public :
        grid<T>* build(const std::vector<std::string> &def, std::vector<cashflow<T>*> &cashflows, const T spot) const
        {
            /* Get the cashflows dimensionality */
            std::vector<bool> dim_spec;
            for (unsigned int i = 0; i < cashflows.size(); i++)
            {
                cashflows[i]->dimensions(&dim_spec);
            }
        
            /* Check for crazies */
            assert((dim_spec.size() < 4) || ((dim_spec.size() < 5) && dim_spec[4]) || !"Error: Are you nuts 4d PDE, go use MC");
            assert((dim_spec.size() < 2) || ((dim_spec.size() < 3) && dim_spec[1]) || !"Error: Only 1 dimension (plus asianing) is supported so far.");
        
            /* Parse a mapper per dimension */
            unsigned int input_idx = 0;
            std::vector<boost::shared_ptr<const grid_mapper<T>>> mappers;
            boost::shared_ptr<const grid_mapper<T>> mapper;
            for (unsigned int i = dim_spec.size() - 1; i < dim_spec.size(); i--)
            {
                assert(!dim_spec[i] || (i == (dim_spec.size() - 1)) || !"Error: Only the outer dimension may be disconnected.");

                /* inc in intepreted as a maximum for uniform and an average for sinh */
                const double max_fac    = boost::lexical_cast<T>(def[input_idx + 1]);
                const double min_fac    = boost::lexical_cast<T>(def[input_idx + 2]);
                const double inc        = boost::lexical_cast<T>(def[input_idx + 3]);
                if (def[input_idx].compare("uniform") == 0)
                {
                    assert((def.size() >= (input_idx + 3)) || !"Error: Not enough grids input for payoff.");
                    input_idx += 4;
                    mapper.reset(new uniformed_grid_mapper<double>(max_fac, min_fac, inc));
                    mappers.push_back(mapper);
                }
                else if (def[input_idx].compare("sinh") == 0)
                {
                    assert((def.size() >= (input_idx + 5)) || !"Error: Not enough grids inputs for payoff.");
                    const double grid_ratio     = boost::lexical_cast<T>(def[input_idx + 4]);
                    const double grid_center    = boost::lexical_cast<T>(def[input_idx + 5]);
                    input_idx += 6;
                    mapper.reset(new sinh_grid_mapper<double>(grid_center, grid_ratio, max_fac, min_fac, inc));
                    mappers.push_back(mapper);
                }
                else
                {
                    std::cout << "Error: Unknown grid type " << def[0] << std::endl;
                    assert(false);
                }
            }
            
            return new grid<T>(mappers, cashflows, spot, dim_spec.back());
        }
};

#endif
