#ifndef __CS_ADI_2D_H__
#define __CS_ADI_2D_H__

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <assert.h>

#include "pde_solver.h"
#include "correctors.h"

#include "cashflow.h"

#include "scenario_point.h"


template<class T>
class cs_adi_pde_solver : public pde_solver<T>
{
    public :
        cs_adi_pde_solver() { };
         
        void solve(scenario_point<T> &p);
        
    private :
        std::vector<T>  _tmp;
};

    
template<class T>    
void cs_adi_2d_pde_solver<T>::solve(scenario_point<T> &p)
{
    grid<T> &grid = p.get_grid();
    assert((grid.dimensions() > 2) || (grid.dimensions() == 2 && !grid.outer_grid_disconnected()));
    
    /* Resize tmp storage */
    if (_tmp.size() < static_cast<unsigned int>(grid.total_geometric_size()))
    {
        _tmp.resize(grid.total_geometric_size());
    }
    
    const T t_inc   = p.get_time_increment();
    const T theta   = 0.5; 
    
    const int dim_m1      = grid.dimensions() - 1;
    const int nr_of_grids = grid.outer_grid_disconnected() ? grid.size(0) : 1;
    
    explicit_predictor<T> pred(theta, t_inc, sigma_s1, sigma_s2, r_s1, r_s2, this->_corr);
    implicit_corrector<T> x_corr(r_s1, theta, s1_inc, t_inc, sigma_s1, this->_nas2, this->_nas1);
    implicit_corrector<T> y_corr(r_s2, theta, s2_inc, t_inc, sigma_s2, 1, this->_nas2);
    cross_term_corrector<T> xy_corr(0.5, t_inc, s1_inc, s2_inc, sigma_s1, sigma_s2, this->_corr, this->_nas2);
    for (int i = 0; i < nr_of_grids; i++)
    {
        /* Explicit predictor step */
        pred(grid.read_v_grid(), grid.write_v_grid());

        /* S1 direction corrector */
        for (int j = 0; j < grid.size(0); j++)
        {
            x_corr(&grid.read_v_grid()[j], &grid.write_v_grid()[j], &_tmp[j]);
        }
                
        /* S2 direction corrector */
        for (int j = 0; j < grid.size(1); j++)
        {
            y_corr(&grid.read_v_grid()[(this->_nas2 * j)], &_tmp[(this->_nas2 * j)], &_tmp[(this->_nas2 * j)]);
        }
        
        /* Explicit cross term corrector */
        for (int j = 1; j < grid.size(1) - 1; j++)
        {
            xy_corr(&_tmp[0], grid.read_v_grid(), grid.write_v_grid(), j);
        }
    
        /* Second S1 direction corrector */
        for (int j = 0; j < grid.size(0); j++)
        {
            x_corr(&grid.read_v_grid()[j], &grid.write_v_grid()[j], &grid.write_v_grid()[j]);
        }
                       
        /* Second S2 direction corrector */
        for (int j = 0; j < grid.size(0); j++)
        {
            y_corr(&grid.read_v_grid()[(this->_nas2 * j)], &grid.write_v_grid()[(this->_nas2 * j)], &grid.write_v_grid()[(this->_nas2 * j)]);
        }
    }
}

#endif
