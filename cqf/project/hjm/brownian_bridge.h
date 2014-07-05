#ifndef __BROWNIAN_BRIDGE_H__
#define __BROWNIAN_BRIDGE_H__


template<class T>
class brownian_bridge
{
    public :
        /* CTOR */
        brownian_bridge(const unsigned int points) 
        : l_wgt(new T [points]), r_wgt(new T [points]), sigma(new T [points]), tmp(new T [points]), 
          l_idx(new unsigned int [points]), idx(new unsigned int [points]), 
          r_idx(new unsigned int [points]), points(points)
        {
            unsigned int *mapped = new unsigned int [points];
            memset(&mapped[0], 0, points * sizeof(unsigned int));
            
            mapped[points - 1] = 1;
            idx[0]      = points - 1;
            l_idx[0]    = 0;
            r_idx[0]    = 0;
            sigma[0]    = std::sqrt(points);
            l_wgt[0]    = 0.0;
            r_wgt[0]    = 0.0;
            
            unsigned int j = 0;
            for (unsigned int i = 1; i < points; i++)
            {
                /* Find unmapped point */
                while (mapped[j])
                {
                    ++j;
                }
                
                /* Find the next mapped points */
                unsigned int k = j;
                while (!mapped[k])
                {
                    ++k;
                }

                /* Save indexes of left, right and center */                
                unsigned int c = j + ((k - 1 - j) >> 1);
                mapped[c]   = i;
                l_idx[i]    = j;
                idx[i]      = c;
                r_idx[i]    = k;
                
                /* Save weightings for left and right and variance of current */
                const T c_p1_mj    = static_cast<T>(c + 1 - j);
                const T k_p1_mj    = static_cast<T>(k + 1 - j);
                const T k_mc       = static_cast<T>(k - c);
                l_wgt[i] = k_mc / k_p1_mj;
                r_wgt[i] = c_p1_mj / k_p1_mj;
                sigma[i] = std::sqrt((c_p1_mj * k_mc) / k_p1_mj);

                /* Advance to look for unmapped points in the next space */
                j = k + 1;
                if (j >= points)
                {
                    j = 0;
                }
            }
            
            /* Clean up */
            delete [] mapped;
        }


        /* CTOR */
        brownian_bridge(const std::vector<double> &dates, const double vol, const unsigned int points) 
        : l_wgt(new T [points]), r_wgt(new T [points]), sigma(new T [points]), tmp(new T [points]),
          l_idx(new unsigned int [points]), idx(new unsigned int [points]), 
          r_idx(new unsigned int [points]), points(points)
        {
            unsigned int *mapped = new unsigned int [points];
            memset(&mapped[0], 0, points * sizeof(unsigned int));
            
            double *v = new double [dates.size()];
            const double vol_sq = vol * vol;
            for (unsigned int i = 0; i < dates.size(); i++)
            {
                v[i] = dates[i] * vol_sq;
            }
            
            mapped[points - 1] = 1;
            idx[0]      = points - 1;
            l_idx[0]    = 0;
            r_idx[0]    = 0;
            sigma[0]    = std::sqrt(v[points - 1]);
            l_wgt[0]    = 0.0;
            r_wgt[0]    = 0.0;
            
            unsigned int j = 0;
            for (unsigned int i = 1; i < points; i++)
            {
                /* Find unmapped point */
                while (mapped[j])
                {
                    ++j;
                }
                
                /* Find the next mapped points */
                unsigned int k = j;
                while (!mapped[k])
                {
                    ++k;
                }

                /* Save indexes of left, right and center */                
                unsigned int c = j + ((k - 1 - j) >> 1);
                mapped[c]   = i;
                l_idx[i]    = j;
                idx[i]      = c;
                r_idx[i]    = k;
                
                /* Save weightings for left and right and variance of current */
//                const T c_p1_mj    = static_cast<T>(c + 1 - j) * 0.01;
//                const T k_p1_mj    = static_cast<T>(k + 1 - j) * 0.01;
//                const T k_mc       = static_cast<T>(k - c) * 0.01;
                l_wgt[i] = (v[k + 1] - v[c + 1]) / (v[k + 1] - v[j]);//k_mc / k_p1_mj;
                r_wgt[i] = (v[c + 1] - v[j]) / (v[k + 1] - v[j]);//c_p1_mj / k_p1_mj;
                sigma[i] = std::sqrt(((v[c + 1] - v[j]) * (v[k + 1] - v[c + 1])) / (v[k + 1] - v[j]));//std::sqrt((c_p1_mj * k_mc) / k_p1_mj);

                /* Advance to look for unmapped points in the next space */
                j = k + 1;
                if (j >= points)
                {
                    j = 0;
                }
            }
            
            /* Clean up */
            delete [] mapped;
            delete [] v;
        }
        
        /* Copy CTOR */
        brownian_bridge(const brownian_bridge &b) 
        : l_wgt(new T [b.points]), r_wgt(new T [b.points]), sigma(new T [b.points]), tmp(new T [b.points]),
          l_idx(new unsigned int [b.points]), idx(new unsigned int [b.points]), 
          r_idx(new unsigned int [b.points]), points(b.points)
        {
            memcpy(&this->l_wgt[0], &b.l_wgt[0], (points * sizeof(T)));
            memcpy(&this->r_wgt[0], &b.r_wgt[0], (points * sizeof(T)));
            memcpy(&this->sigma[0], &b.sigma[0], (points * sizeof(T)));
            memcpy(&this->l_idx[0], &b.l_idx[0], (points * sizeof(unsigned int)));
            memcpy(&this->idx[0],   &b.idx[0],   (points * sizeof(unsigned int)));
            memcpy(&this->r_idx[0], &b.r_idx[0], (points * sizeof(unsigned int)));
        }
        
        /* DTOR */
        ~brownian_bridge()
        {
            delete [] l_wgt;
            delete [] r_wgt;
            delete [] sigma;
            delete [] tmp;
            delete [] l_idx;
            delete [] idx;
            delete [] r_idx;
        }
        
        
        void map_randoms(const T *const rands, T *const mapped_rands) const
        {
            tmp[points - 1] = sigma[0] * rands[0];
            for (unsigned int i = 1; i < points; i++)
            {
                const unsigned int j = l_idx[i];
                const unsigned int l = idx[i];
                const unsigned int k = r_idx[i];
                if (j > 0)
                {
                    tmp[l] = l_wgt[i] * tmp[j - 1] + r_wgt[i] * tmp[k] + sigma[i] * rands[i];
                }
                else
                {
                    tmp[l] = r_wgt[i] * tmp[k] + sigma[i] * rands[i];
                }
            }
            
            std::adjacent_difference(&tmp[0], &tmp[points - 1], &mapped_rands[0]);
        }
        
    private :
        brownian_bridge& operator=(const brownian_bridge &)
        {
            return *this;
        }
        
        T               *const  l_wgt;
        T               *const  r_wgt;
        T               *const  sigma;
        T               *const  tmp;
        unsigned int    *const  l_idx;
        unsigned int    *const  idx;
        unsigned int    *const  r_idx;
        const unsigned int      points;
};


#endif
