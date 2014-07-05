#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <assert.h>

#include "quasi_newton_solver.h"
#include "option.h"
#include "fd_solver.h"

#include "matrix.h"


bool quasi_newton_solver::line_search(std::vector<option*> &xold, const double fold, double *g, double *p,
	std::vector<option*> &x, double &f, const double stpmax, int n, double s)
{
    /* Find the magnitude of the Newton direction */
	double sum = 0.0;
	for (int i = 0; i < n; i++)
    {
        sum += p[i] * p[i];
    }
	sum = sqrt(sum);
    
    /* If the magnitude is too long scale the Newton direction */
	if (sum > stpmax)
    {
		for (int i = 0; i < n; i++)
        {
            p[i] *= stpmax / sum;
        }
    }
    
	double slope = 0.0;
	for (int i = 0; i < n; i++)
    {
		slope += g[i] * p[i];
	}
//    std::cout << "slope: " << slope << std::endl;
    
    /* Algorithm must be moving up hill */
    assert(slope >= 0.0);
    
    /* Find the minimum scale */
	double test = 0.0;
	for (int i = 0; i < n; i++)
    {
		const double tmp = fabs(p[i]) / std::max(fabs(xold[i]->get_notional()), 1.0);
        test = std::max(tmp, test);
	}
    const double alamin = std::numeric_limits<double>::epsilon() / test;
	
    /* Start by trying to move the full distance in the Newton direction */
    /* This gives quadratic convergence if it works */
    double alam     = 1.0;
	double alam2    = 0.0;
    double f2       = 0.0;
	while (true)
    {
//        std::cout << "Line search with lambda " << alam << std::endl;
        /* Update by moving in the Newton direction */
		for (int i = 0; i < n; i++)
        {
            x[i]->set_notional(xold[i]->get_notional() + alam * p[i]);
        }
        
		f = solver.solve(x, s);
//        std::cout << "target f: " << (fold + alpha * alam * slope) << std::endl;
//        std::cout << "f this iteration: " << f << std::endl;
        
        /* The set size has become too small, the algorithm has converged */
		double tmplam;
        if (alam < alamin)
        {
			for (int i = 0; i < n; i++)
            {
                x[i]->set_notional(xold[i]->get_notional());
            }
            
//            std::cout << "Convergence" << std::endl;
			return true;
		}
        /* Function has increased enough */
        else if (f >= (fold + (alpha * alam * slope)))
        {
//            std::cout << "moved far enough" << std::endl;
            return false;
        }
		else
        {
            /* First iteration so we can only use a quadratic to make the next guess */
			if (alam == 1.0)
            {
				tmplam = -slope / (2.0 * (f - fold - slope));
			}
            /* Use cubic approximation to make the next guess */
            else
            {
				const double rhs1 = f - fold - alam * slope;
				const double rhs2 = f2 - fold - alam2 * slope;
				const double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
				const double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
				if (a == 0.0)
                {
                    tmplam = -slope / (2.0 * b);
                }
				else
                {
					const double disc = b * b - 3.0 * a * slope;
					if (disc < 0.0)
                    {
                        tmplam = 0.5 * alam;
                    }
					else if (b <= 0.0)
                    {
                        tmplam = (-b + sqrt(disc)) / (3.0 * a);
                    }
					else
                    {
                        tmplam = -slope / (b + sqrt(disc));
                    }
				}
			}
		}
		alam2 = alam;
		f2 = f;
        
        /* Guard against too faster decrease of the step size */
		alam = std::min(0.5 * alam, std::max(tmplam, 0.1 * alam));
	}
}


/* Uses finite differencing to calculate the partial derivative in all dimensions */
void quasi_newton_solver::direction(const std::vector<option*> &p, double *df, const double f, const double s, const int n)
{
    /* In each free direction */
    for (int i = 0; i < n; i++)
    {
        /* Move a small positive distance */
        const double cur_not = p[i]->get_notional();
        double h = finite_diff_step * fabs(cur_not);
        if (h == 0.0)
        {
            h = finite_diff_step;
        }
        p[i]->set_notional(cur_not + h);

        /* Use finite difference to calculate partial derivative */
        /* Use the actual movement to reduce error */
        const double fh = solver.solve(p, s);
        df[i] = (fh - f) / (p[i]->get_notional() - cur_not);
        p[i]->set_notional(cur_not);
    }
}


double quasi_newton_solver::maximise(const std::vector<option*> &hedge, option &target, const double s)
{
    std::vector<option*> p(hedge.begin(), hedge.end());
    p.push_back(&target);

	double *dg      = new double [hedge.size()];
    double *g       = new double [hedge.size()];
    double *hdg     = new double [hedge.size()];
    double *xi      = new double [hedge.size()];
	double *hessin  = new double [hedge.size() * hedge.size()];
    memset(&hessin[0], 0, (hedge.size() * hedge.size() * sizeof(double)));

    std::vector<option*> pnew(p.size());
    for (unsigned int i = 0; i < p.size(); i++)
    {
        pnew[i] = p[i]->clone();
    }
    
    /* Initialise the hessin to the unit matrix and get the initial point and gradient */
	double fp = solver.solve(p, s);
	direction(p, g, fp, s, hedge.size());
    double sum = 0.0;
	for (unsigned int i = 0; i < hedge.size(); i++)
    {
        hessin[(i * hedge.size()) + i] = 1.0;
		xi[i]               = g[i];
		sum                += p[i]->get_notional() * p[i]->get_notional();
	}
	double stpmax = max_step * std::max(sqrt(sum), static_cast<double>(hedge.size()));

    /* Iterate up the gradient */
    int iter = 0;
    double fret;
    double x_mov;
    double g_siz;
	do
    {
        /* Check the iteration limit */
        assert(iter++ <= max_iter);
        
        /* Maximise along the current gradient */
		line_search(p, fp, g, xi, pnew, fret, stpmax, hedge.size(), s);
		fp = fret;
		for (unsigned int i = 0; i < hedge.size(); i++)
        {
			xi[i] = pnew[i]->get_notional() - p[i]->get_notional();
			p[i]->set_notional(pnew[i]->get_notional());
//            std::cout << "new notional: " << p[i]->get_notional() << std::endl;
		}
        
        /* Check for convergence in the movement of the point */
		x_mov = 0.0;
		for (unsigned int i = 0; i < hedge.size(); i++)
        {
			const double tmp = fabs(xi[i]) / std::max(fabs(p[i]->get_notional()), 1.0);
            x_mov = std::max(tmp, x_mov);
		}
        
		
        /* Test for gradient is zero */
        for (unsigned int i = 0; i < hedge.size(); i++)
        {
            dg[i] = g[i];
        }
        
        direction(p, g, fp, s, hedge.size());
        g_siz = 0.0;
		const double den = std::max(fret, 1.0);
		for (unsigned int i = 0; i < hedge.size(); i++)
        {
            const double tmp = fabs(g[i]) * std::max(fabs(p[i]->get_notional()), 1.0) / den;
            g_siz = std::max(tmp, g_siz);
        }
        
        /* Difference in the new and old gradient */
        for (unsigned int i = 0; i < hedge.size(); i++)
        {
            dg[i] = g[i] - dg[i];
        }
        
        for (unsigned int i = 0; i < hedge.size(); i++)
        {
            hdg[i] = 0.0;
            for (unsigned int j = 0; j < hedge.size(); j++)
            {
                hdg[i] += hessin[(i * hedge.size()) + j] * dg[j];
            }
        }
        
        /* Dot product of denominators */
        double fac      = 0.0;
        double fae      = 0.0; 
        double sumdg    = 0.0;
        double sumxi    = 0.0;
        for (unsigned int i = 0; i < hedge.size(); i++)
        {
            fac += dg[i] * xi[i];
            fae += dg[i] * hdg[i];
            sumdg += dg[i] * dg[i];
            sumxi += xi[i] * xi[i];
        }
        
        /* Skip BFGS update if not big enough */
        if (fac > std::sqrt(std::numeric_limits<double>::epsilon() * sumdg * sumxi))
        {
            fac = 1.0 / fac;
            const double fad = 1.0 / fae;
            for (unsigned int i = 0; i < hedge.size(); i++)
            {
                dg[i] = fac * xi[i] - fad * hdg[i];
            }
            
            /* Update the hessin */
            for (unsigned int i = 0; i < hedge.size(); i++)
            {
                for (unsigned int j = i; j < hedge.size(); j++)
                {
                    hessin[(i * hedge.size()) + j] += fac * xi[i] * xi[j]
                        - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
                    hessin[(j * hedge.size()) + i] = hessin[(i * hedge.size() ) + j];
                }
            }
        }
        
        /* Next conjugate direction to go */
        for (unsigned int i = 0; i < hedge.size(); i++)
        {
			xi[i] = 0.0;
			for (unsigned int j = 0; j < hedge.size(); j++)
            {
                xi[i] += hessin[(i * hedge.size()) + j] * g[j];
            }
		}
	} while ((g_siz >= grad_tol) && (x_mov >= tolerance));
    
    /* Clean up */
    delete [] hessin;
	delete [] dg;
    delete [] g;
    delete [] hdg;
    delete [] xi;
    
    std::cout << fret;
    for (unsigned int i = 0; i < pnew.size(); i++)
    {
        std::cout << " " << pnew[i]->get_notional();
        delete pnew[i];
    }
    std::cout << std::endl;
    
    return fret;
}
