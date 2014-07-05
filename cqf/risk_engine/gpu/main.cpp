#include <iostream>
#include <math.h>
#include <memory>

#include "boost/chrono.hpp"

#include "cyclic_reduction.cuh"
#include "volatility_inversion.cuh"


void volatility_inversion_host(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num)
{
	for (int option = 0; option < num; ++option)
	{
		float error;
		int i = 0;
		do
		{
			const float price = call_price(s, r, v[option], t[option], k[option]);
			const float vega = call_vega(s, r, v[option], t[option], k[option]);
		
			error = (p[option] - price);
			v[option] += error / vega;
		} while ((fabs(error) > tol) && (i++ < iter));
	}
}


template<int GUESSES>
void parallel_volatility_inversion_host(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num)
{
	/* Build guesses */
	float v_ladder[GUESSES];
	float err_ladder[GUESSES];
	float ferr_ladder[GUESSES];

	for (int option = 0; option < num; ++option)
	{
		std::cout << std::endl << "Option: " << option << std::endl;
		float ladder_span = 0.04f;	/* Span the guesses over 4% */
		float v_min = v[option] - (ladder_span * 0.5f);
		float v_max = v[option] + (ladder_span * 0.5f);	

		int i = 0;
		do
		{
			/* Do some guessing */
			for (int guess = 0; guess < GUESSES; ++guess)
			{
				/* Work out guess */
				v_ladder[guess] = v_min + ((v_max - v_min) * (guess / static_cast<float>(GUESSES - 1)));

				/* Price */
				err_ladder[guess] = call_price(s, r, v_ladder[guess], t[option], k[option]) - p[option];
				ferr_ladder[guess] = fabs(err_ladder[guess]);
				std::cout << "Guess: " << v_ladder[guess] << " at: " << guess << " gave error: " << err_ladder[guess] << std::endl;
			}
		
			/* Find minimum error */
			int min_err_pos = 0;
			for (int guess = 1; guess < GUESSES; ++guess)
			{
				if (ferr_ladder[guess] < ferr_ladder[min_err_pos])
				{
					min_err_pos = guess;
				}
			}

			if (ferr_ladder[min_err_pos] < tol)
			{
				v[option] = v_ladder[min_err_pos];
				break;
			}

			/* Pick the span for the next ladder */
			/* Doesnt matter if v_min is actually higher than v_max so long as 0 is crossed */
			if ((err_ladder[0] * err_ladder[GUESSES - 1]) >= 0.0f) /* Root not found (or two roots found) */
			{
				ladder_span *= 2.0f;
				if (abs(err_ladder[0] - err_ladder[GUESSES - 1]) < tol) /* Ladder is very flat so no direction */
				{
					v_min -= ladder_span * 0.5f;
					v_max += ladder_span * 0.5f;
					ladder_span *= 2.0f;
				}
				else if (ferr_ladder[0] < ferr_ladder[GUESSES - 1]) /* Lower end is closer to root so expand it */
				{
					v_min = v_ladder[0];
					v_max = v_min - ladder_span;
				}
				else /* Upper end is closer to root so expand it */
				{
					v_min = v_ladder[GUESSES - 1];
					v_max = v_min + ladder_span;
				}
			}
			else if (min_err_pos == 0) /* Root found at lower extreme of ladder */
			{
				v_min = v_ladder[min_err_pos];
				if ((err_ladder[0] * err_ladder[1]) >= 0.0f)
				{
					ladder_span *= 2.0f;
					v_max = (v_min - ladder_span);
				}
				else
				{
					ladder_span *= 0.5f;
					v_max = v_ladder[1];
				}
			}
			else if (min_err_pos == (GUESSES - 1)) /* Root found at upper extreme of ladder */
			{
				v_min = v_ladder[min_err_pos];
				if ((err_ladder[GUESSES - 1] * err_ladder[GUESSES - 2]) >= 0.0f)
				{
					ladder_span *= 0.5f;
					v_max = v_ladder[GUESSES - 2];
				}
				else
				{
					ladder_span *= 2.0f;
					v_max = (v_min + ladder_span);
				}
			}
			else if ((err_ladder[min_err_pos] * err_ladder[min_err_pos - 1]) < 0.0f) /* Root in bin below min error */
			{
				ladder_span *= 1.0f / GUESSES;
				v_min = v_ladder[min_err_pos];
				v_max = v_ladder[min_err_pos - 1];
			}
			else /* Root in bin above min error */
			{
				ladder_span *= 1.0f / GUESSES;
				v_min = v_ladder[min_err_pos];
				v_max = v_ladder[min_err_pos + 1];
			}
			std::cout << std::endl;
		} while (i++ < iter);
	}
}


void vega_guided_parallel_volatility_inversion_host(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num, const int guesses)
{
	/* Build guesses */
	float v_ladder[32];
	float vega_ladder[32];
	float err_ladder[32];
	float ferr_ladder[32];

	float ladder_span = 0.04f;	/* Span the guesses over 4% */
	for (int option = 0; option < num; ++option)
	{
		float v_mid = v[option];

		int i = 0;
		int min_err_pos;
		do
		{
			/* Do some guessing */
			for (int guess = 0; guess < guesses; ++guess)
			{
				/* Work out guess */
				v_ladder[guess] = v_mid + (ladder_span * (guess / static_cast<float>(guesses - 1))) - (ladder_span * 0.5f);

				/* Price */
				err_ladder[guess] = p[option] - call_price(s, r, v_ladder[guess], t[option], k[option]);
				vega_ladder[guess] = call_vega(s, r, v_ladder[guess], t[option], k[option]);
				ferr_ladder[guess] = abs(err_ladder[guess]);
				std::cout << "Guess: " << v_ladder[guess] << " Error: " << err_ladder[guess] << " Vega: " << vega_ladder[guess] << std::endl;
			}
		
			/* Find minimum error */
			min_err_pos = 0;
			for (int guess = 1; guess < guesses; ++guess)
			{
				if (ferr_ladder[guess] < ferr_ladder[min_err_pos])
				{
					min_err_pos = guess;
				}
			}

			/* Pick the span for the next ladder */
			v_mid = v_ladder[min_err_pos] + (err_ladder[min_err_pos] / vega_ladder[min_err_pos]);
			std::cout << "Min Err Pos: " << min_err_pos << " New Mid: " << v_mid << std::endl;
			ladder_span *= 0.5f;
		
		} while ((ferr_ladder[min_err_pos] > tol) && (i++ < iter));

		v[option] = v_ladder[min_err_pos];
	}
}


void prepare_data(float *const guess_near, float *const guess_far, const int num)
{
	for (int i = 0; i < num; ++i)
	{
		guess_far[i] = 0.211f;
		guess_near[i] = 0.291f;
	}
}


void check_results(const float *const v_exp, float *const v_guess, const int num, const float tol)
{
	/* Check result */
	for (int i = 0; i < num; ++i)
	{
		if (fabs(v_guess[i] - v_exp[i]) > tol)
		{
			std::cout << "Error at option: " << i << ", Got: " << v_guess[i] << ", Expected: " << v_exp[i] << std::endl;
		}
	}

	/* Clear guess */
	for (int i = 0; i < num; ++i)
	{
		v_guess[i] = 0.20f;
	}
}


void volatility_inversion_test()
{
	/* Build some test instruments */
	const int num_t = 12;
	const int num_k = 4;
	const int num = num_k * num_t;

	const float s = 100.0f;
	const float r = 0.05f;

	const int iter = 20;
	const float tol = 0.00005f;
	std::unique_ptr<float []> guess_near(new float [num]);
	std::unique_ptr<float []> guess_far(new float [num]);

	const float t_base = 0.5f;
	const float t_period = 0.5f;
	std::unique_ptr<float []> t(new float [num]);
	for (int i = 0; i < num_t; ++i)
	{
		for (int j = 0; j < num_k; ++j)
		{
			t[(i * num_k) + j] = t_base + i * t_period;
		}
	}

	const float k_base = 0.5f;
	const float k_period = 0.1f;
	std::unique_ptr<float []> k(new float [num]);
	for (int i = 0; i < num_t; ++i)
	{
		for (int j = 0; j < num_k; ++j)
		{
			k[(i * num_k) + j] = (k_base + j * k_period) * s;
		}
	}

	std::unique_ptr<float []> v_exp(new float [num]);
	for (int i = 0; i < num; ++i)
	{
		v_exp[i] = 0.3f;
	}

	std::unique_ptr<float []> p(new float [num]);
	for (int i = 0; i < num; ++i)
	{
		p[i] = call_price(s, r, v_exp[i], t[i], k[i]);
	}

	/* Time on host */
	prepare_data(guess_near.get(), guess_far.get(), num);
	const boost::chrono::high_resolution_clock::time_point host_t0 = boost::chrono::high_resolution_clock::now();
	volatility_inversion_host(s, r, guess_near.get(), t.get(), k.get(), p.get(), tol, iter, num);
	const boost::chrono::high_resolution_clock::time_point host_t1 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_near.get(), num, tol);
	std::cout << "    Near took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t1 - host_t0) << std::endl;

	const boost::chrono::high_resolution_clock::time_point host_t2 = boost::chrono::high_resolution_clock::now();
	volatility_inversion_host(s, r, guess_far.get(), t.get(), k.get(), p.get(), tol, iter, num);
	const boost::chrono::high_resolution_clock::time_point host_t3 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_far.get(), num, tol);
	std::cout << "    Far took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t3 - host_t2) << std::endl;

	/* Time on device */
	std::cout << "Serial volatility inversion" << std::endl;
	prepare_data(guess_near.get(), guess_far.get(), num);
	const boost::chrono::high_resolution_clock::time_point device_t0 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_near.get(), t.get(), k.get(), p.get(), tol, iter, num, serial);
	const boost::chrono::high_resolution_clock::time_point device_t1 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_near.get(), num, tol);
	std::cout << "    Near took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t1 - device_t0) << std::endl;

	const boost::chrono::high_resolution_clock::time_point device_t2 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_far.get(), t.get(), k.get(), p.get(), tol, iter, num, serial);
	const boost::chrono::high_resolution_clock::time_point device_t3 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_far.get(), num, tol);
	std::cout << "    Far took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t3 - device_t2) << std::endl;

	std::cout << "Parallel volatility inversion" << std::endl;
	prepare_data(guess_near.get(), guess_far.get(), num);
	const boost::chrono::high_resolution_clock::time_point device_t4 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_near.get(), t.get(), k.get(), p.get(), tol, iter, num, parallel);
	//parallel_volatility_inversion_host<4>(s, r, guess_near.get(), t.get(), k.get(), p.get(), tol, iter, num);
	const boost::chrono::high_resolution_clock::time_point device_t5 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_near.get(), num, tol);
	std::cout << "    Near took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t5 - device_t4) << std::endl;

	const boost::chrono::high_resolution_clock::time_point device_t6 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_far.get(), t.get(), k.get(), p.get(), tol, iter, num, parallel);
	//parallel_volatility_inversion_host<4>(s, r, guess_far.get(), t.get(), k.get(), p.get(), tol, iter, num);
	const boost::chrono::high_resolution_clock::time_point device_t7 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_far.get(), num, tol);
	std::cout << "    Far took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t7 - device_t6) << std::endl;

	std::cout << "Vega guided parallel volatility inversion" << std::endl;
	prepare_data(guess_near.get(), guess_far.get(), num);
	const boost::chrono::high_resolution_clock::time_point device_t8 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_near.get(), t.get(), k.get(), p.get(), tol, iter, num, vega_parallel);
	const boost::chrono::high_resolution_clock::time_point device_t9 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_near.get(), num, tol);
	std::cout << "    Near took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t9 - device_t8) << std::endl;

	const boost::chrono::high_resolution_clock::time_point device_t10 = boost::chrono::high_resolution_clock::now();
	volatility_inversion(s, r, guess_far.get(), t.get(), k.get(), p.get(), tol, iter, num, vega_parallel);
	const boost::chrono::high_resolution_clock::time_point device_t11 = boost::chrono::high_resolution_clock::now();
	check_results(v_exp.get(), guess_far.get(), num, tol);
	std::cout << "    Far took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t11 - device_t10) << std::endl;
}


void parallel_cyclic_reduction_host(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	std::unique_ptr<float []> lower_tmp(new float [dim]);
	std::unique_ptr<float []> upper_tmp(new float [dim]);
	std::unique_ptr<float []> result_tmp(new float [dim]);
	
	for (int span = 1 ; span < dim; span *= 2)
	{
		for (int rank = 0; rank < dim; ++rank)
		{
			if (rank - span >= 0)
			{
				lower_tmp[rank] = (diagonal[rank - span] != 0.0f) ? (-lower[rank] / diagonal[rank - span]) : 0.0f;
			}
			else
			{
				lower_tmp[rank] = 0.0f;
			}

			if (rank + span < dim)
			{
				upper_tmp[rank] = (diagonal[rank + span] != 0.0f) ? (-upper[rank] / diagonal[rank + span]) : 0.0f;
			}
			else
			{
				upper_tmp[rank] = 0.0f;
			}

			result_tmp[rank] = equal[rank];
		}

		for (int rank = 0; rank < dim; ++rank)
		{
			if (rank - span >= 0)
			{
				diagonal[rank] += lower_tmp[rank] * upper[rank - span];
				result_tmp[rank] += lower_tmp[rank] * equal[rank - span];
				lower_tmp[rank] *= lower[rank - span];
			}

			if (rank + span < dim)
			{
				diagonal[rank] += upper_tmp[rank] * lower[rank + span];
				result_tmp[rank] += upper_tmp[rank] * equal[rank + span];
				upper_tmp[rank] *= upper[rank + span];
			}
		}

		for (int rank = 0; rank < dim; ++rank)
		{
			lower[rank] = lower_tmp[rank];
			upper[rank] = upper_tmp[rank];
			equal[rank] = result_tmp[rank];
		}
	}

	for (int rank = 0; rank < dim; ++rank)
	{
		equal[rank] /= diagonal[rank];
	}
}


int cyclic_reduction_forward_reduction_host(float *lower, float *diagonal, float *upper, float *equal, const int dim, const int from, const int to)
{
	/* Forward reduction */
	int step = from;
	for (; (step * to * 3) <= dim; step <<= 1)
	{
		for (int addr = (step << 1) - 1; addr < dim; addr += (step << 1))
		{
			if (addr - step >= 0)
			{
				const float alpha	= -lower[addr] / diagonal[addr - step];
				equal[addr]    += (alpha * equal[addr - step]);
				diagonal[addr] += (alpha * upper[addr - step]);
				lower[addr]		= alpha * lower[addr - step];
			}

			if (addr + step < dim)
			{
				const float gamma	= -upper[addr] / diagonal[addr + step];
				equal[addr]	   += (gamma * equal[addr + step]);
				diagonal[addr] += (gamma * lower[addr + step]);
				upper[addr]		= gamma * upper[addr + step];
			}
		}
	}

	return step;
}


void cyclic_reduction_back_substitution_host(float *lower, float *diagonal, float *upper, float *equal, const int dim, const int from, const int to)
{
	/* Backward Substitution */
	int step = from;
	for (; step > 0; step >>= 1)
	{
		for (int addr = step - 1; addr < dim; addr += (step << 1))
		{
			if (addr - step >= 0)
			{
				equal[addr] -= (lower[addr] * equal[addr - step]);
			}

			if (addr + step < dim)
			{
				equal[addr] -= (upper[addr] * equal[addr + step]);
			}

			equal[addr] = equal[addr] / diagonal[addr];
		}
	}
}


void cyclic_reduction_host(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	/* Forward reduction */
	int step = cyclic_reduction_forward_reduction_host(lower, diagonal, upper, equal, dim, 1, 1);

	/* Solve base system */
	if ((dim / step) == 2) /* Solve simultaneous equations */
	{
		const int equal_addr = (step << 1) - 1;
		const float a0 = diagonal[equal_addr - step];
		const float a1 = lower[equal_addr];
		const float b0 = upper[equal_addr - step];
		const float b1 = diagonal[equal_addr];
		const float c0 = equal[equal_addr - step];
		const float c1 = equal[equal_addr];

		equal[equal_addr] = (c0 * a1 - a0 * c1) / (a1 * b0 - a0 * b1);
		equal[equal_addr - step] = (c0 - b0 * equal[equal_addr]) / a0;
	}
	else /* blk_size == 1, equations are already solved */
	{
		const int equal_addr = step - 1;
		equal[equal_addr] = equal[equal_addr] / diagonal[equal_addr];
	}
	step >>= 1;

	/* Backward Substitution */
	cyclic_reduction_back_substitution_host(lower, diagonal, upper, equal, dim, step, 0);
}


void hybrid_cyclic_reduction_host(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	/* Cyclic forward reduction */
	int step = cyclic_reduction_forward_reduction_host(lower, diagonal, upper, equal, dim, 1, 4);

	/* Parallel cyclic reduction to solve system */
	std::unique_ptr<float []> lower_tmp(new float [dim]);
	std::unique_ptr<float []> upper_tmp(new float [dim]);
	std::unique_ptr<float []> result_tmp(new float [dim]);
	for (int span = step; span < dim; span *= 2)
	{
		for (int rank = step - 1; rank < dim; rank += step)
		{
			if (rank - span >= 0)
			{
				lower_tmp[rank] = (diagonal[rank - span] != 0.0f) ? (-lower[rank] / diagonal[rank - span]) : 0.0f;
			}
			else
			{
				lower_tmp[rank] = 0.0f;
			}

			if (rank + span < dim)
			{
				upper_tmp[rank] = (diagonal[rank + span] != 0.0f) ? (-upper[rank] / diagonal[rank + span]) : 0.0f;
			}
			else
			{
				upper_tmp[rank] = 0.0f;
			}

			result_tmp[rank] = equal[rank];
		}		

		for (int rank = step - 1; rank < dim; rank += step)
		{
			if (rank - span >= 0)
			{
				diagonal[rank] += lower_tmp[rank] * upper[rank - span];
				result_tmp[rank] += lower_tmp[rank] * equal[rank - span];
				lower_tmp[rank] *= lower[rank - span];
			}

			if (rank + span < dim)
			{
				diagonal[rank] += upper_tmp[rank] * lower[rank + span];
				result_tmp[rank] += upper_tmp[rank] * equal[rank + span];
				upper_tmp[rank] *= upper[rank + span];
			}
		}

		for (int rank = step - 1; rank < dim; rank += step)
		{
			lower[rank] = lower_tmp[rank];
			upper[rank] = upper_tmp[rank];
			equal[rank] = result_tmp[rank];
		}
	}

	for (int rank = step - 1; rank < dim; rank += step)
	{
		equal[rank] /= diagonal[rank];
	}

	/* Cyclic backward substitution */
	cyclic_reduction_back_substitution_host(lower, diagonal, upper, equal, dim, step >> 1, 0);
}


void gaussian_elimination_host(float *lower, float *diagonal, float *upper, float *equal, const int dim)
{
	/* Forward reduction */
	for (int i = 1; i < dim; ++i)
	{
		const float alpha = lower[i] / diagonal[i - 1];
		diagonal[i] -= alpha * upper[i - 1];
		equal[i] -= alpha * equal[i - 1];
	}

	/* Backward substitution */
	equal[dim - 1] /= diagonal[dim - 1];
	for (int i = dim - 2; i >= 0; --i)
	{
		equal[i] -= upper[i] * equal[i + 1];
		equal[i] /= diagonal[i];
	}
}


void prepare_data(float *const lower, float *const diagonal, float *const upper, float *const equal, const float *const expected, const int dim)
{
	/* Build a test matrix */
	const float lower_base = 1.0f;
	const float lower_step = 1.0f;
	//std::cout << "Lower: ";
	for (int i = 0; i < dim; ++i)
	{
		lower[i] = std::sin(lower_base + (i * lower_step));
		//std::cout << lower[i] << ", ";
	}
	//std::cout << std::endl;

	const float diagonal_base = 1.0f;
	const float diagonal_step = 2.5f;
	//std::cout << "Diagonal: ";
	for (int i = 0; i < dim; ++i)
	{
		diagonal[i] = std::cos(diagonal_base + (i * diagonal_step));
		//std::cout << diagonal[i] << ", ";
	}
	//std::cout << std::endl;

	const float upper_base = 1.0f;
	const float upper_step = 1.5f;
	//std::cout << "Upper: ";
	for (int i = 0; i < dim; ++i)
	{
		upper[i] = std::sin(upper_base + (i * upper_step));
		//std::cout << upper[i] << ", ";
	}
	//std::cout << std::endl;

	//std::cout << "Equal: ";
	for (int i = 0; i < dim; ++i)
	{
		if (i == 0)
		{
			equal[i] = (expected[i] * diagonal[i]) + (expected[i + 1] * upper[i]);
		}
		else if (i == (dim - 1))
		{
			equal[i] = (expected[i - 1] * lower[i]) + (expected[i] * diagonal[i]);
		}
		else
		{
			equal[i] = (expected[i - 1] * lower[i]) + (expected[i] * diagonal[i]) + (expected[i + 1] * upper[i]);
		}
		//std::cout << equal[i] << ", ";
	}
	//std::cout << std::endl << std::endl;
}


void check_tridiagonal_results(const float *const equal, const float *const expected, const int dim)
{
	/* Check result */
	for (int i = 0; i < dim; ++i)
	{
		if (fabs(equal[i] - expected[i]) > 0.1f)
		{
			std::cout << "Error at row: " << i << ", Got: " << equal[i] << ", Expected: " << expected[i] << std::endl;
		}
	}
}


void tridiagonal_test(const int dim)
{
	const float expected_base = 1.0f;
	const float expected_step = 1.0f;
	std::unique_ptr<float []> expected(new float [dim]);
	for (int i = 0; i < dim; ++i)
	{
		expected[i] = std::sin(expected_base + (i * expected_step));
	}

	std::unique_ptr<float []> lower(new float [dim]);
	std::unique_ptr<float []> diagonal(new float [dim]);
	std::unique_ptr<float []> upper(new float [dim]);
	std::unique_ptr<float []> equal(new float [dim]);

	/* Cyclic reduction */
	/* Time on host */
	std::cout << "Cyclic reduction: " << std::endl;
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t0 = boost::chrono::high_resolution_clock::now();
	cyclic_reduction_host(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t1 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Host took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t1 - host_t0) << std::endl;
	
	/* Time on device */
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t0 = boost::chrono::high_resolution_clock::now();
	cyclic_reduction(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t1 = boost::chrono::high_resolution_clock::now();
	std::cout << "    Device took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t1 - device_t0) << std::endl;
	
	/* Parallel cyclic reduction */
	/* Time on host */
	std::cout << "Parallel cyclic reduction: " << std::endl;
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t2 = boost::chrono::high_resolution_clock::now();
	parallel_cyclic_reduction_host(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t3 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Host took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t3 - host_t2) << std::endl;
	
	/* Time on device */
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t2 = boost::chrono::high_resolution_clock::now();
	parallel_cyclic_reduction(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t3 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Device took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t3 - device_t2) << std::endl;

	/* Hybrid cyclic reduction */
	/* Time on host */
	std::cout << "Hybrid cyclic reduction: " << std::endl;
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t4 = boost::chrono::high_resolution_clock::now();
	hybrid_cyclic_reduction_host(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t5 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Host took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t5 - host_t4) << std::endl;

	/* Time on device */
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t4 = boost::chrono::high_resolution_clock::now();
	hybrid_cyclic_reduction(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point device_t5 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Device took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(device_t5 - device_t4) << std::endl;
	
	/* Gaussian elimination */
	std::cout << "Gaussian elimination: " << std::endl;
	prepare_data(lower.get(), diagonal.get(), upper.get(), equal.get(), expected.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t6 = boost::chrono::high_resolution_clock::now();
	gaussian_elimination_host(lower.get(), diagonal.get(), upper.get(), equal.get(), dim);
	const boost::chrono::high_resolution_clock::time_point host_t7 = boost::chrono::high_resolution_clock::now();
	check_tridiagonal_results(equal.get(), expected.get(), dim);
	std::cout << "    Host took: " << boost::chrono::duration_cast<boost::chrono::microseconds>(host_t7 - host_t6) << std::endl;
}


int main()
{
	volatility_inversion_test();

//	for (int i = 3; i < 513; ++i)
//	{
//		std::cout << std::endl << "Running test: " << i << std::endl;
//		tridiagonal_test(i);
//	}
//	tridiagonal_test(512);

	return 0;
}