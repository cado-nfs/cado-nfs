#ifndef UTILS_RANDOM_DISTRIBUTIONS_HPP_
#define UTILS_RANDOM_DISTRIBUTIONS_HPP_

#include <gmp.h>

double random_uniform(gmp_randstate_t rstate);

double random_normal_standard(gmp_randstate_t rstate);

double random_normal(gmp_randstate_t rstate, double mean, double sdev);

/* return the expected maximum among n picks for a normal low with given
 * mean and sdev */
double extreme_normal(double n, double mean, double sdev);

double random_normal_constrained(gmp_randstate_t rstate, double mean, double sdev, double a, double b);

/* given a probability mass function which gives the gaussian with mean
 * and sdev given my mx[0] and mx[1], but truncated to the interval
 * [a,b[, return the mean and sdev of the resulting distribution.
 *
 * This is just an illustration, which can be used to witness how the
 * normal approximation can end up being catastrophic if we're more
 * Poisson-like.
 */
void accuracy_of_normal_approximation_to_binomial(double * my, double *mx, unsigned long a, unsigned long b);

double random_poisson(gmp_randstate_t rstate, double lambda);

/* This is the random variable associated to the *size* of the sample */
double random_binomial(gmp_randstate_t rstate, unsigned long n, double p);


#endif	/* UTILS_RANDOM_DISTRIBUTIONS_HPP_ */
