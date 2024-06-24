#ifndef UTILS_RANDOM_DISTRIBUTIONS_HPP_
#define UTILS_RANDOM_DISTRIBUTIONS_HPP_

#include <gmp.h>

/* Some generic functions for random picking along distributions
 */

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


/* Let X be an integer random variable on [0,N) with probability mass
 * function p(i) and cumulative distribution function c(i), so that
 * Pr(X=k) = \int_{i}^{i+1} p(t) dt = c(i+1) - c(i)
 *
 * Let r(x) (from [0,1) to [0,N)) be the reciprocal of the cumulative
 * distribution function c.
 *
 * We want to pick a certain number of samples of this random variable,
 * _without repetition_.  To do so, the algorithm that we use works as
 * follows:
 *  - pick a uniform random integer x in [0,A). Initially, think of A=1.
 *  - compute i as the floor of r(x). The sample is i.
 *  - mark the interval c(i), c(i+1) as "forbidden": further picks will
 *  avoid it. This means that for the next pick, A is decreased by
 *  c(i+1)-c(i), and we need to adjust x according to what's on the left
 *  and on the right.
 *
 * The data structure that we use to run this algorithm is
 * punched_interval.
 *
 * A punched interval has an upper and a lower bound, and records the
 * number of holes it contains. If it has holes, it is subdivided in two
 * punched interval, on each side of one hole (and these intervals may,
 * in turn, contain further holes).
 */
struct punched_interval_s {
    double b0, b1;
    double holes;
    int has_left;
    /* free blocks use the "left" pointer below for the next argument in
     * the free list */
    struct punched_interval_s * left;
    struct punched_interval_s * right;
};
typedef struct punched_interval_s * punched_interval_ptr;

void punched_interval_free(punched_interval_ptr c, punched_interval_ptr * pool);

void punched_interval_set_full(punched_interval_ptr x, double b0, double b1);

punched_interval_ptr punched_interval_alloc(punched_interval_ptr * pool, double b0, double b1);

void punched_interval_free_pool(punched_interval_ptr * pool);

void punched_interval_pre_free_pool(punched_interval_ptr * pool, int max, int print);

void punched_interval_print_rec(FILE * f, punched_interval_ptr c);

void punched_interval_print(FILE * f, punched_interval_ptr c);

unsigned long punched_interval_pick(punched_interval_ptr * pool, punched_interval_ptr c,
        double (*dist_q)(const void *, double), 
        double (*dist_qrev)(const void *, double), 
        const void * f,
        gmp_randstate_ptr rstate);


#endif	/* UTILS_RANDOM_DISTRIBUTIONS_HPP_ */
