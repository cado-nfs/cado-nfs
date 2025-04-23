#ifndef CADO_UTILS_RANDOM_DISTRIBUTIONS_HPP
#define CADO_UTILS_RANDOM_DISTRIBUTIONS_HPP

#include <cstdio>

#include <memory>

#include "gmp_aux.h"

/* Some generic functions for random picking along distributions
 */

double random_uniform(gmp_randstate_t rstate);

double random_normal_standard(cxx_gmp_randstate & rstate);

double random_normal(cxx_gmp_randstate & rstate, double mean, double sdev);

/* return the expected maximum among n picks for a normal low with given
 * mean and sdev */
double extreme_normal(double n, double mean, double sdev);

double random_normal_constrained(cxx_gmp_randstate & rstate, double mean, double sdev, double a, double b);

/* given a probability mass function which gives the gaussian with mean
 * and sdev given my mx[0] and mx[1], but truncated to the interval
 * [a,b[, return the mean and sdev of the resulting distribution.
 *
 * This is just an illustration, which can be used to witness how the
 * normal approximation can end up being catastrophic if we're more
 * Poisson-like.
 */
void accuracy_of_normal_approximation_to_binomial(double * my, const double *mx, unsigned long a, unsigned long b);

double random_poisson(cxx_gmp_randstate & rstate, double lambda);

/* This is the random variable associated to the *size* of the sample */
double random_binomial(cxx_gmp_randstate & rstate, unsigned long n, double p);

struct matrix_column_distribution {
    // virtual double p(double x) const = 0;
    virtual double q(double x) const = 0;
    virtual double qrev(double x) const = 0;
    // virtual double qq(double x) const = 0;
    matrix_column_distribution() = default;
    matrix_column_distribution(matrix_column_distribution const&) = default;
    matrix_column_distribution(matrix_column_distribution &&) = default;
    matrix_column_distribution& operator=(matrix_column_distribution const&) = default;
    matrix_column_distribution& operator=(matrix_column_distribution &&) = default;
    virtual ~matrix_column_distribution() = default;
};


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
struct punched_interval {
    double b0 = 0, b1 = 0;
    double holes = 0;
    int has_left = 0;

    typedef std::unique_ptr<punched_interval> node_t;
    std::unique_ptr<punched_interval> left, right;

    /* free blocks use the "left" pointer below for the next argument in
     * the free list */
    typedef node_t pool_t;

    static void recycle(node_t &&, pool_t &);
    static node_t alloc(pool_t &, double b0, double b1);

    void print_rec(FILE *) const;
    void print(FILE *) const;

    unsigned long pick(pool_t & pool,
        matrix_column_distribution const & D,
        cxx_gmp_randstate & rstate);

    private:
    void punch_inner(pool_t & pool, double x0, double x1);
    unsigned long pick_inner(pool_t & pool,
        matrix_column_distribution const & D,
        double x);
};

#endif	/* CADO_UTILS_RANDOM_DISTRIBUTIONS_HPP */
