#ifndef CADO_BEST_POLYNOMIALS_QUEUE_H
#define CADO_BEST_POLYNOMIALS_QUEUE_H

#include <stddef.h>
#include "cado_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

struct best_polynomials_queue_s {
    /* actual implementation is hidden in C++ code */
    void * pimpl;
};
typedef struct best_polynomials_queue_s best_polynomials_queue[1];
typedef struct best_polynomials_queue_s * best_polynomials_queue_ptr;
typedef const struct best_polynomials_queue_s * best_polynomials_queue_srcptr;

void best_polynomials_queue_init(best_polynomials_queue_ptr b, int queue_lenth);
void best_polynomials_queue_clear(best_polynomials_queue_ptr b);
double best_polynomials_queue_get_best_score(best_polynomials_queue_srcptr b);
double best_polynomials_queue_get_worst_score(best_polynomials_queue_srcptr b);
size_t best_polynomials_queue_get_count(best_polynomials_queue_srcptr b);
void best_polynomials_queue_try_push(best_polynomials_queue_ptr b, cado_poly_srcptr cpoly, double score);
void best_polynomials_queue_print(best_polynomials_queue_srcptr b, FILE *, const char * prefix);
void best_polynomials_queue_do(best_polynomials_queue_srcptr b, void (*f)(int, double, cado_poly_ptr, void *), void *arg);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_BEST_POLYNOMIALS_QUEUE_H */
