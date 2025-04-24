#ifndef CADO_FM_HPP
#define CADO_FM_HPP

#include <cstdio>

typedef struct fm {
    // facul_method_code
    // ec_parameterization_t
    // unsigned long
    // unsigned long
    unsigned long * method; // contain: METHOD, CURVE, B1, B2
    double * proba;
    double * time;
    unsigned int len_method; // lenght of the method (default:4)
    unsigned int len_proba;  // index of array proba
    unsigned int len_time;   // index of array time
    unsigned int len_p_min;
    /*
       The prime number such that : proba[i] equals to the probability
       to find a prime number of len_p_min+i bits with our nmethod.
     */
} fm_t;

#ifdef __cplusplus
extern "C" {
#endif

fm_t * fm_create(void);

void fm_free(fm_t * t);

unsigned long const * fm_get_method(fm_t const * t);

double const * fm_get_proba(fm_t const * t);

double const * fm_get_time(fm_t const * t);

unsigned int fm_get_len_method(fm_t const * t);

unsigned int fm_get_len_proba(fm_t const * t);

unsigned int fm_get_len_time(fm_t const * t);

unsigned int fm_get_len_p_min(fm_t const * t);

void fm_set_method(fm_t * t, unsigned long const * value, unsigned int len);

void fm_set_proba(fm_t * t, double const * value, unsigned int len,
                  unsigned int len_p_min);

void fm_set_time(fm_t * t, double const * value, unsigned int len);

void fm_put_zero(fm_t * t);

bool fm_is_zero(fm_t const * t);

fm_t * fm_copy(fm_t const * t);

int fm_is_equal(fm_t const * c1, fm_t const * c2);

int fm_print(fm_t const *);

int fm_fprint(FILE * output_file, fm_t const * t);

#ifdef __cplusplus
}
#endif

#endif /* CADO_FM_HPP */
