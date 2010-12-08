#ifndef CADO_POLY_H_
#define CADO_POLY_H_

#include <stdio.h>
#include <gmp.h>

#include "params.h"

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAXDEGREE]") */
#define MAXDEGREE 10

typedef struct
{
  mpz_t n;        /* number to factor */
  double skew;    /* skewness */
  int degree;     /* (algebraic) degree */
  mpz_t *f;       /* algebraic coefficients */
  mpz_t *g;       /* rational coefficients */
  int degreeg;    /* degree of polynomial g */
  mpz_t m;        /* common root of f and g mod n */
  char type[256]; /* type (gnfs or snfs) */
  unsigned long rlim; /* rational  factor base bound */
  unsigned long alim; /* algebraic factor base bound */
  int lpbr;           /* rational  large prime bound is 2^lpbr */
  int lpba;           /* algebraic large prime bound is 2^lpba */
  int mfbr;           /* bound for rational  residuals is 2^mfbr */
  int mfba;           /* bound for algebraic residuals is 2^mfba */
  double rlambda;     /* lambda sieve parameter on the rational  side */
  double alambda;     /* lambda sieve parameter on the algebraic side */
  int qintsize;       /* sieve block range */
} cado_poly_struct;

typedef cado_poly_struct cado_poly[1];
typedef cado_poly_struct * cado_poly_ptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void fprint_polynomial (FILE *, mpz_t *, const int);

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int cado_poly_read (cado_poly, const char *filename);
extern int cado_poly_read_stream (cado_poly, FILE *);
extern int cado_poly_set_plist(cado_poly poly, param_list pl);

extern void cado_poly_init (cado_poly);
extern void cado_poly_clear (cado_poly);

/* sanity check */
extern void check_polynomials (cado_poly);


// extern int cado_poly_read(cado_poly, char *filename) __attribute__((deprecated));
// extern int cado_poly_read_stream(cado_poly, FILE *f) __attribute__((deprecated));

/* legacy -- to be removed */
static inline int read_polynomial(cado_poly c, char *f)
    __attribute__((deprecated));
static inline int read_polynomial(cado_poly c, char *f)
{
    cado_poly_init(c);
    return cado_poly_read(c,f);
}


#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLY_H_ */
