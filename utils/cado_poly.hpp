#ifndef CADO_POLY_H
#define CADO_POLY_H

#include <stdio.h>      // FILE
#include <gmp.h>

#ifdef __cplusplus
#include "params.h"
#endif
#include "mpz_poly.h"

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAX_DEGREE]") */
#define MAX_DEGREE 10

struct cado_poly_s {
  mpz_t n;        /* number to factor */
  double skew;    /* skewness from poly file, if given, otherwise 0. */

  int nb_polys;   /* number of polynomials used, 2 in most cases */
  mpz_poly * pols;
};
typedef struct cado_poly_s cado_poly[1];
typedef struct cado_poly_s * cado_poly_ptr;
typedef const struct cado_poly_s * cado_poly_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int cado_poly_read (cado_poly_ptr, const char *filename);
extern int cado_poly_read_stream (cado_poly_ptr, FILE *);
extern int cado_poly_read_next_poly_from_stream (cado_poly_ptr, FILE *);
extern void cado_poly_set (cado_poly_ptr p, cado_poly_srcptr q);
extern void cado_poly_swap (cado_poly_ptr p, cado_poly_ptr q);

extern void cado_poly_fprintf (FILE *, const char * prefix, cado_poly_srcptr);

extern void
cado_poly_fprintf_MurphyE (FILE *fp, const char * prefix, int side,
        double MurphyE, double bound_f, double bound_g, double area);

/* More functions for printing cado_poly are defined in polyselect/ as only
 * binaries in polyselect/ use them and some functions (like L2_skewness, ...)
 * only defined in polyselect/ are needed.
 */

extern void cado_poly_init (cado_poly_ptr);
extern void cado_poly_provision_new_poly(cado_poly_ptr);
extern void cado_poly_clear (cado_poly_ptr);
extern void cado_poly_reset (cado_poly_ptr);

extern int cado_poly_check_mapping(mpz_poly_ptr G, cado_poly_srcptr cpoly,
        mpz_srcptr N);

// Compute m as the common root of f and g mod N.
// N is taken as a third argument; it can be a strict factor of the N
// stored in the polynomial.
// If this fails, then in most of the cases we have found a factor of N
// that is given instead of m.
// The return value tells whether it worked (and then m is the common
// root) or it failed (and then m is a factor of N).
extern int cado_poly_getm(mpz_ptr, cado_poly_srcptr, mpz_srcptr);

/* Return the rational side or -1 if two algebraic side */
extern int cado_poly_get_ratside (cado_poly_srcptr);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* Same idea as for cxx_mpz and friends */
struct cxx_cado_poly {
    cado_poly x;
    cxx_cado_poly() { cado_poly_init(x); }
    cxx_cado_poly(cado_poly_srcptr f) { cado_poly_init(x); cado_poly_set(x, f); }
    static void configure_switches(cxx_param_list &) {}
    static void configure_aliases(cxx_param_list & pl) {
        param_list_configure_alias(pl, "skew", "S");
    }
    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "poly", "polynomial file");
        param_list_decl_usage(pl, "skew", "skewness");
    }

    cxx_cado_poly(cxx_param_list & pl) {
        cado_poly_init(x);
        const char *tmp;
        if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
            fprintf(stderr, "Error: -poly is missing\n");
            param_list_print_usage(pl, NULL, stderr);
            cado_poly_clear(x);
            exit(EXIT_FAILURE);
        }
        if (!cado_poly_read(x, tmp)) {
            fprintf(stderr, "Error reading polynomial file %s\n", tmp);
            cado_poly_clear(x);
            exit(EXIT_FAILURE);
        }
        /* -skew (or -S) may override (or set) the skewness given in the
         * polynomial file */
        param_list_parse_double(pl, "skew", &(x->skew));
        if (x->skew <= 0.0) {
            fprintf(stderr, "Error, please provide a positive skewness\n");
            cado_poly_clear(x);
            exit(EXIT_FAILURE);
        }
    }

    ~cxx_cado_poly() { cado_poly_clear(x); }
    cxx_cado_poly(cxx_cado_poly const & o) {
        cado_poly_init(x);
        cado_poly_set(x, o.x);
    }
    cxx_cado_poly & operator=(cxx_cado_poly const & o) {
        cado_poly_set(x, o.x);
        return *this;
    }
    cxx_cado_poly(cxx_cado_poly && o) {
        cado_poly_init(x);
        cado_poly_swap(x, o.x);
    }
    cxx_cado_poly& operator=(cxx_cado_poly && o) {
        cado_poly_swap(x, o.x);
        return *this;
    }
    operator cado_poly_ptr() { return x; }
    operator cado_poly_srcptr() const { return x; }
    cado_poly_ptr operator->() { return x; }
    cado_poly_srcptr operator->() const { return x; }

    cxx_mpz_poly const & operator[](int i) const {
        /* This is really a gross hack. Since an mpz_poly is really just
         * the same thing as an mpz_poly in memory, just return the
         * mpz_poly by pretending that it's a cxx_mpz_poly. I confess
         * it's ugly, and I'd of course like to make cxx_mpz_poly the
         * first class citizen instead, but I can't because of
         * polyselect.
         */
        ASSERT_ALWAYS(i < x->nb_polys);
        mpz_poly_srcptr p = x->pols[i];
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return *reinterpret_cast<cxx_mpz_poly const *>(p);
    }
};
#if GNUC_VERSION_ATLEAST(4,3,0)
extern void cado_poly_init(cxx_cado_poly & pl) __attribute__((error("cado_poly_init must not be called on a cado_poly reference -- it is the caller's business (via a ctor)")));
extern void cado_poly_clear(cxx_cado_poly & pl) __attribute__((error("cado_poly_clear must not be called on a cado_poly reference -- it is the caller's business (via a dtor)")));
#endif



#endif
#endif	/* CADO_POLY_H */
