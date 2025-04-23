#ifndef CADO_POLYSELECT_GFPKDLPOLYSELECT_IMPL_H
#define CADO_POLYSELECT_GFPKDLPOLYSELECT_IMPL_H

#include "gfpkdlpolyselect.h"

#ifdef __cplusplus
extern "C" {
#endif

// table structure, old version.

typedef struct {
  long int t;     // parameter t
  long int PY[DEG_PY + 1]; // no --> use one of the poly stuct of cado-nfs!
  long int f[MAX_DEGREE + 1];  // polynomial f of degree at most MAX_DEGREE 
                         // set to 10 at the moment in utils/cado_poly.h
} tPyf_t;

// tables containing polynomials f
// how to encode phi form ?
typedef struct {
  int deg_f;
  int deg_Py;
  int deg_phi;
  int size; // nb of elements f in table
  const tPyf_t* tab; // la table de {f, Py, phi}
  long int phi[MAX_DEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
  //(a root of PY) modulo PY so of degree at most DEG_PY-1 --> of DEG_PY coeffs.
} tPyf_poly_t;

typedef tPyf_poly_t* tPyf_poly_ptr_t;

// table structure for CONJ with Fp2 (and JLSV1 Fp4 ?).

typedef struct {
  long int f[MAX_DEGREE + 1];  // polynomial f of degree at most MAX_DEGREE 
                              // set to 10 at the moment in utils/cado_poly.h
  long int Py[DEG_PY + 1];    // no --> use one of the poly stuct of cado-nfs!
  long int phi[MAX_DEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
} fPyphi_t;

typedef struct {
  int deg_f;
  int deg_Py;
  int deg_phi;
  int size;
  const fPyphi_t* tab;
} fPyphi_poly_t;

typedef fPyphi_poly_t* fPyphi_poly_ptr_t;

// for keeping parameters of each poly (each side) along the process of polyselect

typedef struct {
  unsigned int deg;
  unsigned int sgtr[2];
  unsigned int deg_subfield;
  unsigned int sgtr_subfield[2];
  unsigned int smexp; // from 1 to deg
  int nb_max_easybadideals;
  int nb_max_verybadideals;
} polyselect_parameter_poly_t;

typedef polyselect_parameter_poly_t ppf_t[1];

typedef struct {
  unsigned int n; // also written k in earlier versions
  mpz_srcptr p;
  /* This field seems to be mostly unused */
  mpz_srcptr ell ATTRIBUTE_DEPRECATED;
  unsigned int mnfs;
} polyselect_parameters_t;

typedef polyselect_parameters_t pp_t[1];

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_GFPKDLPOLYSELECT_IMPL_H */
