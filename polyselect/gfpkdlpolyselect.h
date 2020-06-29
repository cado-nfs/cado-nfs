/**
 * \file cado-nfs/polyselect/gfpkdlpolyselect.h
 * 
 * \date 21/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief header file of gfpkdlpolyselect.c
 *        Compute two polynomials f, g suitable for discrete logarithm in 
 *        extension fields, with the conjugation method.
 *
 * \test TODO
 */

#ifndef GFPKDLPOLYSELECT_H
#define GFPKDLPOLYSELECT_H

#define ERROR_UNSPECIFIED 0
#define POLY_OK 1
#define ERROR_DEGREE 2
#define ERROR_NOT_IRREDUCIBLE 3
#define ERROR_NO_DEGREE_FACTOR 4
#define ERROR_SIGNATURE 5
#define ERROR_SUBFIELD 6
#define ERROR_SMEXP 7
#define ERROR_EASYBADIDEALS 8
#define ERROR_VERYBADIDEALS 9
#define ERROR_MAX 10

#include <stdio.h>      // FILE
#include <stdbool.h>    // for bool (in C)
#include <gmp.h>
#include "cado_poly.h"   // for MAX_DEGREE
#include "mpz_poly.h"

#define DEG_PY 2
#if DEG_PY > 2
#error "the code works only for Py of degree <= 2, sorry."
#else
//#define PHI_COEFF_INT 1 // ?

// table structure, old version.

typedef struct {
  long int t;     // parameter t
  long int PY[DEG_PY + 1]; // no --> use one of the poly stuct of cado-nfs!
  long int f[MAX_DEGREE + 1];  // polynomial f of degree at most MAX_DEGREE 
                         // set to 10 at the moment in utils/cado_poly.h
}tPyf_t;

// tables containing polynomials f
// how to encode phi form ?
typedef struct {
  unsigned int deg_f;
  unsigned int deg_Py;
  unsigned int deg_phi;
  int size; // nb of elements f in table
  const tPyf_t* tab; // la table de {f, Py, phi}
  long int phi[MAX_DEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
  //(a root of PY) modulo PY so of degree at most DEG_PY-1 --> of DEG_PY coeffs.
}tPyf_poly_t;

typedef tPyf_poly_t* tPyf_poly_ptr_t;

// table structure for CONJ with Fp2 (and JLSV1 Fp4 ?).

typedef struct {
  long int f[MAX_DEGREE + 1];  // polynomial f of degree at most MAX_DEGREE 
                              // set to 10 at the moment in utils/cado_poly.h
  long int Py[DEG_PY + 1];    // no --> use one of the poly stuct of cado-nfs!
  long int phi[MAX_DEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
}fPyphi_t;

typedef struct {
  unsigned int deg_f;
  unsigned int deg_Py;
  unsigned int deg_phi;
  int size;
  const fPyphi_t* tab;
}fPyphi_poly_t;

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
}polyselect_parameter_poly_t;

typedef polyselect_parameter_poly_t ppf_t[1];

typedef struct {
  unsigned int n; // also written k in earlier versions
  mpz_srcptr p;
  mpz_srcptr ell;
  unsigned int mnfs;
}polyselect_parameters_t;

typedef polyselect_parameters_t pp_t[1];

#ifdef __cplusplus
extern "C" {
#endif

// to convert a poly from int poly[] to mpz_poly format
void mpz_poly_set_si(mpz_poly f, const int * h, int deg_h);
void mpz_poly_setcoeff_sli(mpz_poly_ptr f, int i, long int z);
void mpz_poly_set_sli(mpz_poly f, const long int * h, int deg_h);

// works only if PY is of degree 2
void eval_mpz_phi_mpz_uv(mpz_poly g, mpz_t** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v);
void eval_si_phi_mpz_y(mpz_poly g, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], unsigned int deg_phi, mpz_t y);
void eval_si_phi_mpz_uv(mpz_poly g, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], unsigned int deg_phi, mpz_t u, mpz_t v);
// works only if PY is of degree 2

void init_eval_mpz_phi_mpz_uv(mpz_poly g, mpz_t** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v);
void init_eval_si_phi_mpz_y(mpz_poly g, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], unsigned int deg_phi, mpz_t y);
void init_eval_si_phi_mpz_uv(mpz_poly g, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], unsigned int deg_phi, mpz_t u, mpz_t v);


bool is_irreducible_ZZ(mpz_poly_srcptr phi);
bool is_irreducible_mod_p_deg2(mpz_poly_srcptr phi, mpz_srcptr p, mpz_t Discr);
bool is_irreducible_mod_p(mpz_poly_srcptr phi, mpz_srcptr p);
bool is_irreducible_mod_p_si(const long int* Py, int deg_Py, mpz_srcptr p);

// same as is_good_poly_check_all in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
int is_good_poly(pp_t pp, ppf_t ppf, mpz_poly f);
// I need a function that works also for polynomials of small coefficients (signed long int)
bool is_good_phi(mpz_poly phi, unsigned int n, mpz_t p);
bool is_good_f_PY(fPyphi_t* fPyphi, mpz_poly** phi);

// set f_id l'indice de la bonne ligne du tableau de {f, Py, phi}
bool get_f_CONJ(int* f_id, mpz_t * tab_roots_Py, int* nb_roots_Py, const fPyphi_poly_t * ff, mpz_srcptr p);
// ff->tab_f[i] is the line with a right f.

// case MNFS: tab of g_i, 1 <= i <= mnfs.
bool get_g_CONJ(mpz_poly g[], mpz_poly phi, ppf_t params_g, int f_id, mpz_t * tab_roots_Py, int nb_roots_Py, const fPyphi_poly_t * ff, pp_t params);
// ff->tab[i] is the line with a right f.
// ff->tab[i].f

/* print functions */

void mpz_poly_fprintf_cado_format_line (FILE *fp, mpz_poly f, 
					const int j, const char* label_poly);
void mpz_phi_poly_fprintf_cado_format_line (FILE *fp, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], unsigned int deg_phi, unsigned int deg_Py, int j, const char *label_poly);
void fprintf_gfpn_poly_info (FILE* fp, mpz_poly f, const char *label_poly);
void gfpk_print_params(unsigned int n, mpz_srcptr p, mpz_srcptr ell);

// the function to call for generating a .poly file.
// works only for n=2 at the moment.
// , mpz_t ell, unsigned int mnfs
int gfpkdlpolyselect( unsigned int n, mpz_srcptr p, mpz_srcptr ell, unsigned int mnfs, const char* out_filename);

#ifdef __cplusplus
}
#endif



#endif // DEG_PY > 2 is not supported
#endif // GFPKDLPOLYSELECT_H
