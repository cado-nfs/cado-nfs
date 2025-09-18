/**
 * \file cado-nfs/polyselect/gfpkdlpolyselect.c
 *
 * \date 08/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief Compute two polynomials f, g suitable for discrete logarithm in
 *        extension fields, with the conjugation method.
 *
 * \test TODO
 * \version 1: dedicated for GF(p^2).
 */

#include "cado.h" // IWYU pragma: keep

#include <math.h> // log
#include <stdio.h>
#include <stdlib.h> // abort free

#include <gmp.h>

#include "auxiliary.h"
#include "cado_poly.h"
#include "gfpkdlpolyselect.h"
#include "gmp_aux.h"
#include "lll.h"    // mat_Z, LLL
#include "macros.h" // ASSERT
#include "mpz_poly.h"
#include "polyselect_alpha.h"
#include "polyselect_norms.h"
#include "rootfinder.h"

// #define PHI_COEFF_INT 1 // ?

#include "gfpkdlpolyselect_impl.h"

#include "table_f_Py_phi__f_deg4_s02_C4V4_h1_Py_s20_f01.h"
// contains: fPyphi_poly_t ff4

// works only if PY is of degree 2
// static void eval_mpz_phi_mpz_uv(mpz_poly_ptr g, mpz_t** phi_coeff, int
// deg_phi, mpz_t u, mpz_t v);
static void eval_si_phi_mpz_y(mpz_poly_ptr g,
                              long int const phi_coeff[MAX_DEGREE + 1][DEG_PY],
                              int deg_phi, mpz_t y);
static void eval_si_phi_mpz_uv(mpz_poly_ptr g,
                               long int const phi_coeff[MAX_DEGREE + 1][DEG_PY],
                               int deg_phi, mpz_t u, mpz_t v);
// works only if PY is of degree 2

// static void init_eval_mpz_phi_mpz_uv(mpz_poly_ptr g, mpz_t** phi_coeff, int
// deg_phi, mpz_t u, mpz_t v); static void init_eval_si_phi_mpz_y(mpz_poly_ptr
// g, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], int deg_phi, mpz_t y);
// static void init_eval_si_phi_mpz_uv(mpz_poly_ptr g, const long int
// phi_coeff[MAX_DEGREE + 1][DEG_PY], int deg_phi, mpz_t u, mpz_t v);

// static bool is_irreducible_ZZ(mpz_poly_srcptr phi);
static bool is_irreducible_mod_p_deg2(mpz_poly_srcptr phi, mpz_srcptr p,
                                      mpz_t Discr);
static bool is_irreducible_mod_p(mpz_poly_srcptr phi, mpz_srcptr p);
static bool is_irreducible_mod_p_si(long int const * Py, int deg_Py,
                                    mpz_srcptr p);

// same as is_good_poly_check_all in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
// static int is_good_poly(pp_t pp, ppf_t ppf, mpz_poly f);
// I need a function that works also for polynomials of small coefficients
// (signed long int) static bool is_good_phi(mpz_poly phi, unsigned int n, mpz_t
// p); static bool is_good_f_PY(fPyphi_t* fPyphi, mpz_poly** phi);

// set f_id l'indice de la bonne ligne du tableau de {f, Py, phi}
static bool get_f_CONJ(int * f_id, mpz_t * tab_roots_Py, int * nb_roots_Py,
                       fPyphi_poly_t const * ff, mpz_srcptr p);
// ff->tab_f[i] is the line with a right f.

// case MNFS: tab of g_i, 1 <= i <= mnfs.
static bool get_g_CONJ(mpz_poly g[], mpz_poly phi, ppf_t params_g, int f_id,
                       mpz_t * tab_roots_Py, int nb_roots_Py,
                       fPyphi_poly_t const * ff, pp_t params);
// ff->tab[i] is the line with a right f.
// ff->tab[i].f

/* print functions */

static void mpz_poly_fprintf_cado_format_line(FILE * fp, mpz_poly f, int j,
                                              char const * label_poly);
// static void mpz_phi_poly_fprintf_cado_format_line (FILE *fp, const long int
// phi_coeff[MAX_DEGREE + 1][DEG_PY], int deg_phi, unsigned int deg_Py, int j,
// const char *label_poly);
static void fprintf_gfpn_poly_info(FILE * fp, mpz_poly_srcptr f,
                                   char const * label_poly);

/**
 * \brief Return appropriate degrees of polynomials f and g for CONJUGATION
 * method according to extension degree field n > 1.
 *
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[out] deg_f degree of polynomial f, this is 2*n
 * \param[out] deg_g degree of polynomial g, this is n
 *
 * note that deg_f >= deg_g by general convention in all this cado-nfs library.
 */
static void get_degree_CONJ_f_g(int n, int * deg_f, int * deg_g)
{
    *deg_f = 2 * n;
    *deg_g = n;
}

/**
 * \brief evaluate phi at (u,v) and outputs a polynomial g,
 *        assuming that Py is of degree 2.
 *
 * \param[in] phi_coeff double entry array. The poly phi is expressed as
 * (ZZ[y])[X] i.e. sum_{i=0}^{deg phi} ((a_{i0} + a_{i1}*y)*X^i)
 * \param[in] deg_phi degree of polynomial phi, for the moment this must be less
 * than 2.
 * \param[in] u  mpz_t integer, of size half the size of p
 * \param[in] v  mpz_t integer, of size half the size of p
 *               u/v = y mod p with y a root of PY mod p.
 * \param[out] g polynomial g of same degree as phi (already initialized)
 *
 * NOTE: maybe input phi and deg_phi directly ? Or a poly structure ?
 */
// works only if PY is of degree 2
/* 20240905 -- unused, removing
void eval_mpz_phi_mpz_uv(mpz_poly_ptr g, mpz_t** phi_coeff, int deg_phi, mpz_t
u, mpz_t v)
{
  unsigned int i;
  for (i=0; i <= deg_phi; i++){
    mpz_mul(mpz_poly_coeff(g, i), v, phi_coeff[i][0]);    // gi <- phi_i0 * v
    mpz_addmul(mpz_poly_coeff(g, i), u, phi_coeff[i][1]); // gi <- gi + phi_i1 *
u
  }
  mpz_poly_cleandeg(g, deg_phi); // set deg to deg_phi at most (check that coeff
is not 0)
}
*/

// t[0] + t[1]*X + ... + t[deg_phi]*X^deg_phi
// with t[i] = t[i][0] + t[i][1]*Y + ... t[i][deg_Py - 1]*Y^(deg_Py-1)
// t[i][j] are int
// u, v are mpz_t (multi precision integers)
// ff.phi[i][0] * v +
// ff.phi[i][1] * u
// phi_i = phi_i0 + phi_i1 * Y
// the function does not use modular arithmetic but exact integer arithmetic
static void eval_si_phi_mpz_uv(mpz_poly_ptr g,
                               long int const phi_coeff[MAX_DEGREE + 1][DEG_PY],
                               int deg_phi, mpz_t u, mpz_t v)
{
    for (int i = 0; i <= deg_phi; i++) {
        // phi has coefficients of type (signed) long int
        mpz_mul_si(mpz_poly_coeff(g, i), v,
                   phi_coeff[i][0]); // gi <- phi_i0 * v
        mpz_addmul_si(mpz_poly_coeff(g, i), u,
                      phi_coeff[i][1]); // gi <- gi + phi_i1 * u
    }
    mpz_poly_cleandeg(g, deg_phi);
}

void eval_si_phi_mpz_y(mpz_poly_ptr g,
                       long int const phi_coeff[MAX_DEGREE + 1][DEG_PY],
                       int deg_phi, mpz_t y)
{
    for (int i = 0; i <= deg_phi; i++) {
        // phi has coefficients of type (signed) long int
        // mpz_add_si does not exist.
        mpz_set_si(mpz_poly_coeff(g, i), phi_coeff[i][0]); // gi <- phi_i0
        mpz_addmul_si(mpz_poly_coeff(g, i), y,
                      phi_coeff[i][1]); // gi <- phi_i0 + phi_i1 * y
    }
    mpz_poly_cleandeg(g, deg_phi);
}

/**
 * \brief return a pointer to a table [{t, PY, f}] with f of degree n.
 *
 * \param[in] deg_f integer greater than 1 ( deg_f = 4 or 6 implemented for now)
 * \param[out] ff (pointer to a) table of polynomials f (the one of higher
 * degree) of degree deg_f, NULL if there is no such table for deg_f.
 *
 * the same table is used for GF(p^2) with CONJ and GF(p^4) with JLSV1.
 */
static bool polygen_CONJ_get_tab_f(unsigned int const deg_f,
                                   fPyphi_poly_t const ** ff)
{
    if (deg_f == 4) {
        *ff = &(fPyphi_4);
        /*
      } else if (deg_f == 6) {
        ff = TAB_f6_CYCLIC;
        */
    } else {
        ff = NULL;
    }
    return ff != NULL;
}

#if 0 /* unused code */
/**
 * \brief test irreducibility over integers of polynomial phi
 * 
 * \param[in] phi: polynomial
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 * 
 * 1st version : for Degree(phi) <= 2.
 */
static bool is_irreducible_ZZ(mpz_poly_srcptr phi)
{
  /* if deg <= 1, return true
     if deg == 2, test discriminant 
        --> if D < 0, phi is irreducible
        --> if D = 0, phi is a square and is not irreducible
        --> if D > 0 is a square, then phi has two degree 2 factors
	             otherwise phi is irreducible
     if deg >= 3: maybe use the method of finding a factorization modulo 
     small primes, then conclude according to that
  */

  int degree, sign_D, exact_sqrt;
  bool is_irreducible = false;
  mpz_t D, tmp;
  degree = phi->deg;
  if (degree <= 1){
    is_irreducible = true;
  }
  else{
    if (degree == 2){
      //compute discriminant, phi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_mul(D, mpz_poly_coeff_const(phi, 1), mpz_poly_coeff_const(phi, 1)); //   D <- b^2
      mpz_mul(tmp, mpz_poly_coeff_const(phi, 0), mpz_poly_coeff_const(phi, 2)); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and phi is reducible
	is_irreducible = true;
      }else{
	if (sign_D < 0){
	  is_irreducible = false;
	}else{
	  mpz_sub(D, D, tmp); //   D <- b^2 - 4*a*c and is > 0
	  // is_square ?
	  exact_sqrt = mpz_perfect_square_p(D);
	  if (exact_sqrt == 0){
	    is_irreducible = false;
	  }else{
	    is_irreducible = true;
	  }
	}// sign_D < 0
      }// sign_D == 0
      mpz_clear(D);
      mpz_clear(tmp);
    }// degree == 2
  }// degree <= 1
  return is_irreducible;
}
#endif

/**
 * \brief test irreducibility modulo a prime p of polynomial phi
 *
 * \param[in] phi: polynomial with large coefficients
 * \param[in] p: (large) prime number
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 *
 * 1st version : for Degree(phi) <= 2.
 */
bool is_irreducible_mod_p_deg2(mpz_poly_srcptr phi, mpz_srcptr p, mpz_t Discr)
{
    /* if deg <= 1, return true
       if deg == 2, test discriminant
          --> if D is a square mod p, (compute Legendre symbol)
              then phi has two degree 2 factors
              otherwise phi is irreducible
       if deg >= 3: use MPFQ library ?
    */
    int degree, Legendre_D_p;
    bool is_irreducible = false;
    mpz_t tmp;
    degree = phi->deg;
    if (degree <= 1) {
        is_irreducible = true;
    } else {
        if (degree == 2) {
            // compute discriminant, phi = a*x^2 + b*x + c
            // mpz_init(Discr); assume already done in calling function
            mpz_init(tmp);
            mpz_mul(Discr, mpz_poly_coeff_const(phi, 1),
                    mpz_poly_coeff_const(phi, 1)); //   D <- b^2
            mpz_mul(tmp, mpz_poly_coeff_const(phi, 0),
                    mpz_poly_coeff_const(phi, 2)); // tmp <- a*c
            mpz_mul_ui(tmp, tmp, 4);               // tmp <- 4*a*c
            mpz_sub(Discr, Discr, tmp); //   D <- b^2 - 4*a*c and is != 0
            if (mpz_sgn(Discr) == 0) {
                // very easy, D = 0 and phi is reducible
                is_irreducible = false;
            } else {
                // is_square mod p ?
                // compute Legendre symbol
                Legendre_D_p = mpz_legendre(Discr, p);
                if (Legendre_D_p >=
                    0) { // indeed, (D/p) can be 0 if D is a multiple of p.
                    is_irreducible = false;
                } else { // (D/p) = -1
                    is_irreducible = true;
                }
            }
            mpz_clear(tmp);
        } else {
            // degree > 2, use a function that tests irreducibility. But the
            // function does not exist yet.
            fprintf(stderr, "Implement me !\n");
            abort();
        }
    }
    return is_irreducible;
}

/**
 * \brief test irreducibility modulo a prime p of polynomial phi
 *
 * \param[in] phi: polynomial with large coefficients
 * \param[in] p: (large) prime number
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 *
 * 1st version : for Degree(phi) <= 2.
 */
bool is_irreducible_mod_p(mpz_poly_srcptr phi, mpz_srcptr p)
{
    /* if deg <= 1, return true
       if deg == 2, test discriminant
          --> if D is a square mod p, (compute Legendre symbol)
              then phi has two degree 2 factors
              otherwise phi is irreducible
       if deg >= 3: use MPFQ library ?
    */
    int degree;
    bool is_irreducible = false;
    mpz_t D;

    degree = phi->deg;
    if (degree <= 1) {
        is_irreducible = true;
    } else if (degree == 2) {
        mpz_init(D);
        is_irreducible = is_irreducible_mod_p_deg2(phi, p, D);
        mpz_clear(D);
    } else {
        // degree > 2, use a function that tests irreducibility. But the
        // function does not exist yet.
        fprintf(stderr, "Implement me !\n");
        abort();
    } // degree == 2
    return is_irreducible;
}

// same as above but with a poly with "small" coeffs given in a tab, with size
// <-> deg.
bool is_irreducible_mod_p_si(long int const * f, int deg_f, mpz_srcptr p)
{
    int Legendre_D_p;
    bool is_irreducible = false;
    mpz_t D, tmp;
    mpz_t a, b, c;
    if (deg_f <= 1) {
        is_irreducible = true;
    } else {
        if (deg_f == 2) {
            // compute discriminant, phi = a*x^2 + b*x + c
            mpz_init(D);
            mpz_init(tmp);
            mpz_init_set_si(a, f[0]);
            mpz_init_set_si(b, f[1]);
            mpz_init_set_si(c, f[2]);
            mpz_mul(D, b, b);        //   D <- b^2
            mpz_mul(tmp, a, c);      // tmp <- a*c
            mpz_mul_ui(tmp, tmp, 4); // tmp <- 4*a*c
            mpz_sub(D, D, tmp);      //   D <- b^2 - 4*a*c and is != 0
            if (mpz_sgn(D) == 0) {
                // very easy, D = 0 and phi is reducible
                is_irreducible = false;
            } else {
                // is_square mod p ?
                // compute Legendre symbol
                Legendre_D_p = mpz_legendre(D, p);
                if (Legendre_D_p >=
                    0) { // indeed, (D/p) can be 0 if D is a multiple of p.
                    is_irreducible = false;
                } else { // (D/p) = -1
                    is_irreducible = true;
                }
            }
            mpz_clear(D);
            mpz_clear(tmp);
            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(c);
        } else {
            // degree > 2, use a function that tests irreducibility. But the
            // function does not exist yet.
            fprintf(stderr, "Implement me !\n");
            abort();
        } // degree == 2
    } // degree <= 1
    return is_irreducible;
}

// [almost] same as is_good_poly in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
// int is_good_poly_light(pp_t params, ppf_t poly_params, mpz_poly f){}

/* check:
   + deg f
   + f irreducible over Q
   + f has a degree params->n factor mod params->p
   o [ Galois group of f: assume that the table is already checked.]
   o [ h=1 ]
   + signature(f) : how to do that in C ??? [look at the Cohen, then ask
   Pierrick and Paul.] o signature subfield: how to find the subfield ? [id] o
   badideals: call the magma function... o smexp: does this already exist in C ?
   call magma maybe. [ yes, look up sm_side_info_init in sm_utils ]
*/

// same as is_good_poly_check_all in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
int is_good_poly_deg2(pp_t params, ppf_t poly_params, mpz_poly_srcptr g)
{
    /* check:
       - deg f
       - f irreducible over Q
       - f has a degree params->n factor mod params->p
         for CONJ: check
       - [ Galois group of f: assume that the table is already checked.]
       - [ h=1 ]
       - signature(f) : how to do that in C ??? [look at the Cohen, then ask
       Pierrick and Paul.]
       - signature subfield: how to find the subfield ? [id]
       - badideals: call the magma function...
       - smexp: does this already exists in C ? call magma maybe.
          [ yes, look up sm_side_info_init in sm_utils ]
    */

    bool sgtr_ok = false;
    mpz_t D;
    mpz_init(D);
    bool is_irreducible = is_irreducible_mod_p_deg2(g, params->p, D);
    if (!is_irreducible) {
        mpz_clear(D);
        return false;
    }

    // check signature of g --> D < 0 or D > 0 ?
    // get signature:
    unsigned int * sgtr_g = poly_params->sgtr;
    if ((sgtr_g[0] == 2) && (sgtr_g[1] == 0)) {
        sgtr_ok = mpz_sgn(D) >= 0;
    } else {
        if ((sgtr_g[0] == 0) && (sgtr_g[1] == 1)) {
            sgtr_ok = mpz_sgn(D) < 0;
        }
    }
    mpz_clear(D);
    return (is_irreducible && sgtr_ok);
}

/**
 * \brief test wether phi is a suitable candidate for computing g from phi
 *
 * \param[in] p multiprecision integer, prime
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[in] phi polynomial with coefficients mod p (so very large)
 *
 * \return bool true if phi is of degree n and irreducible mod p
 * \return false otherwise
 */
// bool is_good_phi(mpz_poly phi, unsigned int n, mpz_t p){}
//     return (Degree(phi_p) eq n) and IsIrreducible(phi_p);

/**
 * \brief select suitable polynomial f
 *
 * \param[in]  p            multiprecision integer, prime
 * \param[in]  ff           ptr to struct fPyphi_poly_t (ff->tab[i].f for i-th
 * f)
 * \param[out] f_id         index of first good f in ff->tab. If end of tab is
 * reached, then f_id is out of tab size. (e.g. ff->size)
 * \param[out] tab_phi      table of possible phi (at most ff->deg_Py) for the
 * f_id-th f in ff->tab
 * \param[out] tab_roots_Py tab of roots y of Py corresponding to f
 *
 * phi is a degree n factor of f, with coefficients in y with y a root of Py.
 * \return bool
 */

bool get_f_CONJ(int * f_id, mpz_t * tab_roots_Py, int * nb_roots_Py,
                fPyphi_poly_t const * ff, mpz_srcptr p)
{
    int i = 0;
    int k = 0;
    bool found_good_f = false;
    mpz_t * y;
    mpz_poly phiy;         // phi evaluated at a given y
    mpz_poly Pyi_mpz_poly; // Py with mpz_t coeffs

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (ff != NULL) {
        // Now, find an appropriate poly f in that table.
        // start with the first one, etc because the polynomials are sorted in
        // decreasing order of interest (decreasing order of Murphy E value)
        mpz_poly_init(Pyi_mpz_poly, ff->deg_Py);
        mpz_poly_init(phiy, ff->deg_Py);
        y = (mpz_t *)malloc(ff->deg_Py * sizeof(mpz_t));
        for (int j = 0; j < ff->deg_Py; j++) {
            mpz_init(y[j]);
        }
        i = *f_id;
        while ((found_good_f != true) && (i < ff->size)) { // nothing found

            if (!is_irreducible_mod_p_si((ff->tab[i]).Py, ff->deg_Py, p)) { //
                // Py is not irreducible i.e. if deg=2, means has roots mod p.
                // careful:: if one day, Py will be of degree > 2, then will
                // need here find_root(Py) instead.
                mpz_poly_setcoeffs_si(Pyi_mpz_poly, (ff->tab[i]).Py,
                                      ff->deg_Py); // long int -> int bof
                // compute the two roots of Pyi,
                // --> compute a square root mod p
                // HOW TO DO THAT ?
                // we obtain y a root of Py mod p.
                // 1st arg: tab of roots (ptr)
                // 2nd: the poly to find the roots
                // 3rd: p (prime)
                unsigned int nb_roots_y =
                    mpz_poly_roots(y, Pyi_mpz_poly, p, rstate);
                // this function calls either a ulong version or an mpz_t
                // version. nb_roots_y = mpz_poly_roots_ulong (y, Pyi_mpz_poly,
                // p); nb_roots_y = mpz_poly_roots_mpz   (y, Pyi_mpz_poly, p);

                // compute the possible phi(s) foreach y and store them (good
                // phi, good corresponding y)
                k = 0;
                for (unsigned int l = 0; l < nb_roots_y; l++) {
                    eval_si_phi_mpz_y(phiy, (ff->tab[i]).phi, ff->deg_phi,
                                      y[l]);
                    // phiy is now a poly in one variable with large coeffs.
                    // test if at least one phiy is irreducible
                    if (is_irreducible_mod_p(phiy, p)) {
                        found_good_f = true;
                        mpz_init_set(tab_roots_Py[k], y[l]);
                        k++;
                    } // end if
                } // end for loop over the roots y of Py mod p
            } // Py had roots mod p
            i++; // loop over tab ff->tab of polys {f, Py, phi}
        } // either poly f found or end of ff->tab reached.
        mpz_poly_clear(Pyi_mpz_poly);
        mpz_poly_clear(phiy);
        for (int j = 0; j < ff->deg_Py; j++) {
            mpz_clear(y[j]);
        }
        free(y);
    }
    gmp_randclear(rstate);
    if (found_good_f) {
        *f_id = i - 1;
    } else {
        *f_id = i;
    }
    *nb_roots_Py = k;
    return found_good_f;
}

// add mpz_poly phi maybe in other versions, for n>2.
// no optimization on g0, g1 here.
// g0, g1 are output of LLL.
static bool get_g_phi_y(mpz_poly g[], int * nb_found, ppf_t params_g, int f_id,
                        mpz_t y, fPyphi_poly_t const * ff, pp_t params)
{

    mpz_srcptr p = params->p;
    unsigned int n = params->n;
    unsigned long int skewness = 12;
    unsigned long int gcd_uv = 1;
    mpz_t u, v, tmp;
    mat_Z M;
    bool gi_ok = false, found = false;
    int i = 0, nb_g_found = 0;

    if (n == 2) {
        // f = x^4 + 1, phi = x^2 + y*x + 1, g= v*x^2 + u*x + v, we want u^2 -
        // 4v^2 < 0
        /* case 1: phi = x^2 + y*x + 1,
           we need u^2 - 4v^2 < 0 i.e. |u| < 2 |v|
           case 2: phi = x^2 + x + y, g = v*x^2 + v*x + u
           we need v^2 - 4uv < 0 <=> v(v-4u) < 0
           so v>0, u>0, v-4u < 0 -> v < 4u
           etc.
         */
        // TODO: add skewness until g0 has right signature params_g->sgtr
        // will skewness change anything if g is not irreducible ?

        LLL_init(&M, 2, 2);
        mpz_init(tmp);
        /// INDICES START AT 1, NOT AT 0.
        mpz_set(M.coeff[1][1], p);
        mpz_set_si(M.coeff[1][2], 0);

        mpz_mul_ui(M.coeff[2][1], y, skewness);
        mpz_set_si(M.coeff[2][2], 1);
        /* [ p   0
             s*y 1 ] with s*y = skewness*y  */
        // a/b = 3/4 for LLL, use u/v to store them.
        mpz_init_set_si(u, 3);
        mpz_init_set_si(v, 4);
        /* reduce with LL the number skewness*y -> u/v mod p */
        LLL(tmp, M, NULL, u, v);
        /* [ u0 v0
             u1 v1 ]  with gcd(ui, vi) = 1. */
        // g0 <- v0 * phi(y <- u0/v0)
        // g1 <- v1 * phi(y <- u1/v1)

        for (i = 1; i <= 2; i++) { // careful, M coeffs start at i=1, not i=0.
            // u <-> M.coeff[i][1], v <-> M.coeff[i][2]*sy
            gcd_uv = mpz_gcd_ui(tmp, M.coeff[i][1], skewness);
            // (unsigned long int) gcd_uv <- gcd(v*s, u) = gcd(s, u) since
            // gcd(v,u) = 1
            ASSERT((skewness % gcd_uv) == 0);
            if (gcd_uv > 1) {
                skewness = skewness / gcd_uv;
                mpz_tdiv_q_ui(u, M.coeff[i][1], gcd_uv);
            } else {
                mpz_set(u, M.coeff[i][1]);
            }
            if (skewness > 1) {
                mpz_mul_ui(v, M.coeff[i][2], skewness);
            } else {
                mpz_set(v, M.coeff[i][2]);
            }
            // we need v > 0
            if (mpz_sgn(v) < 0) {
                mpz_neg(v, v);
                mpz_neg(u, u);
            }
            eval_si_phi_mpz_uv(g[i - 1], ((ff->tab)[f_id]).phi, ff->deg_phi, u,
                               v);

            gi_ok = is_good_poly_deg2(params, params_g, g[i - 1]);
            if (gi_ok) {
                nb_g_found++;
            }
        }
        *nb_found = nb_g_found;
        found = nb_g_found > 0;
        LLL_clear(&M);
        mpz_clear(tmp);
        mpz_clear(u);
        mpz_clear(v);

    } // n == 2
    else {
        *nb_found = 0;
        found = false; // n> 2 not suported yet
    }
    return found;
}
bool get_g_CONJ(mpz_poly g[], mpz_poly phi, ppf_t params_g, int f_id,
                mpz_t * tab_roots_Py, int nb_roots_Py, fPyphi_poly_t const * ff,
                pp_t params)
{
    bool found_g = false;
    int i, nb_found;
    mpz_poly * gi;
    int nb_LLL_poly = 2;
    /* for each phiy:
       reduce y with lll (add skewness if needed to get the right signature)
       -> (u,v)
       evaluate phi (ff->tab_f[f_id].phi ) at (u,v)

     */
    if (ff->deg_Py == 2) {
        gi = (mpz_poly *)malloc(nb_LLL_poly * sizeof(mpz_poly));
        for (i = 0; i < nb_LLL_poly; i++) {
            mpz_poly_init(gi[i], ff->deg_phi);
        }
        // case where we can simply run LLL(y) -> (u,v) then evaluate phi at
        // (u,v). one pair of possible g per phi. (u0,v0), (u1,v1). for each
        // phi_y:
        //     compute the pair (g0,g1) --> can add skewness in LLL to force a
        //     good signature.
        // v1: keep the 1st couple (g0,g1) found.
        i = 0;
        while ((i < nb_roots_Py) && (!found_g)) {
            found_g = get_g_phi_y(gi, &nb_found, params_g, f_id,
                                  tab_roots_Py[i], ff, params);
            i++;
        }
        // now, we have a couple (g0, g1) for a given root y of Py, of index i--
        i--;
        // search for good gi, 1 <= i <= mnfs. --> NEED A SEPARATE FUNCTION.
        if (found_g) {
            eval_si_phi_mpz_y(phi, ((ff->tab)[f_id]).phi, ff->deg_phi,
                              tab_roots_Py[i]);
            for (i = 0; (i < nb_found) && (i < (int)params->mnfs); i++) {
                mpz_poly_set(g[i], gi[i]);
            }
        }
        for (i = 0; i < nb_LLL_poly; i++) {
            mpz_poly_clear(gi[i]);
        }
        free(gi);
    } else {
        // use LLL over phi directly to get g.
        printf("# ff->deg_Py = %d but only ff->deg_Py = 2 is implemented for "
               "the moment.\n",
               ff->deg_Py);
        found_g = false; // not yet implemented
    }
    return found_g;
}

/**
 * \brief The main function that computes f and g and print them in the .poly
 * file.
 *
 * \param[in] p multiprecision integer, prime
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[in] label tring to label the output .poly file
 *
 */

// , mpz_t ell, unsigned int mnfs
int gfpkdlpolyselect(int n, mpz_srcptr p, mpz_srcptr ell MAYBE_UNUSED,
                     int const mnfs, char const * out_filename)
{
    int done = 0;
    int deg_f = 2 * n, deg_g = 2 * n;
    // take the largest possibility as default init for deg_f and deg_g.
    mpz_poly f;

    fPyphi_poly_t const * ff;
    mpz_poly phi;
    mpz_t tab_roots_Py[DEG_PY]; // a tab for the roots of Py
    int nb_roots_Py;
    pp_t params;
    ppf_t params_g = {{2, {0, 1}, 0, {0, 0}, 2, -1, 0}};
    params->n = n;
    params->p = p;
    // params->ell = ell;
    params->mnfs = mnfs;

    if (n != 2) {
        /* not supported */
        return 0;
    }

    get_degree_CONJ_f_g(n, &deg_f, &deg_g);
    polygen_CONJ_get_tab_f(deg_f, &ff);

    if (ff == NULL) {
        printf("\n# pb ff == NULL.\n");
        // NO ff available for this n.
        return 0;
    }

    // ff pointe sur une structure qui contient tout ce qu'il faut : une
    // table de {f, Py, phi} plus quelques autres parametres
    mpz_poly_init(f, deg_f);
    mpz_poly * g = (mpz_poly *)malloc(mnfs * sizeof(mpz_poly));
    FATAL_ERROR_CHECK(g == NULL,
            "not enough memory to allocate for table of g.");
    for (int i = 0; i < mnfs; i++) {
        mpz_poly_init(g[i], ff->deg_phi);
    }
    mpz_poly_init(phi, n);

    for(int f_id = 0 ; f_id < ff->size && !done ; f_id++) {
        // 1. find a good f in table. If no f is found, return failed and
        // exit 0. tab_roots_Py will be mpz_init() inside the get_f_CONJ
        // function.
        if (!get_f_CONJ(&f_id, tab_roots_Py, &nb_roots_Py, ff, p))
            continue;

        // 2. find a good g with LLL.
        if (get_g_CONJ(g, phi, params_g, f_id, tab_roots_Py,
                    nb_roots_Py, ff, params))
        {

            mpz_poly_setcoeffs_si(f, ff->tab[f_id].f, ff->deg_f);
            FILE * outputpoly;
            if (out_filename != NULL) {
                outputpoly = fopen(out_filename, "w");
                ASSERT_ALWAYS(outputpoly != NULL);
            } else {
                outputpoly = stdout;
            }
            gmp_fprintf(outputpoly, "n: %Zd\n", p);
            fprintf(outputpoly, "skew: 1.0\n");
            mpz_poly_fprintf_cado_format_line(outputpoly, f, 0,
                    "f");
            fprintf_gfpn_poly_info(outputpoly, f, "f");

            for (int i = 0; i < mnfs; i++) {
                mpz_poly_fprintf_cado_format_line(outputpoly, g[i],
                        i + 1, "g");
                fprintf_gfpn_poly_info(outputpoly, g[i], "g");
            }
            fprintf(outputpoly, "# gcd(f, g) = phi = ");
            for (int i = 0; i <= phi->deg; i++) {
                gmp_fprintf(outputpoly, "%Zd",
                        mpz_poly_coeff_const(phi, i));
                if (i != phi->deg) {
                    fprintf(outputpoly, ",");
                } else {
                    fprintf(outputpoly, "\n");
                }
            }

            fclose(outputpoly);
            done = 1;
        }
        for (int k = 0; k < nb_roots_Py; k++)
            mpz_clear(tab_roots_Py[k]);
    }

    mpz_poly_clear(f);
    mpz_poly_clear(phi);
    for (int i = 0; i < mnfs; i++) {
        mpz_poly_clear(g[i]);
    }
    free(g);

    return done;
}

// print functions

/* Print f of degree d with the following format
    poly<j>: f0,f1,f2,...,fd\n
    (new format decided in 2015 for DL in GF(p^n))
*/
void mpz_poly_fprintf_cado_format_line(FILE * fp, mpz_poly f, int const j,
                                       char const * label_poly)
{
    if (label_poly != NULL) {
        fprintf(fp, "# ");
        fputs(label_poly, fp);
        fprintf(fp, "\n");
    }
    fprintf(fp, "poly%d: ", j);
    for (int i = 0; i < f->deg; i++) {
        gmp_fprintf(fp, "%Zd,", mpz_poly_coeff_const(f, i));
    }
    // i = f->deg;
    if (f->deg >= 0) {
        gmp_fprintf(fp, "%Zd\n", mpz_poly_lc(f));
    }
}

#if 0

/* This code is unused, and deserves a close look. What is the semantics
 * of deg_Py exactly ? Is it 2 ? Is it a **degree** ? Is it a length
 * (degree-1) ? The loop that goes up to deg_Py-1 below is suspicious,
 * unless deg_Py is actually only a length. And the simultaneous info
 * DEG_PY // deg_Py is a source for confusion, too. Do we have the
 * implicit assert that deg_Py <= DEG_PY ? (then it should probably be an
 * ASSERT_ALWAYS */

static void print_coeff_Y_phi(FILE *fp, const long int phi_coeff_y[DEG_PY], unsigned int deg_Py){
  fprintf (fp, "(");
  fprintf (fp, "%ld", phi_coeff_y[0]);
  if (deg_Py >= 1){
    fprintf (fp, " + %ld*Y", phi_coeff_y[1]);
  }
  for (unsigned int k=2; k < deg_Py; k++){
    fprintf (fp, " + %ld*Y^%d", phi_coeff_y[k], k);
  }
  fprintf (fp, ")");
}


void mpz_phi_poly_fprintf_cado_format_line (FILE *fp, const long int phi_coeff[MAX_DEGREE + 1][DEG_PY], int deg_phi, unsigned int deg_Py, int j, const char *label_poly)
{
  if (label_poly != NULL){
    fprintf (fp, "# ");
    fputs (label_poly, fp);
    fprintf (fp, "\n");
  }
  fprintf (fp, "poly-phi-%d:", j);

  print_coeff_Y_phi(fp, phi_coeff[0], deg_Py);
  if (deg_phi >= 1){
    fprintf (fp, " + ");    
    print_coeff_Y_phi(fp, phi_coeff[1], deg_Py);
    fprintf (fp, "*X");
  }
  for (unsigned int i = 2; i <= deg_phi; i++)
  {
    fprintf (fp, " + ");    
    print_coeff_Y_phi(fp, phi_coeff[i], deg_Py);
    fprintf (fp, "*X^%d", i);
  }
  fprintf(fp, "\n");
}
#endif

void fprintf_gfpn_poly_info(FILE * fp, mpz_poly_srcptr f,
                            const char * label_poly)
{
    // unsigned int i;
    double skew, logmu, alpha;
    double const exp_rot[] = {
        0,   0.5, 1.0, 2.0,
        3.0, 4.0, 5.0, 0}; // ??? copy paste from dlpolyselect.c:26
    skew = L2_skewness(f);
    logmu = L2_lognorm(f, skew);
    alpha = get_alpha(f, get_alpha_bound());
    fprintf(fp, "# ");
    if (label_poly != NULL) {
        fputs(label_poly, fp);
    }
    fprintf(fp,
            " lognorm %1.2f, skew %1.2f, alpha %1.2f, E %1.2f, "
            "exp_E %1.2f\n",
            logmu, skew, alpha, logmu + alpha,
            logmu + expected_alpha(exp_rot[f->deg] * log(skew)));
}

void gfpk_print_params(int n, mpz_srcptr p, mpz_srcptr ell)
{
    printf("p = ");
    gmp_printf("%Zd\n", p);
    printf("ell = ");
    gmp_printf("%Zd\n", ell);
    printf("n = %d\n", n);
}
