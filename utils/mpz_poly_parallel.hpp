#ifndef MPZ_POLY_PARALLEL_HPP_
#define MPZ_POLY_PARALLEL_HPP_

#include <gmp.h>
#include "mpz_poly.h"
#include "mpz_polymodF.h"
#include "rootfinder.h"

/* We duplicate here the parallel interface of mpz_poly. The functions
 * are exactly the same, they're only member functions of some objects.
 *
 * The implementation in the cpp file arranges so that ::mpz_poly_mul
 * resolves to mpz_poly_notparallel_info::mpz_poly_mul
 */

/* The interface below equips both mpz_poly_notparallel_info and
 * mpz_poly_parallel_info
 */
template<typename T>
struct mpz_poly_parallel_interface {
    /* part of the mpz_poly API can use a limited amount of parallelism.
     * To use these functions, call them as member functions of the
     * parallel_info object.
     */
    void mpz_poly_sub(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
    void mpz_poly_sub_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h, mpz_srcptr m);
    void mpz_poly_mul(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
    void mpz_poly_mul_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
    void mpz_poly_divexact_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
    int mpz_poly_mod_f_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
    int mpz_poly_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m, mpz_srcptr invm);
    void mpz_poly_mul_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
    void mpz_poly_mul_mod_f (mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f);
    void mpz_poly_reduce_frac_mod_f_mod_mpz (mpz_poly_ptr num, mpz_poly_ptr denom, mpz_poly_srcptr F, mpz_srcptr m);
    void mpz_poly_div_2_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr m);
    void mpz_poly_sqr_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
    void mpz_poly_pow_mod_f_mod_ui(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, unsigned long p);
    void mpz_poly_pow_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, mpz_srcptr p);
    void mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, unsigned long a, mpz_srcptr p);
    mpz_poly* mpz_poly_base_modp_init (mpz_poly_srcptr P0, int p, unsigned long *K, int l);
    void mpz_poly_base_modp_lift(mpz_poly_ptr a, mpz_poly *P, int k, mpz_srcptr pk);

    /* polymodF functions (scaled polynomials modulo a non monic integer
     * polynomial)
     */
    void mpz_polymodF_mul (mpz_polymodF_ptr Q, mpz_polymodF_srcptr P1, mpz_polymodF_srcptr P2, mpz_poly_srcptr F);
    void mpz_poly_reducemodF(mpz_polymodF_ptr P, mpz_poly_srcptr p, mpz_poly_srcptr F);
};

struct mpz_poly_notparallel_info : public mpz_poly_parallel_interface<mpz_poly_notparallel_info> {
    /* This struct needs no fields at all. Conceivably, it could have
     * some inline decider functions that return "false", if the parallel
     * object below gains such functionality.
     */
};

struct mpz_poly_parallel_info : public mpz_poly_parallel_interface<mpz_poly_parallel_info> {
    /* fields can be added if/when the underlying implementations decide
     * to act differently depending on some runtime aspects. Decider
     * functions could be added here, possibly. */
};

extern template struct mpz_poly_parallel_interface<mpz_poly_notparallel_info>;
extern template struct mpz_poly_parallel_interface<mpz_poly_parallel_info>;


#endif	/* MPZ_POLY_PARALLEL_HPP_ */
