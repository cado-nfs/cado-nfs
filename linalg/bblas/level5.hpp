#ifndef BBLAS_LEVEL5_HPP_
#define BBLAS_LEVEL5_HPP_

#include "bblas.hpp"

/**********************************************************************/
/* level 5: polynomials with 64*64 matrix coefficients
 *      m64pol_add
 *      m64pol_mul
 *      m64pol_addmul
 *      m64pol_addmul_kara
 *      m64pol_mul_kara
 *
 *
 */

/* lengths of a1 and a2 are n1 and n2 */
void m64pol_add(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n);

void m64pol_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);

void m64pol_addmul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);

void m64pol_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);

void m64pol_addmul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);

void m64pol_mul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2);

void m64pol_scalmul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);

void m64pol_scalmul_gf2_64_bitslice2(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);

void m64pol_mul_gf2_64_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2);

void m64pol_scalmul_gf2_64_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar);

void m64pol_mul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2);

void m64pol_scalmul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);

void m64pol_mul_gf2_128_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2);

void m64pol_scalmul_gf2_128_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar);

void m64polblock_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K);

void m64polblock_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K);

#endif	/* BBLAS_LEVEL5_HPP_ */
