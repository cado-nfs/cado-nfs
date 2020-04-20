#include "cado.h"
#include "bblas_level3a.hpp"
#include <cstring>
#include "utils.h"

/**********************************************************************/
/* level 3a: basic operations on 64*64 matrices.
 *      mat64_eq
 *      mat64_add
 *      mat64_is_uppertriangular
 *      mat64_is_lowertriangular
 *      mat64_triangular_is_unit
 *      mat64_set_identity
 *      mat64_set_zero
 *      mat64_copy
 *      mat64_transpose     (several variants, one obvious winner)
 */

/* implemented here: see above. names are self-explanatory. */

int mat64_eq(mat64 const & a, mat64 const & b)/*{{{*/
{
    return a == b;
}
/*}}}*/
void mat64_fill_random(mat64 & w, gmp_randstate_t rstate)
{
    memfill_random(w.data(), (64) * sizeof(uint64_t), rstate);
}
void mat64_add_C(mat64 & C, mat64 const & A, mat64 const & B)/*{{{*/
{
    for(unsigned int j = 0 ; j < mat64::width ; j++) {
        C[j] = A[j] ^ B[j];
    }
}
/*}}}*/
int mat64_is_uppertriangular(mat64 const & u)/*{{{*/
{
    mat64::datatype mask = 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        if (u[k]&(mask-1)) return 0;
    }
    return 1;
}/*}}}*/
void mat64_extract_uppertriangular(mat64 & a, mat64 const & b)/*{{{*/
{
    mat64::datatype mask = 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        a[k] = b[k] & ~(mask-1);
    }
}/*}}}*/
int mat64_is_lowertriangular(mat64 const & u)/*{{{*/
{
    mat64::datatype mask = ~(mat64::datatype) 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        if (u[k] & mask) return 0;
    }
    return 1;
}/*}}}*/
void mat64_extract_lowertriangular(mat64 & a, mat64 const & b)/*{{{*/
{
    mat64::datatype mask = ~(mat64::datatype) 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        a[k] = b[k] & ~mask;
    }
}/*}}}*/
int mat64_triangular_is_unit(mat64 const & u)/*{{{*/
{
    mat64::datatype mask = 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        if (!(u[k]&mask)) return 0;
    }
    return 1;
}/*}}}*/
void mat64_triangular_make_unit(mat64 & u)/*{{{*/
{
    mat64::datatype mask = 1;
    for(unsigned int k = 0 ; k < mat64::width ; k++, mask<<=1) {
        u[k] |= mask;
    }
}
/*}}}*/

/* from hacker's delight */
void mat64_transpose_recursive_inplace(mat64 & a)/*{{{*/
{
    uint64_t m = UINT64_C(0x00000000FFFFFFFF);

    for (int j = 32 ; j ; j >>= 1, m ^= m << j) {
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            /* a big-endian-thinking version would be
             *
            uint64_t t = (a[k] ^ (a[k | j] >> j)) & m;
            a[k] ^= t;
            a[k | j] ^= t << j;
            */
            uint64_t t = (a[k] >> j ^ (a[k | j])) & m;
            a[k] ^= t << j;
            a[k | j] ^= t;
        }
    }
}/*}}}*/

/* not in place ! */
void mat64_transpose_simple_and_stupid(mat64 & dst, mat64 const & src)/*{{{*/
{
    ASSERT_ALWAYS(&dst != &src);
    for (int i = 0; i < 64; i++) {
	dst[i] = 0;
	for (int j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & UINT64_C(1)) << j;
	}
    }
}/*}}}*/

void mat64_transpose_recursive(mat64 & dst, mat64 const & src)/*{{{*/
{
    if (&dst != &src) dst = src;
    mat64_transpose_recursive_inplace(dst);
}/*}}}*/

/* {{{ final choices. These are static choices at this point, but it should
 * be the result of some tuning, ideally */
void mat64_add(mat64 & C, mat64 const & A, mat64 const & B)
{
    mat64_add_C(C,A,B);
}

void mat64_transpose(mat64 & dst, mat64 const & src)
{
    mat64_transpose_recursive(dst, src);
}
/*}}}*/

namespace bblas_bitmat_details {
    void bitmat_ops<mat64>::fill_random(mat64 & w, gmp_randstate_t rstate)
    {
        mat64_fill_random(w, rstate);
    }
    void bitmat_ops<mat64>::add(mat64 & C, mat64 const & A, mat64 const & B) {
        mat64_add(C, A, B);
    }
    void bitmat_ops<mat64>::transpose(mat64 & C, mat64 const & A) {
        mat64_transpose(C, A);
    }
    bool bitmat_ops<mat64>::is_lowertriangular(mat64 const & A) {
        return mat64_is_lowertriangular(A);
    }
    bool bitmat_ops<mat64>::is_uppertriangular(mat64 const & A) {
        return mat64_is_uppertriangular(A);
    }
    bool bitmat_ops<mat64>::triangular_is_unit(mat64 const & A) {
        return mat64_triangular_is_unit(A);
    }
    void bitmat_ops<mat64>::extract_uppertriangular(mat64 & a, mat64 const & b) {
        mat64_extract_uppertriangular(a, b);
    }
    void bitmat_ops<mat64>::extract_lowertriangular(mat64 & a, mat64 const & b) {
        mat64_extract_lowertriangular(a, b);
    }
    void bitmat_ops<mat64>::make_uppertriangular(mat64 & u) {
        extract_uppertriangular(u, u);
    }
    void bitmat_ops<mat64>::make_lowertriangular(mat64 & u) {
        extract_lowertriangular(u, u);
    }
    void bitmat_ops<mat64>::make_unit_uppertriangular(mat64 & u) {
        make_uppertriangular(u);
        triangular_make_unit(u);
    }
    void bitmat_ops<mat64>::make_unit_lowertriangular(mat64 & u) {
        make_lowertriangular(u);
        triangular_make_unit(u);
    }
    void bitmat_ops<mat64>::triangular_make_unit(mat64 & u) {
        mat64_triangular_make_unit(u);
    }
}
