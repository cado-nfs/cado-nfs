#include "cado.h"
#include "level3a.hpp"
#include <cstring>

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

int mat64_eq(mat64_srcptr a, mat64_srcptr b)/*{{{*/
{
    return memcmp(a, b, sizeof(mat64)) == 0;
}
/*}}}*/
void mat64_add_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)/*{{{*/
{
    for (int j = 0; j < 64; j++) {
        C[j] = A[j] ^ B[j];
    }
}
/*}}}*/
int mat64_is_uppertriangular(mat64_srcptr u)/*{{{*/
{
    uint64_t mask = 1;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (u[k]&(mask-1)) return 0;
    }
    return 1;
}/*}}}*/
int mat64_is_lowertriangular(mat64_srcptr u)/*{{{*/
{
    uint64_t mask = -2;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (u[k]&mask) return 0;
    }
    return 1;
}/*}}}*/
int mat64_triangular_is_unit(mat64_srcptr u)/*{{{*/
{
    uint64_t mask = 1;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (!(u[k]&mask)) return 0;
    }
    return 1;
}/*}}}*/
void mat64_set_identity(mat64_ptr m)/*{{{*/
{
    uint64_t mask = 1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) m[j]=mask;
}/*}}}*/
void mat64_copy(mat64_ptr b, mat64_srcptr a)/*{{{*/
{
    memcpy(b,a,sizeof(mat64));
}/*}}}*/
void mat64_set_zero(mat64_ptr m)/*{{{*/
{
    memset(m,0,sizeof(mat64));
}/*}}}*/


/* from hacker's delight */
void mat64_transpose_recursive_inplace(mat64_ptr a)/*{{{*/
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
void mat64_transpose_simple_and_stupid(mat64_ptr dst, mat64_srcptr src)/*{{{*/
{
    ASSERT_ALWAYS(dst != src);
    for (int i = 0; i < 64; i++) {
	dst[i] = 0;
	for (int j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & UINT64_C(1)) << j;
	}
    }
}/*}}}*/

void mat64_transpose_recursive(mat64_ptr dst, mat64_srcptr src)/*{{{*/
{
    memcpy(dst, src, sizeof(mat64));
    mat64_transpose_recursive_inplace(dst);
}/*}}}*/

/* {{{ final choices. These are static choices at this point, but it should
 * be the result of some tuning, ideally */
void mat64_add(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    mat64_add_C(C,A,B);
}

void mat64_transpose(mat64_ptr dst, mat64_srcptr src)
{
    mat64_transpose_recursive(dst, src);
}
/*}}}*/
