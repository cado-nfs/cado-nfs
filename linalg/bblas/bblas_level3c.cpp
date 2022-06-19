#include "cado.h" // IWYU pragma: keep
#include <cstring>
#include "bblas_level2a.hpp"  // for addmul_To64_o64
#include "bblas_level2b.hpp"  // for mul_o64_6464, mul_o64_T6464
#include "bblas_level3a.hpp"  // for mat64_transpose
#include "bblas_level3c.hpp"
#include "bblas_mat64.hpp"    // for mat64
#include "bblas_simd.hpp"
#include "macros.h"                        // for MAYBE_UNUSED
// the whole point of bblas_simd is to avoid including these files...
// IWYU pragma: no_include <emmintrin.h>
// IWYU pragma: no_include <immintrin.h>

/**********************************************************************/
/* level 3c: operations on matrices with one arbitrary length (sometimes 2)
 *      copy_N64
 *      cmp_N64
 *      mul_N64_6464   (several variants)
 *      mul_N64_T6464   (several variants)
 *      addmul_N64_6464   (several variants)
 *      mul_TN64_N64
 *      addmul_TN64_N64
 * with varying left length:
 *      mul_TN32_N64
 *      mul_TN64K_N64
 *      addmul_TN64K_N64
 */

/* implemented here:
 *    - mul_N64_6464: Given an N*64 matrix A (N uint64_t's) and a 64*64
 *      matrix B, compute the product (in place).  With respect to
 *      endianness, we match the column of (A[i]&1)'s with B[0].
 *    - mul_N64_T6464: same, but multiply by the transpose of the second
 *      operand.
 *    - mul_TN64_N64: rank-N update. More wordy: this takes, in row major
 *      order, an Nx64 matrix A (transpose of a 64xN matrix), together
 *      with another Nx64 matrix B, and xors the output matrix with the
 *      product transpose(A)*B -- this may as well be seen as the block
 *      dot product of A and B.
 *    - mul_TN32_N64: rank-N update, but creates a matrix of size 32*N
 *    - mul_TN64K_N64: rank-N update, but creates a matrix of size 64K*N
 */

void copy_N64(uint64_t * dst, uint64_t const * src, size_t m)/*{{{*/
{
    memcpy(dst, src, m * sizeof(uint64_t));
}/*}}}*/

int cmp_N64(uint64_t const * dst, uint64_t const * src, size_t m)/*{{{*/
{
    return memcmp(dst, src, m * sizeof(uint64_t));
}/*}}}*/

///////////////////////////////////////////////////////////////////////

/* implements mul_N64_6464 */
/* This can work in place (C==A, or C==B, or both) */
void mul_N64_6464_lookup4(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        uint64_t const * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    /* We don't zero out C before the computation, but rather at the
     * moment we read A[i], so that A==C is supported */
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 15]; aa>>=4;
        C[i]^= Bx[1][aa & 15]; aa>>=4;
        C[i]^= Bx[2][aa & 15]; aa>>=4;
        C[i]^= Bx[3][aa & 15]; aa>>=4;
        C[i]^= Bx[4][aa & 15]; aa>>=4;
        C[i]^= Bx[5][aa & 15]; aa>>=4;
        C[i]^= Bx[6][aa & 15]; aa>>=4;
        C[i]^= Bx[7][aa & 15]; aa>>=4;
        C[i]^= Bx[8][aa & 15]; aa>>=4;
        C[i]^= Bx[9][aa & 15]; aa>>=4;
        C[i]^= Bx[10][aa & 15]; aa>>=4;
        C[i]^= Bx[11][aa & 15]; aa>>=4;
        C[i]^= Bx[12][aa & 15]; aa>>=4;
        C[i]^= Bx[13][aa & 15]; aa>>=4;
        C[i]^= Bx[14][aa & 15]; aa>>=4;
        C[i]^= Bx[15][aa];
    }
}/*}}}*/
void mul_N64_6464_lookup4_blocks(mat64 *C,/*{{{*/
                   mat64 const *A,
                   mat64 const & B,
                   size_t nblocks,
                   size_t Cstride, size_t Astride)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        uint64_t const * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    /* We don't zero out C before the computation, but rather at the
     * moment we read A[i], so that A==C is supported */
    for (size_t b = 0 ; b < nblocks ; b++) {
        mat64 const & AA = A[b * Astride];
        mat64 & CC = C[b * Cstride];
        for (size_t i = 0 ; i < 64; i++) {
            uint64_t aa = AA[i];
            CC[i] = Bx[0][aa & 15]; aa>>=4;
            CC[i]^= Bx[1][aa & 15]; aa>>=4;
            CC[i]^= Bx[2][aa & 15]; aa>>=4;
            CC[i]^= Bx[3][aa & 15]; aa>>=4;
            CC[i]^= Bx[4][aa & 15]; aa>>=4;
            CC[i]^= Bx[5][aa & 15]; aa>>=4;
            CC[i]^= Bx[6][aa & 15]; aa>>=4;
            CC[i]^= Bx[7][aa & 15]; aa>>=4;
            CC[i]^= Bx[8][aa & 15]; aa>>=4;
            CC[i]^= Bx[9][aa & 15]; aa>>=4;
            CC[i]^= Bx[10][aa & 15]; aa>>=4;
            CC[i]^= Bx[11][aa & 15]; aa>>=4;
            CC[i]^= Bx[12][aa & 15]; aa>>=4;
            CC[i]^= Bx[13][aa & 15]; aa>>=4;
            CC[i]^= Bx[14][aa & 15]; aa>>=4;
            CC[i]^= Bx[15][aa];
        }
    }
}/*}}}*/
void addmul_N64_6464_lookup4_blocks(mat64 *C,/*{{{*/
                   mat64 const *A,
                   mat64 const & B,
                   size_t nblocks,
                   size_t Cstride, size_t Astride)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        uint64_t const * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    /* We don't zero out C before the computation, but rather at the
     * moment we read A[i], so that A==C is supported */
    for (size_t b = 0 ; b < nblocks ; b++) {
        mat64 const & AA = A[b * Astride];
        mat64 & CC = C[b * Cstride];
        for (size_t i = 0 ; i < 64; i++) {
            uint64_t aa = AA[i];
            CC[i]^= Bx[0][aa & 15]; aa>>=4;
            CC[i]^= Bx[1][aa & 15]; aa>>=4;
            CC[i]^= Bx[2][aa & 15]; aa>>=4;
            CC[i]^= Bx[3][aa & 15]; aa>>=4;
            CC[i]^= Bx[4][aa & 15]; aa>>=4;
            CC[i]^= Bx[5][aa & 15]; aa>>=4;
            CC[i]^= Bx[6][aa & 15]; aa>>=4;
            CC[i]^= Bx[7][aa & 15]; aa>>=4;
            CC[i]^= Bx[8][aa & 15]; aa>>=4;
            CC[i]^= Bx[9][aa & 15]; aa>>=4;
            CC[i]^= Bx[10][aa & 15]; aa>>=4;
            CC[i]^= Bx[11][aa & 15]; aa>>=4;
            CC[i]^= Bx[12][aa & 15]; aa>>=4;
            CC[i]^= Bx[13][aa & 15]; aa>>=4;
            CC[i]^= Bx[14][aa & 15]; aa>>=4;
            CC[i]^= Bx[15][aa];
        }
    }
}/*}}}*/

void mul_6464lt_6464_lookup4(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B)
{
    constexpr const size_t m = 64;
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        uint64_t const * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    /* We don't zero out C before the computation, but rather at the
     * moment we read A[i], so that A==C is supported */
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = 0;
        for(size_t j = 0 ; 4*j <= i ; j++) {
            C[i] ^= Bx[j][aa & 15]; aa>>=4;
        }
    }
}/*}}}*/

/* implements mul_N64_6464 */
void mul_N64_6464_lookup8(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    uint64_t Bx[8][256];
    for(int j = 0 ; j < 8 ; j++) {
        uint64_t const * bb = B.data() + 8 * j;
        uint64_t w = 0;
        Bx[j][0] = w; w ^= bb[0];
        Bx[j][1] = w; w ^= bb[1];
        Bx[j][3] = w; w ^= bb[0];
        Bx[j][2] = w; w ^= bb[2];
        Bx[j][6] = w; w ^= bb[0];
        Bx[j][7] = w; w ^= bb[1];
        Bx[j][5] = w; w ^= bb[0];
        Bx[j][4] = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9] = w; w ^= bb[0];
        Bx[j][8] = w; w ^= bb[4];
        Bx[j][24] = w; w ^= bb[0];
        Bx[j][25] = w; w ^= bb[1];
        Bx[j][27] = w; w ^= bb[0];
        Bx[j][26] = w; w ^= bb[2];
        Bx[j][30] = w; w ^= bb[0];
        Bx[j][31] = w; w ^= bb[1];
        Bx[j][29] = w; w ^= bb[0];
        Bx[j][28] = w; w ^= bb[3];
        Bx[j][20] = w; w ^= bb[0];
        Bx[j][21] = w; w ^= bb[1];
        Bx[j][23] = w; w ^= bb[0];
        Bx[j][22] = w; w ^= bb[2];
        Bx[j][18] = w; w ^= bb[0];
        Bx[j][19] = w; w ^= bb[1];
        Bx[j][17] = w; w ^= bb[0];
        Bx[j][16] = w; w ^= bb[5];
        Bx[j][48] = w; w ^= bb[0];
        Bx[j][49] = w; w ^= bb[1];
        Bx[j][51] = w; w ^= bb[0];
        Bx[j][50] = w; w ^= bb[2];
        Bx[j][54] = w; w ^= bb[0];
        Bx[j][55] = w; w ^= bb[1];
        Bx[j][53] = w; w ^= bb[0];
        Bx[j][52] = w; w ^= bb[3];
        Bx[j][60] = w; w ^= bb[0];
        Bx[j][61] = w; w ^= bb[1];
        Bx[j][63] = w; w ^= bb[0];
        Bx[j][62] = w; w ^= bb[2];
        Bx[j][58] = w; w ^= bb[0];
        Bx[j][59] = w; w ^= bb[1];
        Bx[j][57] = w; w ^= bb[0];
        Bx[j][56] = w; w ^= bb[4];
        Bx[j][40] = w; w ^= bb[0];
        Bx[j][41] = w; w ^= bb[1];
        Bx[j][43] = w; w ^= bb[0];
        Bx[j][42] = w; w ^= bb[2];
        Bx[j][46] = w; w ^= bb[0];
        Bx[j][47] = w; w ^= bb[1];
        Bx[j][45] = w; w ^= bb[0];
        Bx[j][44] = w; w ^= bb[3];
        Bx[j][36] = w; w ^= bb[0];
        Bx[j][37] = w; w ^= bb[1];
        Bx[j][39] = w; w ^= bb[0];
        Bx[j][38] = w; w ^= bb[2];
        Bx[j][34] = w; w ^= bb[0];
        Bx[j][35] = w; w ^= bb[1];
        Bx[j][33] = w; w ^= bb[0];
        Bx[j][32] = w; w ^= bb[6];
        Bx[j][96] = w; w ^= bb[0];
        Bx[j][97] = w; w ^= bb[1];
        Bx[j][99] = w; w ^= bb[0];
        Bx[j][98] = w; w ^= bb[2];
        Bx[j][102] = w; w ^= bb[0];
        Bx[j][103] = w; w ^= bb[1];
        Bx[j][101] = w; w ^= bb[0];
        Bx[j][100] = w; w ^= bb[3];
        Bx[j][108] = w; w ^= bb[0];
        Bx[j][109] = w; w ^= bb[1];
        Bx[j][111] = w; w ^= bb[0];
        Bx[j][110] = w; w ^= bb[2];
        Bx[j][106] = w; w ^= bb[0];
        Bx[j][107] = w; w ^= bb[1];
        Bx[j][105] = w; w ^= bb[0];
        Bx[j][104] = w; w ^= bb[4];
        Bx[j][120] = w; w ^= bb[0];
        Bx[j][121] = w; w ^= bb[1];
        Bx[j][123] = w; w ^= bb[0];
        Bx[j][122] = w; w ^= bb[2];
        Bx[j][126] = w; w ^= bb[0];
        Bx[j][127] = w; w ^= bb[1];
        Bx[j][125] = w; w ^= bb[0];
        Bx[j][124] = w; w ^= bb[3];
        Bx[j][116] = w; w ^= bb[0];
        Bx[j][117] = w; w ^= bb[1];
        Bx[j][119] = w; w ^= bb[0];
        Bx[j][118] = w; w ^= bb[2];
        Bx[j][114] = w; w ^= bb[0];
        Bx[j][115] = w; w ^= bb[1];
        Bx[j][113] = w; w ^= bb[0];
        Bx[j][112] = w; w ^= bb[5];
        Bx[j][80] = w; w ^= bb[0];
        Bx[j][81] = w; w ^= bb[1];
        Bx[j][83] = w; w ^= bb[0];
        Bx[j][82] = w; w ^= bb[2];
        Bx[j][86] = w; w ^= bb[0];
        Bx[j][87] = w; w ^= bb[1];
        Bx[j][85] = w; w ^= bb[0];
        Bx[j][84] = w; w ^= bb[3];
        Bx[j][92] = w; w ^= bb[0];
        Bx[j][93] = w; w ^= bb[1];
        Bx[j][95] = w; w ^= bb[0];
        Bx[j][94] = w; w ^= bb[2];
        Bx[j][90] = w; w ^= bb[0];
        Bx[j][91] = w; w ^= bb[1];
        Bx[j][89] = w; w ^= bb[0];
        Bx[j][88] = w; w ^= bb[4];
        Bx[j][72] = w; w ^= bb[0];
        Bx[j][73] = w; w ^= bb[1];
        Bx[j][75] = w; w ^= bb[0];
        Bx[j][74] = w; w ^= bb[2];
        Bx[j][78] = w; w ^= bb[0];
        Bx[j][79] = w; w ^= bb[1];
        Bx[j][77] = w; w ^= bb[0];
        Bx[j][76] = w; w ^= bb[3];
        Bx[j][68] = w; w ^= bb[0];
        Bx[j][69] = w; w ^= bb[1];
        Bx[j][71] = w; w ^= bb[0];
        Bx[j][70] = w; w ^= bb[2];
        Bx[j][66] = w; w ^= bb[0];
        Bx[j][67] = w; w ^= bb[1];
        Bx[j][65] = w; w ^= bb[0];
        Bx[j][64] = w; w ^= bb[7];
        Bx[j][192] = w; w ^= bb[0];
        Bx[j][193] = w; w ^= bb[1];
        Bx[j][195] = w; w ^= bb[0];
        Bx[j][194] = w; w ^= bb[2];
        Bx[j][198] = w; w ^= bb[0];
        Bx[j][199] = w; w ^= bb[1];
        Bx[j][197] = w; w ^= bb[0];
        Bx[j][196] = w; w ^= bb[3];
        Bx[j][204] = w; w ^= bb[0];
        Bx[j][205] = w; w ^= bb[1];
        Bx[j][207] = w; w ^= bb[0];
        Bx[j][206] = w; w ^= bb[2];
        Bx[j][202] = w; w ^= bb[0];
        Bx[j][203] = w; w ^= bb[1];
        Bx[j][201] = w; w ^= bb[0];
        Bx[j][200] = w; w ^= bb[4];
        Bx[j][216] = w; w ^= bb[0];
        Bx[j][217] = w; w ^= bb[1];
        Bx[j][219] = w; w ^= bb[0];
        Bx[j][218] = w; w ^= bb[2];
        Bx[j][222] = w; w ^= bb[0];
        Bx[j][223] = w; w ^= bb[1];
        Bx[j][221] = w; w ^= bb[0];
        Bx[j][220] = w; w ^= bb[3];
        Bx[j][212] = w; w ^= bb[0];
        Bx[j][213] = w; w ^= bb[1];
        Bx[j][215] = w; w ^= bb[0];
        Bx[j][214] = w; w ^= bb[2];
        Bx[j][210] = w; w ^= bb[0];
        Bx[j][211] = w; w ^= bb[1];
        Bx[j][209] = w; w ^= bb[0];
        Bx[j][208] = w; w ^= bb[5];
        Bx[j][240] = w; w ^= bb[0];
        Bx[j][241] = w; w ^= bb[1];
        Bx[j][243] = w; w ^= bb[0];
        Bx[j][242] = w; w ^= bb[2];
        Bx[j][246] = w; w ^= bb[0];
        Bx[j][247] = w; w ^= bb[1];
        Bx[j][245] = w; w ^= bb[0];
        Bx[j][244] = w; w ^= bb[3];
        Bx[j][252] = w; w ^= bb[0];
        Bx[j][253] = w; w ^= bb[1];
        Bx[j][255] = w; w ^= bb[0];
        Bx[j][254] = w; w ^= bb[2];
        Bx[j][250] = w; w ^= bb[0];
        Bx[j][251] = w; w ^= bb[1];
        Bx[j][249] = w; w ^= bb[0];
        Bx[j][248] = w; w ^= bb[4];
        Bx[j][232] = w; w ^= bb[0];
        Bx[j][233] = w; w ^= bb[1];
        Bx[j][235] = w; w ^= bb[0];
        Bx[j][234] = w; w ^= bb[2];
        Bx[j][238] = w; w ^= bb[0];
        Bx[j][239] = w; w ^= bb[1];
        Bx[j][237] = w; w ^= bb[0];
        Bx[j][236] = w; w ^= bb[3];
        Bx[j][228] = w; w ^= bb[0];
        Bx[j][229] = w; w ^= bb[1];
        Bx[j][231] = w; w ^= bb[0];
        Bx[j][230] = w; w ^= bb[2];
        Bx[j][226] = w; w ^= bb[0];
        Bx[j][227] = w; w ^= bb[1];
        Bx[j][225] = w; w ^= bb[0];
        Bx[j][224] = w; w ^= bb[6];
        Bx[j][160] = w; w ^= bb[0];
        Bx[j][161] = w; w ^= bb[1];
        Bx[j][163] = w; w ^= bb[0];
        Bx[j][162] = w; w ^= bb[2];
        Bx[j][166] = w; w ^= bb[0];
        Bx[j][167] = w; w ^= bb[1];
        Bx[j][165] = w; w ^= bb[0];
        Bx[j][164] = w; w ^= bb[3];
        Bx[j][172] = w; w ^= bb[0];
        Bx[j][173] = w; w ^= bb[1];
        Bx[j][175] = w; w ^= bb[0];
        Bx[j][174] = w; w ^= bb[2];
        Bx[j][170] = w; w ^= bb[0];
        Bx[j][171] = w; w ^= bb[1];
        Bx[j][169] = w; w ^= bb[0];
        Bx[j][168] = w; w ^= bb[4];
        Bx[j][184] = w; w ^= bb[0];
        Bx[j][185] = w; w ^= bb[1];
        Bx[j][187] = w; w ^= bb[0];
        Bx[j][186] = w; w ^= bb[2];
        Bx[j][190] = w; w ^= bb[0];
        Bx[j][191] = w; w ^= bb[1];
        Bx[j][189] = w; w ^= bb[0];
        Bx[j][188] = w; w ^= bb[3];
        Bx[j][180] = w; w ^= bb[0];
        Bx[j][181] = w; w ^= bb[1];
        Bx[j][183] = w; w ^= bb[0];
        Bx[j][182] = w; w ^= bb[2];
        Bx[j][178] = w; w ^= bb[0];
        Bx[j][179] = w; w ^= bb[1];
        Bx[j][177] = w; w ^= bb[0];
        Bx[j][176] = w; w ^= bb[5];
        Bx[j][144] = w; w ^= bb[0];
        Bx[j][145] = w; w ^= bb[1];
        Bx[j][147] = w; w ^= bb[0];
        Bx[j][146] = w; w ^= bb[2];
        Bx[j][150] = w; w ^= bb[0];
        Bx[j][151] = w; w ^= bb[1];
        Bx[j][149] = w; w ^= bb[0];
        Bx[j][148] = w; w ^= bb[3];
        Bx[j][156] = w; w ^= bb[0];
        Bx[j][157] = w; w ^= bb[1];
        Bx[j][159] = w; w ^= bb[0];
        Bx[j][158] = w; w ^= bb[2];
        Bx[j][154] = w; w ^= bb[0];
        Bx[j][155] = w; w ^= bb[1];
        Bx[j][153] = w; w ^= bb[0];
        Bx[j][152] = w; w ^= bb[4];
        Bx[j][136] = w; w ^= bb[0];
        Bx[j][137] = w; w ^= bb[1];
        Bx[j][139] = w; w ^= bb[0];
        Bx[j][138] = w; w ^= bb[2];
        Bx[j][142] = w; w ^= bb[0];
        Bx[j][143] = w; w ^= bb[1];
        Bx[j][141] = w; w ^= bb[0];
        Bx[j][140] = w; w ^= bb[3];
        Bx[j][132] = w; w ^= bb[0];
        Bx[j][133] = w; w ^= bb[1];
        Bx[j][135] = w; w ^= bb[0];
        Bx[j][134] = w; w ^= bb[2];
        Bx[j][130] = w; w ^= bb[0];
        Bx[j][131] = w; w ^= bb[1];
        Bx[j][129] = w; w ^= bb[0];
        Bx[j][128] = w;
    }
    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 255]; aa>>=8;
        C[i]^= Bx[1][aa & 255]; aa>>=8;
        C[i]^= Bx[2][aa & 255]; aa>>=8;
        C[i]^= Bx[3][aa & 255]; aa>>=8;
        C[i]^= Bx[4][aa & 255]; aa>>=8;
        C[i]^= Bx[5][aa & 255]; aa>>=8;
        C[i]^= Bx[6][aa & 255]; aa>>=8;
        C[i]^= Bx[7][aa & 255];
    }
}/*}}}*/

/* implements mul_N64_6464 */
void mul_N64_6464_vec(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        mul_o64_6464(C++, *A++, B);
    }
}/*}}}*/

/* implements mul_N64_6464 */
void mul_N64_6464_transB(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    mat64 tb;
    mat64_transpose(tb, B);
    mul_N64_T6464(C, A, tb, m);
}/*}}}*/

#if defined(HAVE_AVX2)
/* implements mul_N64_6464 */
void mul_N64_6464_avx2(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
    /* can work in place, so not simply memset0 + addmul (the ^= have been
     * changed to =)
     */
    size_t j;
    constexpr const unsigned int SIMD = 4;
    __m256i *Cw = (__m256i *) C;
    __m256i *Aw = (__m256i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    __m256d zd = _mm256_setzero_pd();

    for (j = 0; j + (SIMD - 1) < m ; j += SIMD*2) {
        __m256i c0 = _mm256_setzero_si256();
        __m256i c1 = _mm256_setzero_si256();
	__m256i a0 = *Aw++;
	__m256i a1 = *Aw++;

	for (int i = 64; i--;) {
            __m256d Bd = _mm256_castsi256_pd(_mm256_set1_epi64x(B[i]));
            __m256d c0d = _mm256_blendv_pd(
                    zd,
                    Bd,
                    _mm256_castsi256_pd(a0));
            __m256d c1d = _mm256_blendv_pd(
                    zd,
                    Bd,
                    _mm256_castsi256_pd(a1));
            c0 = _mm256_xor_si256(c0, _mm256_castpd_si256(c0d));
            c1 = _mm256_xor_si256(c1, _mm256_castpd_si256(c1d));
            a0 = _mm256_slli_epi64(a0, 1);
            a1 = _mm256_slli_epi64(a1, 1);
	}
	*Cw++ = c0;
	*Cw++ = c1;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ = c;
    }
}/*}}}*/
#endif

#if defined(HAVE_SSE2) && ULONG_BITS == 64
/* implements mul_N64_6464 */
void mul_N64_6464_sse(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
    /* can work in place, so not simply memset0 + addmul (the ^= have been
     * changed to =)
     */
    size_t j;
    constexpr const unsigned int SIMD = 2;
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    for (j = 0; j + (SIMD - 1) < m ; j += SIMD) {
        __m128i c = _mm_setzero_si128();
	__m128i a = *Aw++;

        __m128i one = _cado_mm_set1_epi64_c(1);
	for (int i = 0; i < 64; i++) {
	    __m128i bw = _cado_mm_set1_epi64(B[i]);
            // c ^= (bw & -(a & one));
            c = _mm_xor_si128(c, _mm_and_si128(bw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
	    a = _mm_srli_epi64(a, 1);
	}
	*Cw++ = c;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ = c;
    }
}/*}}}*/
#endif

///////////////////////////////////////////////////////////////////////

/* implements mul_N64_T6464 */
void mul_N64_T6464_vec(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        mul_o64_T6464(C++, *A++, B);
    }
}
/*}}}*/

/* implements mul_N64_T6464 */
void mul_N64_T6464_transB(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    mat64 tb;
    mat64_transpose(tb, B);
    mul_N64_6464(C, A, tb, m);
}/*}}}*/

///////////////////////////////////////////////////////////////////////

/* implements addmul_N64_6464 */
/* This can work in place (C==A, or C==B, or both) */
void MAYBE_UNUSED addmul_N64_6464_lookup4(uint64_t *C, /*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        uint64_t const * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i]^= Bx[0][aa & 15]; aa>>=4;
        C[i]^= Bx[1][aa & 15]; aa>>=4;
        C[i]^= Bx[2][aa & 15]; aa>>=4;
        C[i]^= Bx[3][aa & 15]; aa>>=4;
        C[i]^= Bx[4][aa & 15]; aa>>=4;
        C[i]^= Bx[5][aa & 15]; aa>>=4;
        C[i]^= Bx[6][aa & 15]; aa>>=4;
        C[i]^= Bx[7][aa & 15]; aa>>=4;
        C[i]^= Bx[8][aa & 15]; aa>>=4;
        C[i]^= Bx[9][aa & 15]; aa>>=4;
        C[i]^= Bx[10][aa & 15]; aa>>=4;
        C[i]^= Bx[11][aa & 15]; aa>>=4;
        C[i]^= Bx[12][aa & 15]; aa>>=4;
        C[i]^= Bx[13][aa & 15]; aa>>=4;
        C[i]^= Bx[14][aa & 15]; aa>>=4;
        C[i]^= Bx[15][aa];
    }
}/*}}}*/

#if defined(HAVE_AVX2)
/* implements addmul_N64_6464 */
/* C == A seems to work ok */
void addmul_N64_6464_avx2(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
    size_t j;
    __m256i *Cw = (__m256i *) C;
    __m256i *Aw = (__m256i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    constexpr const unsigned int SIMD = 4;
    for (j = 0; j + (SIMD - 1) < m ; j += SIMD) {
        __m256i c = _mm256_setzero_si256();
	__m256i a = *Aw++;

        __m256i one = _mm256_set1_epi64x(INT64_C(1));
	for (int i = 0; i < 64; i++) {
	    __m256i bw = _mm256_set1_epi64x(B[i]);
	    // c ^= (bw & -(a & one));
            c = _mm256_xor_si256(c, _mm256_and_si256(bw, _mm256_sub_epi64(_mm256_setzero_si256(), _mm256_and_si256(a, one))));
	    a = _mm256_srli_epi64(a, 1);
	}
	*Cw = _mm256_xor_si256(*Cw, c);
        Cw++;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ ^= c;
    }
}/*}}}*/
#endif
#if defined(HAVE_SSE2) && ULONG_BITS == 64
/* implements addmul_N64_6464 */
/* C == A seems to work ok */
void addmul_N64_6464_sse(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
    size_t j;
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    constexpr const unsigned int SIMD = 2;
    for (j = 0; j + (SIMD - 1) < m ; j += SIMD) {
        __m128i c = _mm_setzero_si128();
	__m128i a = *Aw++;

        __m128i one = _cado_mm_set1_epi64_c(1);
	for (int i = 0; i < 64; i++) {
	    __m128i bw = _cado_mm_set1_epi64(B[i]);
	    // c ^= (bw & -(a & one));
            c = _mm_xor_si128(c, _mm_and_si128(bw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
	    a = _mm_srli_epi64(a, 1);
	}
	*Cw = _mm_xor_si128(*Cw, c);
        Cw++;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ ^= c;
    }
}/*}}}*/
#endif

///////////////////////////////////////////////////////////////////////

/* implements mul_TN64_N64 */
void mul_TN64_N64_addmul(mat64 & r, uint64_t const *a, uint64_t const *w, size_t n)/*{{{*/
{
    r = 0;
    for (size_t i = 0; i < n; i++) {
        addmul_To64_o64(r, a[i], w[i]);
    }
}/*}}}*/

/* implements mul_TN64_N64 */
void mul_TN64_N64_C(mat64 & b, uint64_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    b = 0;
    addmul_TN64_N64_C(b, A, x, ncol);
}/*}}}*/

///////////////////////////////////////////////////////////////////////

/* implements addmul_TN64_N64 */
void addmul_TN64_N64_C(mat64 & b, uint64_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    uint64_t idx, i, rA;
    uint64_t rx;

    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 64; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}/*}}}*/

///////////////////////////////////////////////////////////////////////

/* implements mul_TN32_N64 */
void mul_TN32_N64_C(uint64_t * b, uint32_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    uint32_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 32 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 32; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}/*}}}*/

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_TN64K_N64_sse2(uint64_t * w, uint64_t const * u, uint64_t const * v, unsigned int n, unsigned int K)/*{{{*/
{
    for(unsigned int i = 0 ; i < n ; i++) {
        __m128i * w0 = (__m128i*) w;
        // TODO: It's possible to expand more, and use a __m128i
        // mb[4][2], or even [4]. This wouldn't change the code much
        // (see the m128 version), and is likely to speed things up a
        // wee bit maybe.
        __m128i mb[4] = {
            _mm_setzero_si128(),
            _cado_mm_setr_epi64(*v, 0 ),
            _cado_mm_setr_epi64(0,  *v),
            _cado_mm_set1_epi64(*v),
        };
        v++;
        __m128i *sw = w0;
        for(unsigned int k = 0 ; k < K ; k++) {
            uint64_t a = *u++;
            for (unsigned int j = 0; j < 64; j += 2) {
                *sw = _mm_xor_si128(*sw, mb[a & 3]);
                a >>= 2;
                sw ++;
            }
        }
    }
}/*}}}*/
#endif

void addmul_TN64K_N64_C(uint64_t * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K)/*{{{*/
{
    uint64_t idx, i, rA;
    uint64_t rx;

    for(idx = 0; idx < ncol; idx++) {
        rx = x[idx];
        uint64_t * pb = (uint64_t *) b;
        for(unsigned int j = 0 ; j < K ; j++) {
            rA = *A++;
            for(i = 0; i < 64; i++) {
                *pb++ ^= rx & -(rA & 1);
                rA >>= 1;
            }
        }
    }
}/*}}}*/

///////////////////////////////////////////////////////////////////////


/* {{{ final choices. These are static choices at this point, but it should
 * be the result of some tuning, ideally */
void mul_N64_6464(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
/* The chosen function is optimal (among the ones here) for N about
 * 20000. At N=2000000, a twice faster version can be obtained. However,
 * it's not critical for cado, so we stick with the slower version.
 *
 * TODO: the lack of proper tuning is now becoming slightly problematic.
 * lookup code seems to be consistently faster, in fact, but lookup8
 * takes over at some point. And it's a bit weird that no sse/avx2
 * version seems to win at this time.
 */
#if defined(HAVE_AVX2)
    mul_N64_6464_avx2(C,A,B,m);
#elif defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_N64_6464_sse(C,A,B,m);
#else
    mul_N64_6464_lookup4(C,A,B,m);
#endif
}/*}}}*/
void addmul_6464_blocks(mat64 *C,
                   mat64 const *A,
                   mat64 const & B, size_t nblocks, size_t Cstride, size_t Astride)
{
    addmul_N64_6464_lookup4_blocks(C, A, B, nblocks, Cstride, Astride);
}

void mul_6464_blocks(mat64 *C,
                   mat64 const *A,
                   mat64 const & B, size_t nblocks, size_t Cstride, size_t Astride)
{
    mul_N64_6464_lookup4_blocks(C, A, B, nblocks, Cstride, Astride);
}

void addmul_N64_6464(uint64_t *C,/*{{{*/
		 uint64_t const *A,
		 mat64 const & B, size_t m)
{
#if defined(HAVE_AVX2)
    addmul_N64_6464_avx2(C,A,B,m);
#elif defined(HAVE_SSE2) && ULONG_BITS == 64
    addmul_N64_6464_sse(C,A,B,m);
#else
    addmul_N64_6464_lookup4(C,A,B,m);
#endif
}/*}}}*/
void mul_N64_T6464(uint64_t *C,/*{{{*/
                   uint64_t const *A,
                   mat64 const & B, size_t m)
{
    mul_N64_T6464_transB(C,A,B,m);
}/*}}}*/
void mul_TN64_N64(mat64 & b, uint64_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    mul_TN64_N64_C(b, A, x, ncol);
}/*}}}*/
void addmul_TN64_N64(mat64 & b, uint64_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    addmul_TN64_N64_C(b, A, x, ncol);
}/*}}}*/
void mul_TN32_N64(uint64_t * b, uint32_t const * A, uint64_t const * x, unsigned int ncol)/*{{{*/
{
    mul_TN32_N64_C(b, A, x, ncol);
}
/*}}}*/
void addmul_TN64K_N64(uint64_t * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K)/*{{{*/
{
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    addmul_TN64K_N64_sse2(b, A, x, ncol, K);
#else
    addmul_TN64K_N64_C(b, A, x, ncol, K);
#endif
}/*}}}*/

void mul_TN64K_N64(uint64_t * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K)/*{{{*/
{
    ASSERT_ALWAYS(b != A && b != x);
    memset((void *) b, 0, 64 * K * sizeof(uint64_t));
    addmul_TN64K_N64(b, A, x, ncol, K);
}/*}}}*/

/*}}}*/



