#include "cado.h" // IWYU pragma: keep
#include <algorithm>                      // for min
#include "macros.h"                       // for ASSERT
#include "bblas_bitmat.hpp"  // for bitmat_ops, bblas_bitmat_de...
#include "bblas_mat8.hpp"
#include "bblas_bitmat_inl.hpp" // IWYU pragma: keep

using namespace bblas_bitmat_details;

template<>
void bitmat_ops<uint8_t>::add(mat8 & C, mat8 const & A, mat8 const & B)
{
    uint64_t & Cx = * (uint64_t *) C.data();
    uint64_t const & Ax = * (uint64_t const *) A.data();
    uint64_t const & Bx = * (uint64_t const *) B.data();
    Cx = Ax ^ Bx;
}

template<>
void bitmat_ops<uint8_t>::transpose(mat8 & C, mat8 const & A)
{
    uint64_t aa = * (uint64_t const *) A.data();

    uint64_t t;
    t = (aa ^ (aa >> 28));
    t = t & UINT64_C(0x00000000F0F0F0F0);
    aa ^= t ^ (t << 28);
    t = (aa ^ (aa >> 14));
    t = t & UINT64_C(0x0000cccc0000cccc);
    aa ^= t ^ (t << 14);
    t = (aa ^ (aa >> 7));
    t = t & UINT64_C(0x00aa00aa00aa00aa);
    aa ^= t ^ (t << 7);
    * (uint64_t *) C.data() = aa;
}

template<>
void bitmat_ops<uint8_t>::mul(mat8 & C, mat8 const & A, mat8 const & B)
{
    mat8 AB = 0;
    addmul(AB, A, B);
    C = AB;
}

template<>
void bitmat_ops<uint8_t>::addmul_blocks(mat8 * C, mat8 const * A, mat8 const& B, size_t nblocks, size_t Cstride, size_t Astride)
{
    for(size_t i = 0 ; i < nblocks ; i++) {
        addmul(C[i * Cstride], A[i * Astride], B);
    }
}
template<>
void bitmat_ops<uint8_t>::mul_blocks(mat8 * C, mat8 const * A, mat8 const& B, size_t nblocks, size_t Cstride, size_t Astride)
{
    for(size_t i = 0 ; i < nblocks ; i++) {
        mul(C[i * Cstride], A[i * Astride], B);
    }
}

void addmul8_naive(mat8 & C,
        mat8 const & A,
        mat8 const & B,
        unsigned int i0,
        unsigned int i1,
        unsigned int yi0,
        unsigned int yi1)
{
    // const uint64_t ones = UINT64_C(0x0101010101010101);
    // uint64_t keepmask = (ones << yi1) - (ones << yi0);
    // uint64_t aa = keepmask & * (uint64_t const *) A.data();
    uint8_t Bx[2][16];
    unsigned int j0 = yi0 / 4;
    unsigned int j1 = (yi1 + 3) / 4;
    for(unsigned int j = j0 ; j < j1 ; j++) {
        const uint8_t * bb = B.data() + 4 * j;
        uint8_t w = 0;
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
    uint8_t mask = (1 << yi1) - (1 << yi0);
    if (yi1 == 8)
        mask = - (1 << yi0);
    for (size_t i = i0; i < i1; i++) {
        uint8_t aa = (A[i] & mask) >> (4 * j0);
        for(unsigned int j = j0 ; j < j1 ; j++) {
            C[i]^= Bx[j][aa & 15]; aa>>=4;
        }
    }
}

template<>
void bitmat_ops<uint8_t>::addmul(mat8 & C,
        mat8 const & A,
        mat8 const & B,
        unsigned int i0,
        unsigned int i1,
        unsigned int yi0,
        unsigned int yi1)
{
    addmul8_naive(C, A, B, i0, i1, yi0, yi1);
}

void trsm8_naive(mat8 const & L,
        mat8 & U,
        unsigned int n0,
        unsigned int n1)
{
    ASSERT(n0 <= n1);
    if (n1 <= n0 + 1) return;
    /* need to determine the very first fragment before we can align */
    if (n0 % 4) {
        uint8_t c[8];
        unsigned int n0b = std::min(n0 + 4 - (n0 % 4), n1);
        uint8_t * uu = U.data() + n0;
        uint8_t const * ll= L.data() + n0;
        unsigned int d = n0b - n0;
        c[0] = 0;
        c[1] = uu[0];
        uint8_t m = 1;
        if (d >= 2) {
            c[2]=uu[1]^=c[(ll[1]>>n0)&1];
            c[3]=uu[1]^uu[0];
            m=3;
        }
        if (d >= 3) {
            c[4]=uu[2]^=c[(ll[2]>>n0)&3];
            c[5]=uu[2]^c[1];
            c[6]=uu[2]^c[2];
            c[7]=uu[2]^c[3];
            m=7;
        }
        for(unsigned int i = n0b ; i < n1 ; i++)
            U[i] ^= c[(L[i] >> n0)&m];
        n0 = n0b;
    }
    if (n1 <= n0 + 1) return;
    unsigned int b = n0;
    for(b = n0 ; b < n1 - (n1 % 4) ; b += 4) {
        uint8_t c[16];
        uint8_t * uu = U.data() + b;
        uint8_t const * ll = L.data() + b;
        c[0]=0;
        c[1]=uu[0];
        c[2]=uu[1]^=c[(ll[1]>>b)&1];
        c[3]=uu[1]^uu[0];
        c[4]=uu[2]^=c[(ll[2]>>b)&3];
        c[5]=uu[2]^c[1];
        c[6]=uu[2]^c[2];
        c[7]=uu[2]^c[3];
        c[8]=uu[3]^=c[(ll[3]>>b)&7];
        /* we might break here if b == 60 */
        c[9]=uu[3]^c[1];
        c[10]=uu[3]^c[2];
        c[11]=uu[3]^c[3];
        c[12]=uu[3]^c[4];
        c[13]=uu[3]^c[5];
        c[14]=uu[3]^c[6];
        c[15]=uu[3]^c[7];
        for(unsigned int i = b + 4 ; i < n1 ; i++)
            U[i] ^= c[(L[i] >> b)&15];
    }
    if (n1 % 4) {
        ASSERT(b < n1);
        unsigned int d = n1 % 4;
        ASSERT(n1 == b + d);
        if (d <= 1) return;
        uint8_t c[4];
        uint8_t * uu = U.data() + b;
        uint8_t const * ll = L.data() + b;
        c[0]=0;
        c[1]=uu[0];
        if (d >= 2) {
            c[2]=uu[1]^=c[(ll[1]>>b)&1];
            c[3]=uu[1]^uu[0];
        }
        if (d >= 3) {
                 uu[2]^=c[(ll[2]>>b)&3];
        }
    }
}

template<>
void bitmat_ops<uint8_t>::trsm(mat8 const & L,
        mat8 & U,
        unsigned int n0,
        unsigned int n1)
{
    trsm8_naive(L, U, n0, n1);
}

void mat8_add_C(mat8 & C, mat8 const & A, mat8 const & B)/*{{{*/
{
    for(unsigned int j = 0 ; j < mat8::width ; j++) {
        C[j] = A[j] ^ B[j];
    }
}
/*}}}*/

template struct bblas_bitmat_details::bitmat_ops<uint8_t>;
