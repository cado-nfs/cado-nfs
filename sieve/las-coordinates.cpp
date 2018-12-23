#include "cado.h"
#include "las-config.h"
#include "las-info.hpp"
#include "las-coordinates.hpp"
#include "gmp_aux.h"

/*  Conversions between different representations for sieve locations:
 *
 * Coordinate systems for identifying an (a,b) pair.
 *
 * (a, b): This is the one from textbooks. We always have b>=0 (only free
 *         relations have b=0, and we don't see them here).
 * (i, j): For a given special q, this is a point in the q-lattice. Given
 *         the lattice basis given by (a0 b0 a1 b1), this corresponds to
 *         the (a,b) pair equal to i*(a0,b0)+j*(a1,b1). By construction
 *         this should lead to one of the norms having doing.p as a factor.
 *         i is within [-I/2, I/2[, and j is within [1, J[
 * (N, x): bucket number N, location x. N is within [0,nb_buckets[
 *         and x within [0,bucket_region[ ; we have:
 *         N*bucket_region+x == (I/2+i)+j*I
 *
 * There are some change of coordinate functions at the end of this file.
 * The coordinate system (N, x) is almost always referred to as just x,
 * where x is a 64-bit value equal to N*bucket_region+x. Thus from this
 * wide x, it is possible to recover both N and (the short, fitting in 32
 * bits) x.
 * Note: the "wide" x might be larger than 32 bits only for I> 16.
 */
void adjustIJsublat(int & i, unsigned int & j, sublat_t const & S)
{
    if (S.m != 0) {
        i = i*S.m + S.i0;
        j = j*S.m + S.j0;
    }
}

void xToIJ(int & i, unsigned int & j, const uint64_t X, int logI)
{
    i = (X & ((1 << logI)-1)) - (1 << (logI - 1));
    j = X >> logI;
}


void NxToIJ(int & i, unsigned int & j, const unsigned int N, const unsigned int x, int logI)
{
    uint64_t X = (uint64_t)x + (((uint64_t)N) << LOG_BUCKET_REGION);
    xToIJ(i, j, X, logI);
}

void IJTox(uint64_t & x, int i, unsigned int j, int logI)
{
    x = (int64_t)i + (((uint64_t)j) << logI) + (uint64_t)(1 << (logI - 1));
}

void IJToNx(unsigned int & N, unsigned int & x, int i, unsigned int j, int logI)
{
    uint64_t xx;
    IJTox(xx, i, j, logI);
    N = xx >> LOG_BUCKET_REGION;
    x = xx & (uint64_t)((1 << LOG_BUCKET_REGION) - 1);
}

void IJToAB(int64_t & a, uint64_t & b, int i, unsigned int j, 
        qlattice_basis const & Q)
{
    adjustIJsublat(i, j, Q.sublat);

    int64_t s = (int64_t)i * (int64_t) Q.a0 + (int64_t)j * (int64_t) Q.a1;
    int64_t t = (int64_t)i * (int64_t) Q.b0 + (int64_t)j * (int64_t) Q.b1;
    if (t >= 0) {
        a = s;
        b = t;
    } else {
        a = -s;
        b = -t;
    }
}

int ABToIJ(int & i, unsigned int & j, const int64_t a, const uint64_t b, qlattice_basis const & Q)
{
    /* Both a,b and the coordinates of the lattice basis can be quite
     * large. However the result should be small.
     */
    mpz_t za,zb,ii,jj,a0,b0,a1,b1;
    mpz_init(ii);
    mpz_init(jj);
    mpz_init(za); mpz_init(zb);
    mpz_init(a0); mpz_init(b0);
    mpz_init(a1); mpz_init(b1);
    mpz_set_int64(za, a); mpz_set_uint64(zb, b);
    mpz_set_int64(a0, Q.a0);
    mpz_set_int64(b0, Q.b0);
    mpz_set_int64(a1, Q.a1);
    mpz_set_int64(b1, Q.b1);
    int ok = 1;
    mpz_mul(ii, za, b1); mpz_submul(ii, zb, a1);
    mpz_mul(jj, zb, a0); mpz_submul(jj, za, b0);
    /*
    int64_t ii =   a * (int64_t) Q.b1 - b * (int64_t)Q.a1;
    int64_t jj = - a * (int64_t) Q.b0 + b * (int64_t)Q.a0;
    */
    if (!mpz_divisible_p(ii, Q.doing.p)) ok = 0;
    if (!mpz_divisible_p(jj, Q.doing.p)) ok = 0;
    mpz_divexact(ii, ii, Q.doing.p);
    mpz_divexact(jj, jj, Q.doing.p);

    if (mpz_sgn(jj) < 0 || (mpz_sgn(jj) == 0 && mpz_sgn(ii) < 0)) {
        mpz_neg(ii, ii);
        mpz_neg(jj, jj);
    }
    i = mpz_get_si(ii);
    j = mpz_get_ui(jj);
    if (Q.sublat.m != 0) {
        int64_t imodm = i % int64_t(Q.sublat.m);
        if (imodm < 0) {
            imodm += Q.sublat.m;
        }
        int64_t jmodm = j % int64_t(Q.sublat.m);
        if (jmodm < 0)
            jmodm += Q.sublat.m;
        if (imodm != Q.sublat.i0 || jmodm != Q.sublat.j0) {
            fprintf(stderr, "# TraceAB: (i,j) does not belong to the right congruence class\n");
            ok = 0;
        } else {
            i = (i - imodm) / int64_t(Q.sublat.m);
            j = (j - jmodm) / int64_t(Q.sublat.m);
        }
    }

    mpz_clear(a0);
    mpz_clear(b0);
    mpz_clear(a1);
    mpz_clear(b1);
    mpz_clear(za);
    mpz_clear(zb);
    mpz_clear(ii);
    mpz_clear(jj);
    return ok;
}

#if 0 /* currently unused */
int ABTox(unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    if (!ABToIJ(i, j, a, b, Q)) return 0;
    IJTox(x, a, b, logI);
    return 1;
}
#endif

#if 0 /* currently unused */
int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    if (!ABToIJ(i, j, a, b, Q)) return 0;
    IJToNx(N, x, a, b, logI);
    return 1;
}
#endif

