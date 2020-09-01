#ifndef BALANCING_H_
#define BALANCING_H_

#include <stdint.h>
#include "macros.h"  // for MIN
#include "mod_ul.h"

// #include "mf.h"

/* This is used only for sanity checking after dispatch. */
/* The primes are 2^64-59 and 2^63-25 */
#define DUMMY_VECTOR_COORD_VALUE(j)     \
                (UINT64_C(0xffffffffffffffc5) / (1 + (j))  \
              + UINT64_C(0x7fffffffffffffe7) * (1 + (j)))
#define DUMMY_VECTOR_COORD_VALUE2(j)     \
        (DUMMY_VECTOR_COORD_VALUE((j) ^ 0xbeef) ^ ((DUMMY_VECTOR_COORD_VALUE((j) ^ 0xcafe) << 32)))

/* here's some magma code for generating these arbitrary vectors, but
 * it's not checked.
  nrows:=65536;
  vv:=func<x|VectorSpace(GF(2),64)!Intseq(x mod 2^64,2,64)>;
  p:=2^64-59;
  q:=2^63-25;
  i2p:=func<a|Polynomial(GF(2),Intseq(a,2))>;
  p2i:=func<p|Seqint(ChangeUniverse(Eltseq(p),Integers()),2)>;
  ixor:=func<a,b|p2i(i2p(a)+i2p(b))>;
  vs:=func<v|vv(Seqint(ChangeUniverse(Eltseq(v),Integers()),2)*2^32)>;
  arbitrary1:=func<j|vv(p div j + q * j)>;
  arbitrary2:=[arbitrary1(1+ixor(j, 0xbeef)) + vs(arbitrary1(1+ixor(j, 0xcafe))):j in [0..nrows-1]];
 */

/* balancing structure */

#define FLAG_COLPERM    1
#define FLAG_ROWPERM    2
#define FLAG_REPLICATE  8       /* work with square matrices (padding if
                                   need be), and replicate one
                                   permutation to the other side, so that
                                   we get conjugated permutations.
                                   */

#define BALANCING_MAGIC UINT32_C(0xba1a0000)

struct balancing_header_s {
    uint32_t zero;      /* previous versions had neither zero nor magic.
                           By enforcing a zero field here, we make sure
                           that older and newer code will choke on
                           incompatible balancing files */
    uint32_t magic;
    uint32_t nh;
    uint32_t nv;
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nzrows;
    uint32_t nzcols;
    uint64_t ncoeffs;
    uint32_t checksum;
    uint32_t flags;
    // pshuf indicates two integers a,b such that the COLUMN i of the input
    // matrix is in fact mapped to column a*i+b mod n in the matrix we work
    // with. pshuf_inv indicates the inverse permutation. a and b do
    // **NOT** depend on nh and nv. Therefore all balancing are
    // compatible even though a permutation is being used.
    // a == rhuf[0] is the smallest integer greater than floor(sqrt(n)) which
    // is coprime to n. b is 42. pshuf_inv[0,1] are computed accordingly.
    //
    // This means that on all occasions, we are in fact working with the
    // matrix M*S, for some permutation S. As long as we are concerned
    // with left nullspace, this does not make any difference, but it
    // does for right nullspace !
    /* Note: pshuf and pshuf_inv have to do with what is called
     * "shuffled-product" elsewhere (and has actually become the default)
     */
    uint32_t pshuf[2];
    uint32_t pshuf_inv[2];
};
typedef struct balancing_header_s balancing_header[1];
typedef struct balancing_header_s * balancing_header_ptr;

struct balancing_s {
    balancing_header h;
    uint32_t trows;     // number of rows including padding.
    uint32_t tcols;     // number of cols including padding.
    uint32_t * rowperm; // row index for new mat. --> row index for old mat.
    uint32_t * colperm; // might be equal to colperm.
};
typedef struct balancing_s balancing[1];
typedef struct balancing_s * balancing_ptr;

#ifdef __cplusplus
extern "C" {
#endif

/* Once the flags and perm[] fields have been provided, the caller must
 * call _finalize() in order to 1) update the trows and tcols fields 2)
 * compute the checksum of the balancing.
 */
extern void balancing_set_row_col_count(balancing_ptr bal);
extern void balancing_finalize(balancing_ptr bal);
extern void balancing_write_inner(balancing_ptr bal, const char *);
extern void balancing_write(balancing_ptr bal, const char * , const char *);
extern void balancing_read(balancing_ptr bal, const char *);
extern void balancing_read_header(balancing_ptr bal, const char * filename);
extern void balancing_clear(balancing_ptr bal);
extern void balancing_init(balancing_ptr bal);

/* helper for the functions below */
static inline unsigned long balancing_index_shuffle_common_(unsigned long r, unsigned long n, uint32_t * shuf)
{
    modulusul_t M;
    modul_initmod_ul(M, n);
    residueul_t x,a,b;
    modul_init_noset0(a, M);
    modul_init_noset0(b, M);
    modul_init_noset0(x, M);
    modul_set_ul(a, shuf[0], M);
    modul_set_ul(b, shuf[1], M);
    modul_set_ul(x, r, M);
    modul_mul(x, x, a, M);
    modul_add(x, x, b, M);
    unsigned long r1 = modul_get_ul(x, M);
    modul_clear(a, M);
    modul_clear(b, M);
    modul_clear(x, M);
    modul_clearmod(M);
    return r1;
}

/* These two relate to the global permutation represented by the rshuf /
 * rshuf_inv arrays */
static inline unsigned long balancing_pre_shuffle(balancing_ptr bal, unsigned long r)
{
    unsigned int K = MIN(bal->h->ncols, bal->h->nrows);
    if (r >= K) return r;
    return balancing_index_shuffle_common_(r, K, bal->h->pshuf);
}
static inline unsigned long balancing_pre_unshuffle(balancing_ptr bal, unsigned long r)
{
    unsigned int K = MIN(bal->h->ncols, bal->h->nrows);
    if (r >= K) return r;
    return balancing_index_shuffle_common_(r, K, bal->h->pshuf_inv);
}


#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_H_ */
