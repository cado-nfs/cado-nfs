#include "cado.h" // IWYU pragma: keep
#define _GNU_SOURCE
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include "gmp_aux.h"    // mpn_rrandom
#include "flint-fft/fft.h"
#include "portability.h"

/*{{{ macros */
#ifndef PTR
#define PTR(x) ((x)->_mp_d)
#endif
#ifndef SIZ
#define SIZ(x) ((x)->_mp_size)
#endif
#ifndef ALLOC
#define ALLOC(x) ((x)->_mp_alloc)
#endif

#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {                                                                  \
    while ((NLIMBS) > 0)                                                \
      {                                                                 \
        if ((DST)[(NLIMBS) - 1] != 0)                                   \
          break;                                                        \
        (NLIMBS)--;                                                     \
      }                                                                 \
  } while (0)
#endif

#define BEGIN_BLOCK	do {
#define END_BLOCK	} while (0)

#define MPZ_GROW_ALLOC(DST, NLIMBS)					\
	BEGIN_BLOCK							\
	if (ALLOC(DST) < (int) (NLIMBS)) {				\
		ALLOC(DST) = (NLIMBS);					\
		PTR(DST)=(mp_limb_t *) realloc(PTR(DST),		\
				(NLIMBS) * sizeof(mp_limb_t));		\
	}								\
	END_BLOCK

#define MPZ_INIT_SET_MPN(DST, SRC, NLIMBS)				\
	BEGIN_BLOCK							\
	ALLOC(DST) = (NLIMBS);						\
	SIZ(DST) = (NLIMBS);						\
	PTR(DST) = (mp_limb_t *)malloc((NLIMBS) * sizeof(mp_limb_t));	\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	
#define MPZ_SET_MPN(DST, SRC, NLIMBS)					\
	BEGIN_BLOCK							\
	MPZ_GROW_ALLOC(DST, NLIMBS);					\
	SIZ(DST) = (NLIMBS);						\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	   
#define MPN_SET_MPZ(DST,NLIMBS,SRC)					\
	BEGIN_BLOCK							\
        mp_size_t r = MIN(ABS(SIZ(SRC)),NLIMBS);                        \
	memcpy((DST),PTR(SRC),r * sizeof(mp_limb_t));   		\
	memset((DST)+SIZ(SRC),0,((NLIMBS)-r) * sizeof(mp_limb_t));      \
	END_BLOCK


#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#ifndef iceildiv
/* unfortunately this fails miserably if x+y-1 overflows */
#define iceildiv(x,y)	(((x)+(y)-1)/(y))
#endif
/*}}}*/

void get_ft_hash(mpz_t h, int bits_per_coeff, void * data, struct fft_transform_info * fti);

#define xxPARI
/*{{{ display macros */
#ifdef PARI
#define ppol(_name, _x, _cx) do {					\
    mpz_t tmp;								\
    mpz_init(tmp);							\
    printf("v"_name"=[");						\
    for(int i = 0 ; i < _nx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + (_nx-1-i) * np, np);			\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		        	\
    }									\
    printf("];\n");							\
    printf(_name"=Pol(vector(#v"_name",i,Mod(v"_name"[i],p)));\n");	\
    mpz_clear(tmp);							\
} while (0)
#define pint(_name, _x, _nx) gmp_printf(_name "=%Nd;\n", _x, _nx)
#else
#define ppol(_name, _x, _cx) do {					\
    mpz_t tmp;								\
    mpz_init(tmp);							\
    printf(_name ":=Polynomial(GF(p),[");				\
    for(int i = 0 ; i < _cx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + i * np, np);				\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		                \
    }									\
    printf("]);\n");							\
    mpz_clear(tmp);							\
} while (0)
#define pint(_name, _x, _cx) gmp_printf(_name ":=%Nd;\n", _x, _cx)
#endif
/*}}}*/

/*{{{ setup operand sizes */
int operand_sizes(int * xbits, int * ybits, int s, gmp_randstate_t rstate)
{
    int base;
    int rs;

    /* We can't handle too small examples */
    do {
        rs = 10 + gmp_urandomm_ui(rstate, 200);
        if (s == 0) s = rs;
        base = 120 * s + gmp_urandomm_ui(rstate, 20 * s);
        *xbits = base + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        *ybits = base + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        if (*xbits < 10) *xbits = 10;
        if (*ybits < 10) *ybits = 10;
        /* Do this so that we don't loop forever on easy cases */
        s++;
    } while (*xbits + *ybits < 4000);
    if (*xbits > *ybits) { int a; a = *xbits; *xbits = *ybits; *ybits = a; }
    return 1;
}

int operand_sizes_fppol(int * cx, int * cy, mpz_t p, int s, gmp_randstate_t rstate)
{
    size_t bits_of_p;
    int n;
    int rs;
    bits_of_p = 32 + gmp_urandomm_ui(rstate, 512);

    mp_size_t np;

    do {
        rs = 10 + gmp_urandomm_ui(rstate, 200);
        if (s == 0) s = rs;
        mpz_ui_pow_ui(p, 2, bits_of_p);
        mpz_sub_ui(p, p, 1);
        for( ; !mpz_probab_prime_p(p, 2) ; mpz_sub_ui(p, p, 2));
        np = mpz_size(p);
        n = 20 * s + gmp_urandomm_ui(rstate, 10 * s);
        *cx = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        *cy = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        if (*cx < 10) *cx = 10;
        if (*cy < 10) *cy = 10;
    } while ((*cx+*cy-1) * np * FLINT_BITS < 4000);

    // finer-grain control can go here. But it does not change the
    // picture much.
    if (*cx > *cy) { int a; a = *cx; *cx = *cy; *cy = a; }

#ifdef PARI
    gmp_printf("p=%Zd;\n", p);
#else
    gmp_printf("p:=%Zd;\n", p);
    printf("KP<x>:=PolynomialRing(GF(p));\n");
#endif

    return 1;
}
/*}}}*/

void fti_disp(struct fft_transform_info* fti)
{
    printf("fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, (size_t) fti->depth);
    printf("fti_trunc0:=%lu;\n", fti->trunc0);
    printf("fti_w:=%lu;\n", fti->w);
    printf("fti_alg:=%d;\n", fti->alg);
}

/*{{{ setting entries to random */
void bitrandom(mp_limb_t * x, int xbits, int longstrings, gmp_randstate_t rstate)
{
    int nx = iceildiv(xbits, FLINT_BITS);
    if (longstrings) {
        mpn_rrandom(x, rstate, nx);
    } else {
        mpn_randomb(x, rstate, nx);
    }
    if (xbits % FLINT_BITS) { x[nx-1] &= (1UL<<(xbits%FLINT_BITS))-1; }
}

void bitrandom_fppol(mp_limb_t * x, int cx, mpz_srcptr p, int longstrings, gmp_randstate_t rstate)
{
    mpz_t tmp;
    mpz_init(tmp);
    int np = mpz_size(p);
    int xbits = cx * mpz_size(p) * FLINT_BITS;
    if (longstrings) {
        mpn_rrandom(x, rstate, cx * np);
    } else {
        mpn_randomb(x, rstate, cx * np);
    }
    if (xbits % FLINT_BITS) { x[cx-1] &= (1UL<<(xbits%FLINT_BITS))-1; }
    for(int i = 0 ; i < cx ; i++) {
        MPZ_SET_MPN(tmp, x + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(x + i * np, np, tmp);
    }
    mpz_clear(tmp);
}
/*}}}*/


/* These are globals, used within this script */

int seed;
int s = 0;
int longstrings = 0;
int xbits;
int ybits;
mp_size_t nx;
mp_size_t ny;
mp_size_t nz;
mp_limb_t * x;
mp_limb_t * y;
mp_limb_t * z;
struct fft_transform_info fti[1];
size_t fft_alloc_sizes[3];
void * tx;
void * ty;
void * tz;
void * tt;
void * qt;

static void alloc_everything()
{
    x = malloc(nx * sizeof(mp_limb_t));
    y = malloc(ny * sizeof(mp_limb_t));
    z = malloc(nz * sizeof(mp_limb_t));
}

static void free_everything() {
    free(tx);
    free(ty);
    free(tz);
    free(tt);
    free(qt);
    free(x);
    free(y);
    free(z);
}

static void prepare_transforms() {
    tx = malloc(fft_alloc_sizes[0]);
    ty = malloc(fft_alloc_sizes[0]);
    tz = malloc(fft_alloc_sizes[0]);
    qt = malloc(fft_alloc_sizes[1]);
    tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_transform_prepare(tx, fti);
    fft_transform_prepare(ty, fti);
    fft_transform_prepare(tz, fti);
}

static void do_renames(const char * step, const char * varname)
{
#ifdef DEBUG_FFT
    char * s, * t;
    int rc;
    rc = asprintf(&s, "/tmp/%s_before_%s.m", varname, step);
    if (rc < 0) abort();
    rc = asprintf(&t, "/tmp/before_%s.m", step);
    if (rc < 0) abort();
    rename(t, s);
    free(t);
    free(s);
    rc = asprintf(&s, "/tmp/%s_after_%s.m", varname, step);
    if (rc < 0) abort();
    rc = asprintf(&t, "/tmp/after_%s.m", step);
    if (rc < 0) abort();
    rename(t, s);
    free(t);
    free(s);
#else
    fprintf(stderr, "// no data file for %s (%s)\n", step, varname);
#endif
}

/* test multiplication of integers */
int test_mul(gmp_randstate_t rstate) /*{{{*/
{
    int xbits, ybits;
    int nacc = 4;
    operand_sizes(&xbits, &ybits, s, rstate);

    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);
    int zbits = xbits + ybits;

    fft_get_transform_info(fti, xbits, ybits, nacc);
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    fti_disp(fti);

    nx = iceildiv(xbits, FLINT_BITS);
    ny = iceildiv(ybits, FLINT_BITS);
    nz = iceildiv(zbits, FLINT_BITS);

    alloc_everything();
    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    pint("A0", x, nx);
    pint("A1", y, ny);
    fft_do_dft(tx, x, nx, tt, fti); do_renames("dft", "P0");
    fft_do_dft(ty, y, ny, tt, fti); do_renames("dft", "P1");
    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    fft_do_ift(z, nz, tz, tt, fti); do_renames("ift", "P2");
    pint("A2", z, nz);

    free_everything();
    return 0;
}/*}}}*/

/* test wrapped product of integers. This computea A*B mod base^n\pm1 */
int test_mulmod(gmp_randstate_t rstate) /*{{{*/
{
    int xbits, ybits;
    int nacc = 4;
    operand_sizes(&xbits, &ybits, s, rstate);
    int wrap = gmp_urandomm_ui(rstate, 128);

    int minwrap = ybits + wrap;

    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);

    fft_get_transform_info_mulmod(fti, xbits, ybits, nacc, minwrap);

    fft_get_transform_allocs(fft_alloc_sizes, fti);
    fti_disp(fti);

    nx = iceildiv(xbits, FLINT_BITS);
    ny = iceildiv(ybits, FLINT_BITS);
    nz = fft_get_mulmod_output_minlimbs(fti);

    alloc_everything();
    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    pint("A0", x, nx);
    pint("A1", y, ny);
    fft_do_dft(tx, x, nx, tt, fti); do_renames("dft", "P0");
    fft_do_dft(ty, y, ny, tt, fti); do_renames("dft", "P1");
    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    fft_do_ift(z, nz, tz, tt, fti); do_renames("ift", "P2");
    pint("A2", z, nz);

    free_everything();
    return 0;
}/*}}}*/

/* multiplication of polynomials */
int test_mul_fppol(gmp_randstate_t rstate) /*{{{*/
{
    mpz_t p;
    mpz_init(p);
    int cx, cy, cz;

    operand_sizes_fppol(&cx, &cy, p, s, rstate);
    mp_size_t np = mpz_size(p);

    fprintf(stderr, "cx:=%d; cy:=%d; bits_of_p:=%zu;\n", cx, cy, mpz_sizeinbase(p, 2));
    cz = cx + cy - 1;

    fft_get_transform_info_fppol(fti, p, cx, cy, 5);
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    fti_disp(fti);

    nx = cx * np; ny = cy * np; nz = cz * np;

    alloc_everything();
    bitrandom_fppol(y, cy, p, longstrings, rstate);
    bitrandom_fppol(x, cx, p, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    ppol("P0", x, cx);
    ppol("P1", y, cy);
    fft_do_dft_fppol(tx, x, cx, tt, fti, p); do_renames("dft", "P0");
    fft_do_dft_fppol(ty, y, cy, tt, fti, p); do_renames("dft", "P1");
    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    /* P0 is always smaller than P1, so the product P0^2 fits within the
     * requested transform length */
    fft_addmul(tz, tx, tx, tt, qt, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    fft_do_ift_fppol(z, cz, tz, tt, fti, p); do_renames("ift", "P2");
    /* beware: after the IFT, coefficient of indices >= trunc are not
     * computed at all -- there's noise in there ! */
    ppol("P2", z, cz);

    free_everything();
    mpz_clear(p);
    return 0;
}/*}}}*/

/* middle product of polynomials */
int test_mp_fppol(gmp_randstate_t rstate)/*{{{*/
{
    mpz_t p;
    mpz_init(p);
    int cx, cy, cz;

    operand_sizes_fppol(&cx, &cy, p, s, rstate);
    mp_size_t np = mpz_size(p);

    /* We're doing the transpose of
     * MUL(cx, cy) == cz ; which is MP(cx, cz) == cy.
     * But we rewrite this as MP(cx, cy) == cz by swapping cy and cz.
     */
    cz = cy;
    cy = cx + cz - 1;
    assert(cy >= cx);
    fprintf(stderr, "/* MP(degree %d, degree %d) -> degree %d */\n",
            cx - 1, cy - 1, cz - 1);

    fft_get_transform_info_fppol_mp(fti, p, cx, cy, 4);
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    fti_disp(fti);

    nx = cx * np; ny = cy * np; nz = cz * np;
    
    alloc_everything();
    bitrandom_fppol(y, cy, p, longstrings, rstate);
    bitrandom_fppol(x, cx, p, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    ppol("P0", x, cx);
    ppol("P1", y, cy);
    fft_do_dft_fppol(tx, x, cx, tt, fti, p); do_renames("dft", "P0");
    fft_do_dft_fppol(ty, y, cy, tt, fti, p); do_renames("dft", "P1");
    fft_mul(tz, tx, ty, tt, fti);
    fft_do_ift_fppol_mp(z, cz, tz, tt, fti, p, cx - 1); do_renames("ift", "P2");
    ppol("P2", z, cz);

    free_everything();
    mpz_clear(p);
    return 0;
}/*}}}*/

int main(int argc, char * argv[])
{
    seed = getpid();
    // s=84; seed=12682; longstrings=0; 
    // s=93; seed=13156; longstrings=0;
    // s=10; seed=22442; longstrings=0;
    // s=84; seed=12682; longstrings=0; 
    // s=93; seed=13156; longstrings=0;
    // s=10; seed=22442; longstrings=0;
    // s=0; seed=8412; longstrings=0;
    // seed=6286; longstrings=0;
    // s=0; seed=16083; longstrings=0;
    // s=0; seed=19066; longstrings=0;
    // s=0; seed=19239; longstrings=0;
    // s=0; seed=19302; longstrings=0;
    // s=0; seed=25058; longstrings=0;
    // s=12; seed=1010; longstrings=0;
    // s=24; seed=6931; longstrings=0;


    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);

#ifdef PARI
    printf("allocatemem(800000000)\n");
#endif

    int done_tests = 0;
    for(int i = 1 ; i < argc ; i++) {
        const char * p = argv[i];
        if (strncmp(p, "test_", 5) == 0) p += 5;
        if (strcmp(p, "mul") == 0) { test_mul(rstate); done_tests++; continue; }
        if (strcmp(p, "mulmod") == 0) { test_mulmod(rstate); done_tests++; continue; }
        if (strcmp(p, "mul_fppol") == 0) { test_mul_fppol(rstate); done_tests++; continue; }
        if (strcmp(p, "mp_fppol") == 0) { test_mp_fppol(rstate); done_tests++; continue; }
        fprintf(stderr, "Unexpected argument: %s\n", p);
        exit(EXIT_FAILURE);
    }
    if (!done_tests) {
        fprintf(stderr, "Please supply at least one test name\n");
        exit(EXIT_FAILURE);
    }
    gmp_randclear(rstate);
    return 0;
}

#if 0

/* This magma script may be run *after* the output of this program */
/* Note that DEBUG_FFT is necessary for this to make any sense.
 */

n:=2^fti_depth;
Z:=Integers();
R:=Integers(2^(n*fti_w) + 1);
RP<T>:=PolynomialRing(R);
/* sqrt(2) in R */
rho:=fti_w mod 2 eq 0 select
        R!2^(fti_w div 2)
    else
        R!(2^(3*u)-2^u)^fti_w
    where u is (n*fti_w div 4);

/* These make sense only in the context of the matrix algorithm (MFA) */
depth1 := fti_depth div 2;
depth2 := fti_depth + 1 - depth1;
n1 := 2^depth1; // for MFA

/* compute the truncation point */
tr:=fti_trunc0;
if tr le 2*n then tr:=2*n+1; end if;
if fti_alg eq 0 then
// trunc must be even and greater than 2n
tr:=fti_trunc0 + (fti_trunc0 mod 2);
else
// trunc must be greater than 2n and multiple of 2*n1
tr:= 2 * n1 * Ceiling(tr / (2 * n1));
end if;

/* transform coefficients are created in bitrev order, at least for the
 * non-MFA algorithm. This script is presently incapable of checking the
 * validity of transformed data for the MFA, but fixing it should not be
 * terribly hard. */
bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;

// mfaorder:=func<x|Seqint([b[1+i]:i in mfabitorder],2) where b is Intseq(x,2,fti_depth+2)>;


/* precompute powers */
seqmatch:=func<S,T,n|S[1..n] eq T[1..n]>;
pows:=[R|1];
for i in [1..2^(fti_depth+2)-1] do Append(~pows, pows[#pows]*rho); end for;
assert rho^(4*n) eq 1;

load "/tmp/P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
load "/tmp/P0_after_dft.m";  tQ0:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_after_dft.m";  tQ1:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P2_after_ift.m";  cQ2:=Polynomial([R!Seqint(x,2^64):x in data]);
cQ2 mod:= T^tr;

ch:=func<Q,tQ|fti_alg eq 0 select  seqmatch([Evaluate(Q,pows[i+1]):i in bitrevseq(fti_depth+2)], tQ, tr) else "skipped (MFA)">;

/* work around magma bugs */
if not assigned P0 then P0:=0; end if;
if not assigned P1 then P1:=0; end if;
if not assigned P2 then P2:=0; end if;
if not assigned A0 then A0:=0; end if;
if not assigned A1 then A1:=0; end if;
if not assigned A2 then A2:=0; end if;

print "Active test: ", check;

if check eq "test_mul" then
    Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
    Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0+Q1) mod (T^(4*n)-1);
    Q2 eq cQ2;
    A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits);
elif check eq "test_mulmod" then
    Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
    Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0+Q1) mod (T^(4*n)-1);
    Q2 eq cQ2;
    A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits) mod (2^(4*n*fti_bits)-1);
elif check eq "test_mul_fppol" then
    zP0:=Evaluate(ChangeRing(P0,Z),2^fti_ks_coeff_bits);
    zP1:=Evaluate(ChangeRing(P1,Z),2^fti_ks_coeff_bits);
    Q0 eq Polynomial(R,Intseq(zP0,2^fti_bits));
    Q1 eq Polynomial(R,Intseq(zP1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0^2+Q0+Q1) mod (T^(4*n)-1);
    Q2 eq cQ2;
    P2 eq 2*P0*P1+P0^2+P0+P1;
elif check eq "test_mp_fppol" then
    zP0:=Evaluate(ChangeRing(P0,Z),2^fti_ks_coeff_bits);
    zP1:=Evaluate(ChangeRing(P1,Z),2^fti_ks_coeff_bits);
    Q0 eq Polynomial(R,Intseq(zP0,2^fti_bits));
    Q1 eq Polynomial(R,Intseq(zP1,2^fti_bits));
    Q2:=(Q0*Q1) mod (T^(4*n)-1);
    Q2 eq cQ2;
    zP2:=Evaluate(ChangeRing(Q2,Z),2^fti_bits);
    P2:=Intseq(zP2,2^fti_ks_coeff_bits)[Degree(P0)+1..Degree(P1)+1];
    MP:=func<P0,P1|[Coefficient(P,i) :i in [Min(Degree(P0),Degree(P1))..Max(Degree(P0),Degree(P1))]] where P is P0*P1>;
    MP(P0,P1) eq P2;
else
    print "check not understood";
end if;

print "checking transform T0: ", ch(Q0, tQ0);
print "checking transform T1: ", ch(Q1, tQ1);
print "checking transform T2: ", ch(Q2, tQ2);


#endif
