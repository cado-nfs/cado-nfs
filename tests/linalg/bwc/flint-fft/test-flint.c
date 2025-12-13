#include "cado.h" // IWYU pragma: keep

#include <string.h>
#include <stdlib.h>

#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include <gmp.h>

#include "gmp_aux.h"

#include "flint-fft/flint.h"
#include "flint-fft/transform_interface.h"
#include "portability.h"
#include "gmp-hacks.h"
#include "macros.h"

void get_ft_hash(mpz_t h, int bits_per_coeff, void * data, struct fft_transform_info * fti);

#define xxPARI
/*{{{ display macros */
#ifdef DEBUG_FFT
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
#else
/* remain silent if we're not working for the external checking script.
 */
#define pint(_name, _x, _nx) /**/
#define ppol(_name, _x, _cx) /**/
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
    // fprintf(stderr, "bits_of_p:=%zu;\n", bits_of_p);

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
        // fprintf(stderr, "cx:=%d; cy:=%d;\n", *cx, *cy);
    } while ((*cx+*cy-1) * np * FLINT_BITS < 4000);

    // finer-grain control can go here. But it does not change the
    // picture much.
    if (*cx > *cy) { int a; a = *cx; *cx = *cy; *cy = a; }

    fprintf(stderr, "s:=%d;\n", s);
#ifdef PARI
    gmp_printf("p=%Zd;\n", p);
#else
    gmp_printf("p:=%Zd;\n", p);
    printf("KP<x>:=PolynomialRing(GF(p));\n");
#endif

    return 1;
}
/*}}}*/

void fti_disp(FILE * f, struct fft_transform_info* fti)
{
    fprintf(f, "fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, (size_t) fti->depth);
    fprintf(f, "fti_trunc0:=%lu;\n", fti->trunc0);
    fprintf(f, "fti_w:=%lu;\n", fti->w);
    fprintf(f, "fti_alg:=%d;\n", fti->alg);
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
mp_size_t nz0;
mp_limb_t * x;
mp_limb_t * y;
mp_limb_t * z;
mp_limb_t * z0;
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
    if (nz0) z0 = malloc(nz0 * sizeof(mp_limb_t));
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
    if (nz0) free(z0);
}

static void prepare_transforms() {
    tx = malloc(fft_alloc_sizes[0]);
    ty = malloc(fft_alloc_sizes[0]);
    tz = malloc(fft_alloc_sizes[0]);
    qt = malloc(fft_alloc_sizes[1]);
    tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_prepare(fti, tx);
    fft_prepare(fti, ty);
    fft_prepare(fti, tz);
}

static void do_renames(const struct fft_transform_info * fti MAYBE_UNUSED, const char * step MAYBE_UNUSED, const char * varname MAYBE_UNUSED)
{
#ifdef DEBUG_FFT
    char * s, * t;
    int rc;
    rc = asprintf(&s, "%s/%s_before_%s.m", fti->tmpdir, varname, step);
    if (rc < 0) abort();
    rc = asprintf(&t, "%s/before_%s.m", fti->tmpdir, step);
    if (rc < 0) abort();
    rc = rename(t, s);
    if (rc < 0) {
        fprintf(stderr, "rename(%s,%s) : %s\n", t, s, strerror(errno));
        exit(EXIT_FAILURE);
    }
    free(t);
    free(s);
    rc = asprintf(&s, "%s/%s_after_%s.m", fti->tmpdir, varname, step);
    if (rc < 0) abort();
    rc = asprintf(&t, "%s/after_%s.m", fti->tmpdir, step);
    if (rc < 0) abort();
    rc = rename(t, s);
    if (rc < 0) {
        fprintf(stderr, "rename(%s,%s) : %s\n", t, s, strerror(errno));
        exit(EXIT_FAILURE);
    }
    free(t);
    free(s);
#else
    // fprintf(stderr, "// no data file for %s (%s)\n", step, varname);
#endif
}

/* test multiplication of integers */
int test_mul0(gmp_randstate_t rstate) /*{{{*/
{
    int xbits, ybits;
    int nacc = 4;
    operand_sizes(&xbits, &ybits, s, rstate);

    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);
    int zbits = xbits + ybits;

    fft_transform_info_init(fti, xbits, ybits, nacc);
    fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
    fti_disp(stdout, fti);
    fti_disp(stderr, fti);

    nx = iceildiv(xbits, FLINT_BITS);
    ny = iceildiv(ybits, FLINT_BITS);
    nz = iceildiv(zbits, FLINT_BITS);
    nz0 = nx + ny;

    alloc_everything();

    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);

    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    pint("A0", x, nx);
    pint("A1", y, ny);
    fft_dft(fti, tx, x, nx, tt); do_renames(fti, "dft", "P0");
    fft_dft(fti, ty, y, ny, tt); do_renames(fti, "dft", "P1");
    fft_compose(fti, tz, tx, ty, tt);
    fft_ift(fti, z, nz, tz, tt); do_renames(fti, "ift", "P2");
    pint("A2", z, nz);

    if (nx >= ny)
        mpn_mul(z0, x, nx, y, ny);
    else
        mpn_mul(z0, y, ny, x, nx);

    if (memcmp(z, z0, nz * sizeof(mp_limb_t)) != 0) {
        fprintf(stderr, "test_mul0 failed!!!\n");
        abort();
    }

    free_everything();
    return 0;
}/*}}}*/
int test_mul(gmp_randstate_t rstate) /*{{{*/
{
    /* compute 2xy+x+y */

    int xbits, ybits;
    int nacc = 4;
    operand_sizes(&xbits, &ybits, s, rstate);

    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);
    int zbits = xbits + ybits;

    fft_transform_info_init(fti, xbits, ybits, nacc);
    fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
    fti_disp(stdout, fti);
    fti_disp(stderr, fti);

    nx = iceildiv(xbits, FLINT_BITS);
    ny = iceildiv(ybits, FLINT_BITS);
    nz = iceildiv(zbits, FLINT_BITS);
    nz0 = 0;

    alloc_everything();

    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);

    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    pint("A0", x, nx);
    pint("A1", y, ny);
    fft_dft(fti, tx, x, nx, tt); do_renames(fti, "dft", "P0");
    fft_dft(fti, ty, y, ny, tt); do_renames(fti, "dft", "P1");
    fft_compose(fti, tz, tx, ty, tt);
    fft_add(fti, tz, tz, tz);
    fft_add(fti, tz, tz, tx);
    fft_add(fti, tz, tz, ty);
    fft_ift(fti, z, nz, tz, tt); do_renames(fti, "ift", "P2");
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

    fft_transform_info_init_mulmod(fti, xbits, ybits, nacc, minwrap);

    fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
    fti_disp(stdout, fti);
    fti_disp(stderr, fti);

    nx = iceildiv(xbits, FLINT_BITS);
    ny = iceildiv(ybits, FLINT_BITS);
    nz = fft_get_mulmod_output_minlimbs(fti);
    nz0 = 0;

    alloc_everything();
    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    pint("A0", x, nx);
    pint("A1", y, ny);
    fft_dft(fti, tx, x, nx, tt); do_renames(fti, "dft", "P0");
    fft_dft(fti, ty, y, ny, tt); do_renames(fti, "dft", "P1");
    fft_compose(fti, tz, tx, ty, tt);
    fft_add(fti, tz, tz, tz);
    fft_add(fti, tz, tz, tx);
    fft_add(fti, tz, tz, ty);
    fft_ift(fti, z, nz, tz, tt); do_renames(fti, "ift", "P2");
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

    fft_transform_info_init_fppol(fti, p, cx, cy, 5);
    fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
    fti_disp(stdout, fti);
    fti_disp(stderr, fti);

    nx = cx * np; ny = cy * np; nz = cz * np;
    nz0 = 0;

    alloc_everything();
    bitrandom_fppol(y, cy, p, longstrings, rstate);
    bitrandom_fppol(x, cx, p, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    ppol("P0", x, cx);
    ppol("P1", y, cy);
    fft_dft(fti, tx, x, cx, tt); do_renames(fti, "dft", "P0");
    fft_dft(fti, ty, y, cy, tt); do_renames(fti, "dft", "P1");
    fft_compose(fti, tz, tx, ty, tt);
    fft_add(fti, tz, tz, tz);
    /* P0 is always smaller than P1, so the product P0^2 fits within the
     * requested transform length */
    fft_addcompose(fti, tz, tx, tx, tt, qt);
    fft_add(fti, tz, tz, tx);
    fft_add(fti, tz, tz, ty);
    fft_ift(fti, z, cz, tz, tt); do_renames(fti, "ift", "P2");
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

    fft_transform_info_init_fppol_mp(fti, p, cx, cy, 4);
    fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
    fti_disp(stdout, fti);
    fti_disp(stderr, fti);

    nx = cx * np; ny = cy * np; nz = cz * np;
    nz0 = 0;
    
    alloc_everything();
    bitrandom_fppol(y, cy, p, longstrings, rstate);
    bitrandom_fppol(x, cx, p, longstrings, rstate);
    prepare_transforms();

    printf("check:=\"%s\";\n", __func__);
    ppol("P0", x, cx);
    ppol("P1", y, cy);
    fft_dft(fti, tx, x, cx, tt); do_renames(fti, "dft", "P0");
    fft_dft(fti, ty, y, cy, tt); do_renames(fti, "dft", "P1");
    fft_compose(fti, tz, tx, ty, tt);
    fft_ift(fti, z, cz, tz, tt); do_renames(fti, "ift", "P2");
    ppol("P2", z, cz);

    free_everything();
    mpz_clear(p);
    return 0;
}/*}}}*/

int main(int argc, char const * argv[])
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


    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

#ifdef PARI
    printf("allocatemem(800000000)\n");
#endif

    int do_mul0 = 0;
    int do_mul = 0;
    int do_mulmod = 0;
    int do_mul_fppol = 0;
    int do_mp_fppol = 0;

    int done_tests = 0;
    for(int i = 1 ; i < argc ; i++) {
        const char * p = argv[i];
        if (strncmp(p, "s=", 2) == 0) {
            s = atoi(p + 2);
            continue;
        }
        if (i + 1 < argc && strcmp(p, "-s") == 0) {
            i++; p=argv[i];
            s = atoi(p);
            continue;
        }
        if (strncmp(p, "seed=", 5) == 0) {
            seed = atol(p + 5);
            continue;
        }
        if (i + 1 < argc && strcmp(p, "-seed") == 0) {
            i++; p=argv[i];
            seed = atoi(p);
            continue;
        }
        if (strncmp(p, "test_", 5) == 0) p += 5;
        if (strcmp(p, "mul0") == 0) { do_mul0++; continue; }
        if (strcmp(p, "mul") == 0) { do_mul++; continue; }
        if (strcmp(p, "mulmod") == 0) { do_mulmod++; continue; }
        if (strcmp(p, "mul_fppol") == 0) { do_mul_fppol++; continue; }
        if (strcmp(p, "mp_fppol") == 0) { do_mp_fppol++; continue; }
        fprintf(stderr, "Unexpected argument: %s\n", p);
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "seed:=%d; longstrings:=%d; */\n", seed, longstrings);
    gmp_randseed_ui(rstate, seed);

    if (do_mul0) { test_mul0(rstate); done_tests++; }
    if (do_mul) { test_mul(rstate); done_tests++; }
    if (do_mulmod) {
        abort();
        /* this test fails, currently. The bug is easy to reproduce by
         * running this test over and over again. The faulty code is
         * within fft_transform_info_adjust_depth, where we "refine bits"
         * -- this might kill minwrap most probably the whole
         * fft_transform_info_adjust_depth should be scrutinized.
         */
        test_mulmod(rstate);
        done_tests++;
    }
    if (do_mul_fppol) { test_mul_fppol(rstate); done_tests++; }
    if (do_mp_fppol) { test_mp_fppol(rstate); done_tests++; }

    if (!done_tests) {
        fprintf(stderr, "Please supply at least one test name\n");
        exit(EXIT_FAILURE);
    }
    gmp_randclear(rstate);
    return 0;
}

/* The magma script in test-flint.m may be run *after* the output of this
 * program */
/* Note that DEBUG_FFT is necessary for this to make any sense.
 */
