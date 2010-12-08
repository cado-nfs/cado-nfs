#include "cado.h"
#define xxxWITH_BARRETT

#include <gmp.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "macros.h"
#include "powers_of_p.h"

double program_starttime;
double print_delay = 0.1;
#define WCT     (wct_seconds() - program_starttime)

struct prime_data {
    unsigned long p;
    // unsigned long * r;
    // unsigned long rj;
    void * powers;      // see .cpp file.

};/* }}} */

int degree;
mpz_t * f_hat;
mpz_t * f_hat_diff;

// {{{ trivial utility
static const char * size_disp(size_t s, char buf[16])
{
    char * prefixes = "bkMGT";
    double ds = s;
    const char * px = prefixes;
    for( ; px[1] && ds > 500.0 ; ) {
        ds /= 1024.0;
        px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}
// }}}

/* {{{ wrappers for some gmp operations, so as to report timings */
// above this threshold, we report each multiplication we do.
#define MUL_REPORT_THRESHOLD    1000000

#define REPORT_THIS(na, nb)     \
    (((na) > 10) && ((nb) > 10) && ((na) + (nb) > MUL_REPORT_THRESHOLD))

static void WRAP_mpz_mul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    double w0 = wct_seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_mul(c,a,b);
    double t1 = seconds();
    double w1 = wct_seconds();
    double rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_invert(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    double w0 = wct_seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_invert(c,a,b);
    double t1 = seconds();
    double w1 = wct_seconds();
    double rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
    if (REPORT_THIS(na, nb)) {
        printf("inv %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_addmul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    double w0 = wct_seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_addmul(c,a,b);
    double t1 = seconds();
    double w1 = wct_seconds();
    double rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_submul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    double w0 = wct_seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_submul(c,a,b);
    double t1 = seconds();
    double w1 = wct_seconds();
    double rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_barrett_mod(mpz_ptr c, mpz_srcptr a, mpz_srcptr p, mpz_srcptr q)
{
    double t0 = seconds();
    double w0 = wct_seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(p);
    barrett_mod(c,a,p,q);
    double t1 = seconds();
    double w1 = wct_seconds();
    double rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
    if (REPORT_THIS(na, nb) && na > nb + 10) {
        printf("mod %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
/* }}} */

/* {{{ evaluation of polynomials (one or two polynomials) */
static void mp_poly_eval_mod(mpz_ptr r, mpz_t * poly, int deg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
        WRAP_mpz_mul(r, r, a);
        // falls back on mpz_mod if qx == NULL
        WRAP_barrett_mod(r, r, q, qx);
        // mpz_mod(r, r, q);
        mpz_add(r, r, poly[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
}

static void mp_2poly_eval_mod(mpz_ptr r, mpz_ptr s, mpz_t * f, mpz_t * g, int degf, int degg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    if (r == NULL) {
        mp_poly_eval_mod(s, g, degg, a, q, qx);
        return;
    }

    if (s == NULL) {
        mp_poly_eval_mod(r, f, degf, a, q, qx);
        return;
    }

    mpz_t w;
    mpz_init(w);
    mpz_set(r,f[0]);
    mpz_set(s,g[0]);
    mpz_set(w, a);
    WRAP_mpz_addmul(r, w, f[1]);
    WRAP_mpz_addmul(s, w, g[1]);
    for(i = 2 ; i <= degf && i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    for( ; i <= degf ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
    }
    for( ; i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
    WRAP_barrett_mod(s, s, q, qx);
    mpz_clear(w);
}
/* }}} */

/* {{{ This wraps around barrett reduction, depending on whether we
 * choose to use it or not */
mpz_ptr my_barrett_init(mpz_srcptr px MAYBE_UNUSED)
{
#ifndef WITH_BARRETT
    return NULL;
#else
    if (mpz_size(px) < 10000)
        return NULL;
    mpz_ptr qx = malloc(sizeof(mpz_t));
    mpz_init(qx);
    barrett_init(qx, px);
    return qx;
#endif
}

void my_barrett_clear(mpz_ptr qx)
{
    if (!qx) return;
    mpz_clear(qx);
    free(qx);
}
/* }}} */

void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/* {{{ */
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1) {
        return;
    }
    int lower = precision - precision / 2;

    // recurse.
    root_lift(p, rx, irx, lower);

    if (WCT > w0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n",
                WCT,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    mpz_srcptr pl = power_lookup_const(p->powers, lower);
    mpz_ptr qk = my_barrett_init(pk);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    mpz_t fprime;
    mpz_init(fprime);

    // we know r to half-precision, 1/f'(r) to quarter-precision. Compute
    // f(r) (to full precision), and obtain 1/f'(r) to half-precision. We
    // compute f'(r) to half-precision as a side-effect of the
    // computation of f(r).
    mp_2poly_eval_mod(tb, fprime, f_hat, f_hat_diff, degree, degree-1, rx, pk, qk);
    /* use irx. only one iteration of newton.  */
    mpz_srcptr ql = NULL;       /* FIXME if we use barrett reduction */
    WRAP_barrett_mod(fprime, fprime, pl, ql);
    WRAP_mpz_mul(ta, irx, fprime);
    WRAP_barrett_mod(ta, ta, pl, ql);
    mpz_sub_ui(ta,ta,1);
    WRAP_mpz_submul(irx, irx, ta);
    WRAP_barrett_mod(irx, irx, pl, ql);

    mpz_clear(fprime);

    WRAP_mpz_mul(tb, irx, tb);
    mpz_sub(rx, rx, tb);
    WRAP_barrett_mod(rx, rx, pk, qk);

    mpz_clear(ta);
    mpz_clear(tb);

    my_barrett_clear(qk);
}/* }}} */

#define EXAMPLE_FROM_RSA768

int main(int argc, char * argv[])
{
    program_starttime = wct_seconds();

#ifdef  EXAMPLE_FROM_C159
    unsigned long p = 9223372036854799379UL;
    unsigned long r = 8197937682795196680UL;
    unsigned int precision = 1908226;
    const char * coeff_strings[] = {
        "33702859150680615562179587684939483",
        "2824683557949990996869965408188",
        "48702062013735864254260849",
        "-969379411138214076748",
        "-3170645766299012",
        "9108387600",
    };
#endif

#ifdef  EXAMPLE_FROM_RSA768
    unsigned long p = 9223372036854824267UL;
    unsigned long r = 3169229453447134476UL;
    unsigned int precision = 29069104;
    const char * coeff_strings[] = {
        "-277565266791543881995216199713801103343120",
        "-18185779352088594356726018862434803054",
        "6525437261935989397109667371894785",
        "-46477854471727854271772677450",
        "-5006815697800138351796828",
        "1276509360768321888",
        "265482057982680",
    };
#endif

    degree = sizeof(coeff_strings) / sizeof(coeff_strings[0]) - 1;

    if (argc >= 3 && strcmp(argv[1], "-prec") == 0) {
        precision = atoi(argv[2]);
    }

    fprintf(stderr, "GMP header: %d.%d.%d, using library %s\n",
            __GNU_MP_VERSION,
            __GNU_MP_VERSION_MINOR,
            __GNU_MP_VERSION_PATCHLEVEL,
            gmp_version);
    fprintf(stderr, "GMP compiled with %s, flags %s\n",
            __GMP_CC,
            __GMP_CFLAGS);
    /* These two are provided by cado_config.h */
    fprintf(stderr, "GMP includes from %s\n", GMP_INCDIR);
    fprintf(stderr, "GMP library from %s\n", GMP_LIBDIR);

#ifdef  WITH_BARRETT
    fprintf(stderr, "This run uses barrett reduction\n");
#else
    fprintf(stderr, "This run does not use barrett reduction\n");
#endif

    struct prime_data prime[1];
    prime->p = p;
    prime->powers = power_lookup_table_init(prime->p);

    f_hat = malloc((degree+1) * sizeof(mpz_t));
    f_hat_diff = malloc((degree) * sizeof(mpz_t));

    for(int i = 0 ; i <= degree ; i++) {
        mpz_init_set_str(f_hat[i], coeff_strings[i], 0);
    }
    {
        mpz_t w;
        mpz_init(w);
        mpz_set_ui(w,1);
        for(int i = degree-2 ; i >= 0 ; i--) {
            mpz_mul(w,w,f_hat[degree]);
            mpz_mul(f_hat[i],f_hat[i],w);
        }
        mpz_set_ui(f_hat[degree],1);
        mpz_clear(w);
    }

    for(int i = 0 ; i <= degree-1 ; i++) {
        mpz_init(f_hat_diff[i]);
        mpz_mul_ui(f_hat_diff[i], f_hat[i+1], i+1);
    }

    {
        char sbuf[32];
        fprintf(stderr, "# [%2.2lf] Lifting to precision l=%d (p^l is approx %s)\n", WCT, precision, size_disp(precision * log(prime->p)/M_LN2 / 8, sbuf));
    }

    fprintf(stderr, "# [%2.2lf] Computing powers of p\n", WCT);
    power_lookup(prime->powers, precision);
    fprintf(stderr, "# [%2.2lf] Computing powers of p: done.\n", WCT);

    program_starttime = wct_seconds();

    mpz_t rx, irx;
    mpz_init(rx);
    mpz_init(irx);

    mpz_set_ui(rx, r);

    mpz_srcptr p1 = power_lookup_const(prime->powers, 1);

    mp_poly_eval_mod(irx, f_hat_diff, degree-1, rx, p1, NULL);
    mpz_invert(irx, irx, p1);

    fprintf(stderr, "# [%2.2lf] Lifting (%lu,x-%lu)\n", WCT, p, r);

    {
        double t0, t1;
        double w0, w1;
        double rate;

        t0 = seconds();
        w0 = WCT;

        root_lift(prime, rx, irx, precision);

        t1 = seconds();
        w1 = WCT;
        rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
        fprintf(stderr, "# [%2.2lf] lift completed. %.2lf, wct %.2lf, %.1f%%\n",
                WCT, t1-t0, w1-w0, rate);
        fprintf(stderr, "# [%2.2lf] limb0 of lifted root is %lu\n",
                WCT, mpz_getlimbn(rx, 0));
    }

    mpz_clear(irx);
    mpz_clear(rx);

    for(int i = 0 ; i <= degree ; i++) {
        mpz_clear(f_hat[i]);
    }
    for(int i = 0 ; i <= degree-1 ; i++) {
        mpz_clear(f_hat_diff[i]);
    }

    free(f_hat);
    free(f_hat_diff);
}

