#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "macros.h"

#include "types.h"
#include "ffspol.h"
#include "fqpol.h"
#include "params.h"
#include "polyfactor.h"

#define MAX_FFS_DEG 10

#define MAX_NPOL 16

// TODO: we use sq_t as the basic type, because this is the one used
// within the fq_t implementation. The latter choice is probably not
// right, but we leave it for the moment.

typedef struct {
    sq_t q;
    sq_t r;
    int n1;
    int n0;
} entry;

typedef struct {
    entry *list;
    int len;
    int alloc;
} entry_list_struct;

typedef entry_list_struct entry_list_t[1];
typedef entry_list_struct * entry_list_ptr;
typedef const entry_list_struct * entry_list_srcptr;

void entry_list_init(entry_list_ptr L) {
    L->list = (entry*) malloc(10*sizeof(entry));
    L->alloc = 10;
    L->len = 0;
}

void entry_list_clear(entry_list_ptr L) {
    free(L->list);
}

void push_entry(entry_list_t L, sq_srcptr q, sq_srcptr r, int n1, int n0) {
    if (L->len == L->alloc) {
        L->alloc += 10;
        L->list = (entry *)realloc(L->list, (L->alloc)*sizeof(entry));
    }
    sq_set(L->list[L->len].q, q);
    sq_set(L->list[L->len].r, r);
    L->list[L->len].n1 = n1;
    L->list[L->len].n0 = n0;
    L->len++;
}

int cmp_entry(const void *A, const void *B) {
    entry a, b;
    a = ((entry *)A)[0];
    b = ((entry *)B)[0];
    if (sq_cmp(a.q, b.q) < 0)
        return -1;
    if (sq_cmp(a.q, b.q) > 0)
        return 1;
    if (a.n1 < b.n1)
        return -1;
    if (a.n1 > b.n1)
        return 1;
    if (a.n0 < b.n0)
        return -1;
    if (a.n0 > b.n0)
        return 1;
    if (sq_cmp(a.r, b.r) < 0)
        return -1;
    if (sq_cmp(a.r, b.r) > 0)
        return 1;
    return 0;
}

void entry_list_out(FILE *file, entry_list_srcptr E) {
    if (E->len == 0)
        return;
    int oldn0=-1, oldn1=-1;
    sq_t oldq;
    sq_set_zero(oldq);
    for (int i = 0; i < E->len; ++i) {
        sq_t q, r;
        sq_set(q, E->list[i].q);
        sq_set(r, E->list[i].r);
        int n1 = E->list[i].n1;
        int n0 = E->list[i].n0;
        if (sq_eq(q, oldq) && n1 == oldn1 && n0 == oldn0) {
            fprintf(file, ",");
            sq_out(file, r);
        } else {
            if (i > 0)
                fprintf(file, "\n");
            sq_set(oldq, q);
            oldn1 = n1; oldn0 = n0;
            if (n1 == 1 && n0 == 0) {
                sq_out(file, q);
                fprintf(file, ": ");
                sq_out(file, r);
            } else {
                sq_out(file, q);
                fprintf(file, ":%d,%d: ", n1, n0);
                sq_out(file, r);
            }
        }
    }
    fprintf(file, "\n");
}

// TODO: share this with sieve.c
static int sq_is_irreducible(sq_srcptr p) {
    fppol_t P;
    fppol_init(P);
    fppol_set_sq(P, p);
    int ret = fppol_is_irreducible(P);
    fppol_clear(P);
    return ret;
}


void lift_root_unramified(sq_ptr rr, ffspol_srcptr F, sq_srcptr r,
        sq_srcptr p, int kmax)
{
    fppol_t aux, aux2, P, R;
    int k = 1;
    fppol_init(aux);
    fppol_init(aux2);
    fppol_init(P);
    fppol_init(R);
    fppol_set_sq(R, r);
    fppol_set_sq(P, p);

    while (k < kmax) {
        if (2*k <= kmax)
            fppol_mul(P, P, P); // p^2k
        else {
            for (int i = k+1; i <= kmax; ++i)
                fppol_mul_sq(P, P, p);
        }
        ffspol_eval(aux, F, R);
        ffspol_eval_diff(aux2, F, R);
        fppol_rem(aux2, aux2, P);
        if (!fppol_invmod(aux2, aux2, P)) {
            fprintf(stderr, "Error in lift_root_unramified: multiple root.\n");
            exit(EXIT_FAILURE);
        }
        fppol_mul(aux, aux, aux2);
        fppol_sub(aux, R, aux);
        fppol_rem(R, aux, P);
        k *= 2;
    }
    sq_set_mp(rr, R);
    fppol_clear(aux);
    fppol_clear(aux2);
    fppol_clear(P);
    fppol_clear(R);
}

// FF(x) := F(phi1 * x + phi0)
void ffspol_linear_comp_sq(ffspol_ptr FF, ffspol_srcptr F,
        sq_t phi1, sq_t phi0)
{
    int d = F->deg;
    ASSERT_ALWAYS(!sq_is_zero(phi1));
    ffspol_t phi, phik, aux;
    ffspol_init2(phi, 2);
    ffspol_init2(phik, d + 1);
    ffspol_init2(aux, d + 1);

    phi->deg = 1;
    fppol_set_sq(phi->coeffs[0], phi0);
    fppol_set_sq(phi->coeffs[1], phi1);

    phik->deg = 0;
    fppol_set_ti(phik->coeffs[0], 0);

    ASSERT_ALWAYS((int)FF->alloc >= (d+1));

    FF->deg = 0;
    fppol_set(FF->coeffs[0], F->coeffs[0]);

    for (int k = 1; k <= d; ++k) {
        ffspol_mul(phik, phik, phi);
        ffspol_smul(aux, phik, F->coeffs[k]);
        ffspol_add(FF, FF, aux);
    }

    ffspol_clear(phi);
    ffspol_clear(phik);
    ffspol_clear(aux);
}

int fppol_pval(fppol_t z, sq_srcptr p) {
    int v = 0;
    if (fppol_is_zero(z))
        return INT_MAX;
    fppol_t zz, pp, r;
    fppol_init(zz);
    fppol_init(pp);
    fppol_init(r);
    fppol_set(zz, z);
    fppol_set_sq(pp, p);
    do {
        fppol_divrem(zz, r, zz, pp);
        if (fppol_is_zero(r))
            v++;
    } while (fppol_is_zero(r) && fppol_deg(zz) > 0);
    fppol_clear(zz);
    fppol_clear(pp);
    fppol_clear(r);
    return v;
}

int ffspol_pvaluation_sq(ffspol_srcptr F, sq_srcptr p)
{
    int val = INT_MAX;
    for (int i = 0; i <= F->deg; ++i) {
        int v = fppol_pval(F->coeffs[i], p);
        if (v == 0)
            return 0;
        val = MIN(val, v);
    }
    return val;
}


// Return 0 if we detect that there were projective roots.
// Otherwise return 1 (success: we finished the job)
// -1 means don't know; should occurr only in recursive calls.
int all_roots_affine(entry_list_ptr E, ffspol_srcptr F,
        sq_t p, int kmax, int k0, int m, sq_t phi1, sq_t phi0)
{
    if (k0 >= kmax)
        return -1;

    int has_proj;

    fq_t roots[MAX_FFS_DEG];
    fq_info_t Fq;
    fq_info_init(Fq, p);
    fqpol_t Fbar;
    fqpol_init(Fbar);
    fqpol_set_ffspol(Fbar, F, Fq);
    has_proj = (Fbar->deg != F->deg);

    int nr = fqpol_roots(roots, Fbar, Fq);
    for (int i = 0; i < nr; ++i) {
        int nmult = fqpol_root_multiplicity(Fbar, roots[i], Fq);
        if (nmult == 1) {
            // the typical case if kmax=1, k0=0.
            // In that case, no lift is needed.
            if (LIKELY(kmax == 1 && k0 == 0))
                push_entry(E, p, roots[i], 1, 0);
            else {
                sq_t rr, phir, pml;
                lift_root_unramified(rr, F, roots[i], p, kmax-k0);
                sq_mul(phir, phi1, rr);
                sq_add(phir, phir, phi0);
                sq_set_ti(pml, 0);
                for (int j = 0; j < m; ++j)
                    sq_mul(pml, pml, p);
                for (int l = 1; l <= kmax-k0; ++l) {
                    sq_mul(pml, pml, p);
                    sq_t phirr;
                    sq_rem(phirr, phir, pml);
                    push_entry(E, pml, phirr, k0+l, k0+l-1);
                }
            }
        } else {
            sq_t r;
            sq_set(r, roots[i]);
            ffspol_t FF;
            ffspol_init2(FF, F->alloc);
            ffspol_linear_comp_sq(FF, F, p, r);
            int val = ffspol_pvaluation_sq(FF, p);

            sq_t pmp1;
            sq_set_ti(pmp1, 0);
            for (int j = 0; j < m+1; ++j)
                sq_mul(pmp1, pmp1, p);
            sq_t phir;
            sq_mul(phir, phi1, r);
            sq_add(phir, phir, phi0);
            sq_rem(phir, phir, pmp1);
            push_entry(E, pmp1, phir, k0+val, k0);
            sq_t nphi0, nphi1;
            sq_mul(nphi1, phi1, p);
            sq_mul(nphi0, phi1, r);
            sq_add(nphi0, nphi0, phi0);

            {
                fppol_t pv;
                fppol_init(pv);
                fppol_set_ti(pv, 0);
                for (int j = 0; j < val; ++j)
                    fppol_mul_sq(pv, pv, p);
                for (int j = 0; j <= FF->deg; ++j)
                    fppol_div(FF->coeffs[j], FF->coeffs[j], pv);
                fppol_clear(pv);
            }
            all_roots_affine(E, FF, p, kmax, k0+val, m+1, nphi1, nphi0);
            ffspol_clear(FF);
        }
    }
    fqpol_clear(Fbar);
    return !has_proj;
}

void all_roots(entry_list_ptr E, ffspol_t F, sq_t p, int powerlim)
{
    int kmax = powerlim / sq_deg(p);
    if (kmax == 0)
        kmax = 1;

    sq_t phi0, phi1;
    sq_set_zero(phi0);
    sq_set_ti(phi1, 0);
    int no_proj = all_roots_affine(E, F, p, kmax, 0, 0, phi1, phi0);

    if (!no_proj) {
        int d = F->deg;
        ffspol_t FF;
        ffspol_init2(FF, d+1);
        fppol_t pk;
        fppol_init(pk);
        fppol_set_ti(pk, 0);
        for (int i = 0; i <= d; ++i) {
            fppol_mul(FF->coeffs[i], F->coeffs[d-i], pk);
            if (i < d)
                fppol_mul_sq(pk, pk, p);
        }
        FF->deg = d;
        int val = ffspol_pvaluation_sq(FF, p);
        ASSERT_ALWAYS(val > 0);
        push_entry(E, p, p, val, 0);

        fppol_set_sq(pk, p);
        for (int i = 1; i < val; ++i)
            fppol_mul_sq(pk, pk, p);
        for (int i = 0; i <= d; ++i)
            fppol_div(FF->coeffs[i], FF->coeffs[i], pk);

        entry_list_t EE;
        entry_list_init(EE);
        all_roots_affine(EE, FF, p, kmax-1, 0, 0, phi1, phi0);
        // convert back the roots.
        for (int i = 0; i < EE->len; ++i) {
            sq_t pp, rr;
            sq_mul(pp, EE->list[i].q, p);
            sq_mul(rr, EE->list[i].r, p);
            sq_add(rr, rr, pp);
            push_entry(E, pp, rr, EE->list[i].n1+val, EE->list[i].n0+val);
        }
        entry_list_clear(EE);

        fppol_clear(pk);
        ffspol_clear(FF);
    }

    qsort((void *)(E->list), E->len, sizeof(entry), cmp_entry);
}

double alpha_p(entry_list_srcptr E, sq_srcptr p, int powerlim) 
{
    static double * pows = NULL;
    if (pows == NULL) {
        pows = (double *) malloc (20*sizeof(double));
        pows[0] = 1.0;
        pows[1] = 1.0/((double)FP_SIZE);
        for (int i = 2; i < 20; ++i) 
            pows[i] = pows[i-1]*pows[1];
    }

    int degp = sq_deg(p);
    double alpha = 0.0;

    // initialize with value for rational polynomial (only up to powerlim)
    alpha += pows[degp];
    for (int i = 2; i*degp < powerlim; ++i)
        alpha += pows[i*degp];

    // subtract contribution of each entry
    if (E->len != 0) {
        double al = 0.0;
        for (int i = 0; i < E->len; ++i)
            al -= (E->list[i].n1 - E->list[i].n0) * pows[sq_deg(E->list[i].q)];
        alpha += al / (pows[degp]+1); // coprimality of (a,b)
    }

    // normalize alpha
    alpha *= degp;  // multiply the valuation by the degree.

    return alpha;
}

void makefb(ffspol_t F, int fbb, const char *filename, int powerlim)
{
    FILE *file;
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error when opening factor base file for writing\n");
        exit(EXIT_FAILURE);
    }

    double alpha = 0.0;
    sq_t p;
    sq_set_ti(p, 0);
    sq_monic_set_next(p, p, 64);
    while (sq_deg(p) <= fbb) {
        entry_list_t E;
        entry_list_init(E);
        E->len = 0;
        all_roots(E, F, p, powerlim);
        alpha += alpha_p(E, p, powerlim);
        entry_list_out(file, E);
        entry_list_clear(E);
        do {
            sq_monic_set_next(p, p, 64);
        } while (!sq_is_irreducible(p));
    }
    fclose(file);
    printf("alpha = %f\n", alpha);
}

void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0             function field polynomial(s) on side 0\n");
    fprintf(stderr, "  pol1             function field polynomial(s) on side 1\n");
    fprintf(stderr, "  fb0              factor base file(s) on side 0\n");
    fprintf(stderr, "  fb1              factor base file(s) on side 1\n");
    fprintf(stderr, "  fbb0             factor base bound(s) on side 0\n");
    fprintf(stderr, "  fbb1             factor base bound(s) on side 1\n");
    fprintf(stderr, "  side             compute fb only for given side\n");
    fprintf(stderr, "  powerlim [1]     put powers up to given degree\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");

    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    int gf = 0;
    int powerlim = 1;
    char *argv0 = argv[0];
    param_list pl;
    param_list_init(pl);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }
//    param_list_print_command_line(stdout, pl);

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
    }

    param_list_parse_int(pl, "powerlim", &powerlim);

    int side = 0, npol[2];

    ffspol_t ffspol[2][MAX_NPOL];

    // read function field polynomials
    for (int i = 0; i < 2; ++i) {
        const char * polstr[MAX_NPOL];
        npol[i] = param_list_lookup_string_list(pl, i ? "pol1" : "pol0",
                                                polstr, MAX_NPOL, ";");
        if (npol[i]) {
          side |= 1<<i;
          for (int j = 0; j < npol[i]; ++j) {
            ffspol_init(ffspol[i][j]);
            ffspol_set_str(ffspol[i][j], polstr[j]);
          }
        }
        param_list_restore_string_list(pl, i ? "pol1" : "pol0",
                                       polstr, npol[i], ";");
    }
    if (side == 0) {
        fprintf(stderr, "Please provide at least one polynomial\n");
        usage(argv0, NULL);
    }
    if (npol[0] > 1 && npol[1] > 1) {
      fprintf(stderr, "Both sides cannot have multiple polynomials\n");
      exit(EXIT_FAILURE);
    }

    // The side parameter can disable one side.
    int sideparam = -1;
    param_list_parse_int(pl, "side", &sideparam);
    if (sideparam == 0) {
        if (!(side & 1))
            usage(argv0, "pol0");
        side = 1;
    }
    if (sideparam == 1) {
        if (!(side & 2))
            usage(argv0, "pol1");
        side = 2;
    }

    // Process each side.
    for (int i = 0; i < 2; ++i) {
      if (side & (1<<i)) {
        int fbb[MAX_NPOL];
        int n = param_list_parse_int_list(pl, i ? "fbb1" : "fbb0",
                                          fbb, npol[i], ";");
        if (!n)
            usage(argv0, i ? "fbb1" : "fbb0");
        for (; n < npol[i]; ++n)
          fbb[n] = fbb[n-1];
        const char *filename[MAX_NPOL];
        n = param_list_lookup_string_list(pl, i ? "fb1" : "fb0",
                                          filename, npol[i], ";");
        if (!n)
            usage(argv0, i ? "fb1" : "fb0");
        if (n != npol[i]) {
          fprintf(stderr, "Not enough factor base files specified "
                  "for side %d\n", i);
          exit(EXIT_FAILURE);
        }
        for (int j = 0; j < npol[i]; ++j)
          makefb(ffspol[i][j], fbb[j], filename[j], powerlim);
        param_list_restore_string_list(pl, i ? "fb1" : "fb0",
                                       filename, npol[i], ";");
      }
    }

    param_list_clear(pl);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < npol[i]; ++j)
        ffspol_clear(ffspol[i][j]);
    return 0;
}
