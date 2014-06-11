#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <math.h> /* for sqrt and floor and log and ceil */
#include "portability.h"
#include "utils.h"


/*
 * Compute g(x) = f(a*x+b), with deg f = d, and a and b are longs.
 * Must have a != 0 and d >= 1
 */
void mp_poly_linear_comp(mpz_t *g, mpz_t *f, int d, long a, long b) {
    ASSERT (a != 0  &&  d >= 1);
    // lazy: use the mpz_poly_t interface of utils/mpz_poly.h
    mpz_poly_t aXpb, aXpbi, G, Aux;
    mpz_poly_init(aXpb, 1);  // alloc sets to zero
    mpz_poly_init(aXpbi, d);
    mpz_poly_init(G, d);
    mpz_poly_init(Aux, d);
    {
        mpz_t aux;
        mpz_init(aux);
        mpz_set_si(aux, a);
        mpz_poly_setcoeff(aXpb, 1, aux);
        mpz_set_si(aux, b);
        mpz_poly_setcoeff(aXpb, 0, aux);
        mpz_clear(aux);
    }
    mpz_poly_copy(aXpbi, aXpb);
    mpz_poly_setcoeff(G, 0, f[0]);
    for (int i = 1; i <= d; ++i) {
        mpz_poly_mul_mpz(Aux, aXpbi, f[i]);
        mpz_poly_add(G, G, Aux);
        if (i < d)
            mpz_poly_mul(aXpbi, aXpbi, aXpb);
    }
    for (int i = 0; i <= d; ++i)
        mpz_poly_getcoeff(g[i], i, G);
    mpz_poly_clear(aXpb);
    mpz_poly_clear(aXpbi);
    mpz_poly_clear(G);
    mpz_poly_clear(Aux);
}

int mpz_p_val(mpz_t z, unsigned long p) {
    int v = 0;
    if (mpz_cmp_ui(z, 0) == 0)
        return INT_MAX;
    mpz_t zz;
    mpz_init(zz);
    mpz_abs(zz, z);
    unsigned long r;
    do {
        r = mpz_tdiv_q_ui(zz, zz, p);
        if (r == 0)
            v++;
    } while (r == 0 && mpz_cmp_ui(zz, 1)!=0);
    mpz_clear(zz);
    return v;
}

int mp_poly_p_val_of_content(mpz_t *f, int d, unsigned long p) {
    int val = INT_MAX;
    for (int i = 0; i <= d; ++i) {
        int v = mpz_p_val(f[i], p);
        if (v == 0)
            return 0;
        val = MIN(val, v);
    }
    return val;
}


void
mp_poly_eval (mpz_t r, mpz_t *poly, int deg, mpz_t a)
{
  int i;

  mpz_set (r, poly[deg]);
  for (i = deg - 1; i >= 0; i--)
    {
      mpz_mul (r, r, a);
      mpz_add (r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
void
mp_poly_eval_diff (mpz_t r, mpz_t *poly, int deg, mpz_t a)
{
  int i;

  mpz_mul_ui (r, poly[deg], (unsigned long)deg);
  for (i = deg - 1; i >= 1; i--)
    {
      mpz_mul (r, r, a);
      mpz_addmul_ui (r, poly[i], (unsigned long)i);
    }
}

// Same function as above, with a slightly different interface.
unsigned long
lift_root_unramified(mpz_t *f, int d, unsigned long r,
        unsigned long p,int kmax) {
    mpz_t aux, aux2, mp_p, mp_r;
    int k = 1;
    mpz_init(aux);
    mpz_init(aux2);
    mpz_init_set_ui(mp_r, r);
    mpz_init_set_ui(mp_p, p);
    while (k < kmax) {
        if (2*k <= kmax)
            mpz_mul(mp_p, mp_p, mp_p); // p^2k
        else {
            for (int i = k+1; i <= kmax; ++i)
                mpz_mul_ui(mp_p, mp_p, p);
        }
        mp_poly_eval(aux, f, d, mp_r);
        mp_poly_eval_diff(aux2, f, d, mp_r);
        if (!mpz_invert(aux2, aux2, mp_p)) {
            fprintf(stderr, "Error in lift_root_unramified: multiple root mod %lu\n", p);
            exit(EXIT_FAILURE);
        }
        mpz_mul(aux, aux, aux2);
        mpz_sub(aux, mp_r, aux);
        mpz_mod(mp_r, aux, mp_p);
        k *= 2;
    }
    r = mpz_get_ui(mp_r);
    mpz_clear(aux);
    mpz_clear(aux2);
    mpz_clear(mp_p);
    mpz_clear(mp_r);
    return r;
}

typedef struct {
    unsigned long q;
    unsigned long r;
    int n1;
    int n0;
} entry;

typedef struct {
    entry *list;
    int len;
    int alloc;
} entry_list;

void entry_list_init(entry_list *L) {
    L->list = (entry*) malloc(10*sizeof(entry));
    L->alloc = 10;
    L->len = 0;
}

void entry_list_clear(entry_list *L) {
    free(L->list);
}

void push_entry(entry_list *L, entry E) {
    if (L->len == L->alloc) {
        L->alloc += 10;
        L->list = (entry *)realloc(L->list, (L->alloc)*sizeof(entry));
    }
    L->list[L->len++] = E;
}

int cmp_entry(const void *A, const void *B) {
    entry a, b;
    a = ((entry *)A)[0];
    b = ((entry *)B)[0];
    if (a.q < b.q)
        return -1;
    if (a.q > b.q)
        return 1;
    if (a.n1 < b.n1)
        return -1;
    if (a.n1 > b.n1)
        return 1;
    if (a.n0 < b.n0)
        return -1;
    if (a.n0 > b.n0)
        return 1;
    if (a.r < b.r)
        return -1;
    if (a.r > b.r)
        return 1;
    return 0;
}

/* see makefb.sage for the meaning of this function and its parameters */
void all_roots_affine(entry_list *L, mpz_t *f, int d, unsigned long p,
        int kmax, int k0, int m, unsigned long phi1, unsigned long phi0) {
    int nroots;
    unsigned long *roots;
    mpz_t aux;

    mpz_poly_t F;
    F->coeff = f;
    F->deg = d;

    if (k0 >= kmax) {
        return;
    }
    roots = (unsigned long*) malloc(d * sizeof(unsigned long));
    mpz_init(aux);

    nroots = mpz_poly_roots_ulong (roots, F, p);
    for (int i = 0; i < nroots; ++i) {
        unsigned long r = roots[i];
        {
            mpz_t mp_r;
            mpz_init_set_ui(mp_r, r);
            mp_poly_eval_diff(aux, f, d, mp_r);
            mpz_clear(mp_r);
        }
        unsigned long dfr = mpz_mod_ui(aux, aux, p);
        if (dfr != 0) {
            unsigned long rr = lift_root_unramified(f, d, r, p, kmax-k0);
            unsigned long phir = phi1 * rr + phi0;
            unsigned long pml = 1;
            for (int j = 0; j < m; ++j)
                pml *= p;
            for (int l = 1; l <= kmax-k0; ++l) {
                pml *= p;
                unsigned long phirr = phir % pml;
                entry E;
                E.q = pml;
                E.r = phirr;
                E.n1 = k0+l;
                E.n0 = k0+l-1;
                push_entry(L, E);
            }
        } else {
            mpz_t *ff;
            ff = (mpz_t *) malloc((d+1)*sizeof(mpz_t));
            for (int j = 0; j <= d; ++j)
                mpz_init(ff[j]);
            mp_poly_linear_comp(ff, f, d, p, r);
            int val = mp_poly_p_val_of_content(ff, d, p);
            unsigned long pmp1 = 1;
            for (int j = 0; j < m+1; ++j)
                pmp1 *= p;
            unsigned long phir = (phi1 * r + phi0) % pmp1;
            entry E;
            E.q = pmp1;
            E.r = phir;
            E.n1 = k0+val;
            E.n0 = k0;
            push_entry(L, E);
            unsigned long nphi1 = phi1*p;
            unsigned long nphi0 = phi0 + phi1*r;
            {
                mpz_t pv;
                mpz_init(pv);
                mpz_set_ui(pv, 1);
                for (int j = 0; j < val; ++j)
                    mpz_mul_ui(pv, pv, p);
                for (int j = 0; j <= d; ++j)
                    mpz_tdiv_q(ff[j], ff[j], pv);
                mpz_clear(pv);
            }
            all_roots_affine(L, ff, d, p, kmax, k0+val, m+1, nphi1, nphi0);
            for (int j = 0; j <= d; ++j)
                mpz_clear(ff[j]);
            free(ff);
        }
    }
    free(roots);
    mpz_clear(aux);
}

/* Compute roots mod powers of p, up to maxbits.
 * TODO:
 * Maybe the number of bits of the power is not really the right measure.
 * We could take into account not only the cost of the updates (which is
 * indeed in k*log(p)), but also the amount of information we gain for
 * each update, which is in log(p). So what would be the right measure
 * ???
 */

/*
 * See makefb.sage for details on this function
 */

entry_list all_roots(mpz_t *f, int d, unsigned long p, int maxbits) {
    int kmax;
    entry_list L;
    entry_list_init(&L);
    kmax = (int)floor(maxbits*log(2)/log(p));
    if (kmax == 0)
        kmax = 1;
    { // handle projective roots first.
        mpz_t *fh;
        mpz_t pk;
        fh = (mpz_t *)malloc ((d+1)*sizeof(mpz_t));
        mpz_init(pk);
        mpz_set_ui(pk, 1);
        for (int i = 0; i <= d; ++i) {
            mpz_init(fh[i]);
            mpz_mul(fh[i], f[d-i], pk);
            if (i < d)
                mpz_mul_ui(pk, pk, p);
        }
        int v = mp_poly_p_val_of_content(fh, d, p);
        if (v > 0) { // We have projective roots only in that case
            {
                entry E;
                E.q = p;
                E.r = p;
                E.n1 = v;
                E.n0 = 0;
                push_entry(&L, E);
            }
            mpz_set_ui(pk, p);
            for (int i = 1; i < v; ++i)
                mpz_mul_ui(pk, pk, p);
            for (int i = 0; i <= d; ++i)
                mpz_tdiv_q(fh[i], fh[i], pk);

            all_roots_affine(&L, fh, d, p, kmax-1, 0, 0, 1, 0);
            // convert back the roots
            for (int i = 1; i < L.len; ++i) {
                entry E = L.list[i];
                E.q *= p;
                E.n1 += v;
                E.n0 += v;
                E.r = E.r*p + E.q;
                L.list[i] = E;
            }
        }
        for (int i = 0; i <= d; ++i)
            mpz_clear(fh[i]);
        mpz_clear(pk);
        free(fh);
    }
    // affine roots are easier.
    all_roots_affine(&L, f, d, p, kmax, 0, 0, 1, 0);

    qsort((void *)(&L.list[0]), L.len, sizeof(entry), cmp_entry);
    return L;
}


void makefb_with_powers(FILE* outfile, mpz_poly_t F, unsigned long alim,
        int maxbits) {
    mpz_t *f = F->coeff;
    int d = F->deg;

    fprintf(outfile, "# Roots for polynomial ");
    mpz_poly_fprintf(outfile, F);
    fprintf(outfile, "# DEGREE: %d\n", d);
    fprintf(outfile, "# alim = %lu\n", alim);
    fprintf(outfile, "# maxbits = %d\n", maxbits);

    entry_list L;
    unsigned long p;
    for (p = 2; p <= alim; p = getprime (p)) {
        L = all_roots(f, d, p, maxbits);
        // print in a compactified way
        int oldn0=-1, oldn1=-1;
        unsigned long oldq = 0;
        for (int i = 0; i < L.len; ++i) {
            unsigned long q = L.list[i].q;
            int n1 = L.list[i].n1;
            int n0 = L.list[i].n0;
            unsigned long r =  L.list[i].r;
            if (q == oldq && n1 == oldn1 && n0 == oldn0)
                fprintf(outfile, ",%lu", r);
            else {
                if (i > 0)
                    fprintf(outfile, "\n");
                oldq = q; oldn1 = n1; oldn0 = n0;
                if (n1 == 1 && n0 == 0)
                    fprintf(outfile, "%lu: %lu", q, r);
                else
                    fprintf(outfile, "%lu:%d,%d: %lu", q, n1, n0, r);
            }
        }
        if (L.len > 0)
            fprintf(outfile, "\n");
        entry_list_clear(&L);
    }
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "polynomial file");
    param_list_decl_usage(pl, "alim", "factor base bound");
    param_list_decl_usage(pl, "maxbits", "(optional) maximal number of "
            "bits of powers");
    param_list_decl_usage(pl, "out", "(optional) name of the output file");
    char str[1024];
    sprintf(str, "(optional) create factor base for given side. "
            "Side must be %d or %d (default is %d, i.e. algebraic).",
            RATIONAL_SIDE, ALGEBRAIC_SIDE, ALGEBRAIC_SIDE);
    param_list_decl_usage(pl, "side", str);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly cpoly;
  FILE * f, *outputfile;
  const char *outfilename = NULL;
  int maxbits = 1;  // disable powers by default
  int side = ALGEBRAIC_SIDE;
  unsigned long alim = 0;
  char *argv0 = argv[0];

  param_list_init(pl);
  declare_usage(pl);
  cado_poly_init(cpoly);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_ulong(pl, "alim", &alim);
  if (alim == 0) {
      fprintf(stderr, "Error: parameter -alim is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "side", &side);
  if (side != RATIONAL_SIDE && side != ALGEBRAIC_SIDE) {
      fprintf(stderr, "Error: side must be %d (for rational) or %d "
              "(for algebraic)\n", RATIONAL_SIDE, ALGEBRAIC_SIDE);
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "maxbits", &maxbits);

  outfilename = param_list_lookup_string(pl, "out");
  if (outfilename != NULL) {
    outputfile = fopen_maybe_compressed(outfilename, "w");
  } else {
    outputfile = stdout;
  }


  if (!cado_poly_read(cpoly, filename))
    {
      fprintf (stderr, "Error reading polynomial file %s\n", filename);
      exit (EXIT_FAILURE);
    }
  if (side == ALGEBRAIC_SIDE)
      makefb_with_powers (outputfile, cpoly->alg, alim, maxbits);
  else
      makefb_with_powers (outputfile, cpoly->rat, alim, maxbits);

  cado_poly_clear (cpoly);
  if (outfilename != NULL) {
    fclose_maybe_compressed(outputfile, outfilename);
  }

  param_list_clear(pl);

  return 0;
}
