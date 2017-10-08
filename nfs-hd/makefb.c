#include "cado.h"
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <inttypes.h>
#include "utils.h"
#include "macros.h"
#include "getprime.h"
#include "makefb.h"
#include "utils_int64.h"
#include "math.h"

/*
 * Mode:
 *   MAIN: export factor bases. Default file name is: <p>,<n>,<#numberfield>.
 *   TIME_MAKEFB: get the time to build factor bases.
 */

//mpz_poly_factor does not work in characteristic 2, so we do it the naive way.

/*
 * Set an ideal_1 at an index.
 *
 * fb: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 * h: the h of the ideal (r, h).
 * fbb: factor base bound for this side.
 */
static void add_ideal_1_part(factor_base_ptr fb, uint64_t * index, uint64_t r,
    mpz_poly_srcptr h, uint64_t fbb, unsigned int t)
{
  ASSERT(h->deg == 1);

  //Verify if the ideal can be added.
  if (r <= fbb) {
    factor_base_set_ideal_1_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
}

/*
 * Set an ideal_u at an index.
 *
 * fb: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 * h: the h of the ideal (r, h).
 * fbb: factor base bound for this side.
 * lpb: large prime bound.
 */
static void add_ideal_u_part(factor_base_ptr fb, uint64_t * index, uint64_t r,
    mpz_poly_srcptr h, uint64_t fbb, mpz_t lpb, unsigned int t)
{
  ASSERT(h->deg > 1);

  //Verify if the ideal can be added.
  mpz_t limit;
  mpz_init(limit);
  mpz_set_d(limit, pow((double)r, (double)h->deg));
  if (mpz_cmp(lpb, limit) >= 0 && r <= fbb) {
    factor_base_set_ideal_u_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
  mpz_clear(limit);
}

/*
 * Set an ideal_pr at an index.
 *
 * fb: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 * fbb: factor base bound for this side.
 */
static void add_ideal_pr_part(factor_base_ptr fb, uint64_t * index, uint64_t r,
    uint64_t fbb, unsigned int t)
{
  //Verify if the ideal can be added.
  if (r <= fbb) {
    factor_base_set_ideal_pr(fb, * index, r, t);
    * index = * index + 1;
  }
}

//WARNING: untested code.
void makefb(factor_base_t * fb, cado_poly_srcptr f, uint64_t * fbb,
    unsigned int t, unsigned int * lpb_bit)
{
  unsigned int V = f->nb_polys;

  ASSERT(V >= 2);
  ASSERT(t >= 2);

  mpz_t * lpb = (mpz_t *) malloc(V * sizeof(mpz_t));
  for (unsigned int i = 0; i < V; i++) {
    mpz_init(lpb[i]);
    mpz_ui_pow_ui(lpb[i], 2, (unsigned long) lpb_bit[i]);
  }

#ifndef NDEBUG
  for (unsigned int i = 0; i < V; i++) {
    ASSERT(f->pols[i]->deg >= 1);
    ASSERT(fbb[i] > 2);
    ASSERT(mpz_cmp_ui(lpb[i], fbb[i]) >= 0);
  }
#endif // NDEBUG

  //Factorise in Fq.
  uint64_t q = 2;
  //Count number of ideal_1, ideal_u and ideal_pr for each sides.
  uint64_t * index1 = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  uint64_t * indexu = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  uint64_t * indexpr = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  gmp_randstate_t state;
  //a = q in mpz. We need a to use mpz_poly_factor.
  mpz_t a;
  mpz_poly_factor_list l;

  mpz_t zero;
  mpz_init(zero);
  mpz_set_ui(zero, 0);

  //Contains the leading coefficient of f->pols[k].
  mpz_t lc;
  mpz_init(lc);

  for (unsigned int k = 0; k < V; k++) {
    index1[k] = 0;
    indexu[k] = 0;
    indexpr[k] = 0;
  }

  gmp_randinit_default(state);
  mpz_init(a);
  mpz_poly_factor_list_init(l);

  prime_info pi;
  prime_info_init(pi);

  //Find the maximum of fbb.
  uint64_t qmax = fbb[0];
  for (unsigned int k = 1; k < V; k++) {
    qmax = MAX(qmax, fbb[k]);
  }

  //Approximative number of enumerated primes.
#ifndef PRINT_INFO
  uint64_t nb_max = (uint64_t) ((double)qmax / (log((double)qmax)));
  uint64_t cpt = 0;
  unsigned int percent = 0;
#endif // PRINT_INFO

  //For all the prime q less than the max of fbb.
  for ( ; q <= qmax; q = getprime_mt(pi)) {
#ifndef PRINT_INFO
    cpt++;
    if ((double)cpt / (double)nb_max * 100.0 >= (double)percent && percent <
        101) {
      printf("# [%3u %%] Building factor bases.\n", percent);
      percent += 10;
    }
#endif // PRINT_INFO
    mpz_set_ui(a, q);
    for (unsigned int k = 0; k < V; k++) {
      //Projective root?
      mpz_set(lc, mpz_poly_lc_const(f->pols[k]));
      if (mpz_congruent_p(lc, zero, a) != 0) {
        add_ideal_pr_part(fb[k], indexpr + k, q, fbb[k], t);
      }
      if (q <= fbb[k]) {
        //Factorization of f->pols[k] mod q.
        mpz_poly_factor(l, f->pols[k], a, state);
        for (int i = 0; i < l->size ; i++) {
          if (l->factors[i]->f->deg == 1) {
            add_ideal_1_part(fb[k], index1 + k, q, l->factors[i]->f, fbb[k], t);
          } else if (l->factors[i]->f->deg < (int)t) {
            add_ideal_u_part(fb[k], indexu + k, q, l->factors[i]->f, fbb[k],
                lpb[k], t);
          }
        }
      }
    }
  }

  mpz_poly_factor_list_clear(l);

  //Realloc the factor base with just the number of each type of ideal.
  for (unsigned int k = 0; k < V; k++) {
    factor_base_realloc(fb[k], index1[k], indexu[k], indexpr[k]);
  }

  mpz_clear(lc);
  mpz_clear(zero);
  free(index1);
  free(indexu);
  free(indexpr);
  gmp_randclear(state);
  mpz_clear(a);

  for (unsigned int i = 0; i < V; i++) {
    mpz_clear(lpb[i]);
  }
  free(lpb);
  prime_info_clear (pi);
}

/* Parse different types. */

static int parse_ulong(unsigned long * x, char ** endptr, char * ptr)
{
  unsigned long xx;
  errno = 0;
  xx = strtoul(ptr, endptr, 10);
  if (errno) {
    // failure
    return 0;
  }
  *x = xx;
  return 1;
}

static int parse_uint64(uint64_t * x, char ** endptr, char * ptr)
{
  return parse_ulong((unsigned long *) x, endptr, ptr);
}

/*
 * Read z from ptr and advance pointer.
 * Assume z has been initialized.
 * Return 0 or 1 for failure / success
 */
static int parse_mpz(mpz_t z, char ** endptr, char * ptr)
{
  int r = gmp_sscanf(ptr, "%Zd", z);
  if (r != 1) {
    *endptr = ptr;
    return 0; // failure
  }
  *endptr = ptr;
  while (isdigit(*endptr[0]) || *endptr[0] == '-') {
    (*endptr)++;
  }
  return 1;
}

/*
 * Read z from ptr and advance pointer.
 * Assume z has been initialized.
 * Return 0 or 1 for failure / success
 */
static int parse_int(int * i, char ** endptr, char * ptr)
{
  int r = sscanf(ptr, "%d", i);
  if (r != 1) {
    *endptr = ptr;
    return 0; // failure
  }
  *endptr = ptr;
  while (isdigit(*endptr[0]) || *endptr[0] == '-') {
    (*endptr)++;
  }
  return 1;
}

/*
 * Read z0,z1,...,zk from ptr and advance pointer.
 * Assume all the zi have been initialized and that the array is large
 * enough to hold everyone.
 * Return the number of zi that are parsed (0 means error or empty list)
 */
static int parse_cs_mpzs(mpz_t *z, char ** endptr, char * ptr)
{
  char *myptr = ptr;
  int cpt = 0;
  for(;;) {
    int ret = parse_mpz(z[cpt], endptr, myptr);
    if (!ret) {
      return 0; // failure or empty list
    }
    // got an mpz
    cpt++;
    myptr = *endptr;
    if (myptr[0] != ',') {
      // finished!
      *endptr = myptr;
      return cpt;
    }
    // prepare for next mpz
    myptr++;
  }
}

static int parse_mpz_poly(mpz_poly_ptr f, char ** endptr, char * str, int n)
{
  int ret;

  mpz_t * coeffs = (mpz_t *) malloc(sizeof(mpz_t) * (n + 1));
  for (int i = 0; i <= n; i++) {
    mpz_init(coeffs[i]);
  }
  ret = parse_cs_mpzs(coeffs, endptr, str);
  mpz_poly_setcoeffs(f, coeffs, n);
  for (int i = 0; i <= n; i++) {
    mpz_clear(coeffs[i]);
  }
  free(coeffs);

  return ret;
}

static int parse_line_mpz_poly(mpz_poly_ptr f, char * str, int n)
{
  char * tmp;
  if (str[0] != 'f' || str[1] != ':')
    return 0;
  str += 2;

  return parse_mpz_poly(f, &tmp, str, n);
}

static int parse_ideal_1(ideal_1_ptr ideal, char *str, unsigned int t)
{
  char *tmp;
  int ret;

  if (str[0] != '1' || str[1] != ':')
    return 0;
  str += 2;
  uint64_t r;
  ret = parse_uint64(&r, &tmp, str);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  mpz_poly h;
  mpz_poly_init(h, 1);
  ret = parse_mpz_poly(h, &tmp, str, 1);
  ASSERT(ret == 2);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  mpz_t * Tr = (mpz_t *) malloc(sizeof(mpz_t) * (t - 1));
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_init(Tr[i]);
  }
  ret = parse_cs_mpzs(Tr, &tmp, str);
  ASSERT(ret == (int)(t - 1));

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;

  double log = atof(str);

  ideal_1_set_element(ideal, r, h, Tr, log, t);

  mpz_poly_clear(h);
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(Tr[i]);
  }
  free(Tr);
  return ret;
}

static int parse_ideal_u(ideal_u_ptr ideal, char * str, unsigned int t)
{
  char *tmp;
  int ret;

  int deg;
  ret = parse_int(&deg, &tmp, str);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;

  uint64_t r;
  ret = parse_uint64(&r, &tmp, str);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  mpz_poly h;
  mpz_poly_init(h, deg);
  ret = parse_mpz_poly(h, &tmp, str, deg);
  ASSERT(ret == deg + 1);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  int size_Tr = ((int)t - deg) * deg;
  mpz_t * Tr = (mpz_t *) malloc(sizeof(mpz_t) * size_Tr);
  for (int i = 0; i < size_Tr; i++) {
    mpz_init(Tr[i]);
  }
  ret = parse_cs_mpzs(Tr, &tmp, str);
  ASSERT(ret == size_Tr);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;

  double log = atof(str);

  ideal_u_set_element(ideal, r, h, Tr, log, t);

  mpz_poly_clear(h);
  for (int i = 0; i < size_Tr; i++) {
    mpz_clear(Tr[i]);
  }
  free(Tr);
  return ret;
}

static int parse_ideal_pr(ideal_pr_ptr ideal, char *str, unsigned int t)
{
  char *tmp;
  int ret;

  if (str[0] != '1' || str[1] != ':')
    return 0;
  str += 2;
  uint64_t r;
  ret = parse_uint64(&r, &tmp, str);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  mpz_poly h;
  mpz_poly_init(h, 1);
  ret = parse_mpz_poly(h, &tmp, str, 1);
  ASSERT(ret == 2);
  ASSERT(mpz_cmp_ui(h->coeff[0], 0) == 0);
  ASSERT(mpz_cmp_ui(h->coeff[1], 1) == 0);

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;
  mpz_t * Tr = (mpz_t *) malloc(sizeof(mpz_t) * (t - 1));
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_init(Tr[i]);
  }
  ret = parse_cs_mpzs(Tr, &tmp, str);
  ASSERT(ret == (int)(t - 1));
#ifndef NDEBUG
  for (unsigned int i = 0; i < t - 1; i++) {
    ASSERT(mpz_cmp_ui(Tr[i], 0) == 0);
  }
#endif

  if (!ret || tmp[0] != ':')
    return 0;
  str = tmp + 1;

  double log = atof(str);

  ideal_pr_set_element(ideal, r, log, t);

  mpz_poly_clear(h);
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(Tr[i]);
  }
  free(Tr);
  return ret;
}

void read_factor_base(FILE * file, factor_base_t * fb, uint64_t * fbb,
    unsigned int * lpb_bit, cado_poly_srcptr f, double * log2_base,
    unsigned int t)
{
  for (unsigned int k = 0; k < (unsigned int) f->nb_polys; k++) {
    int ret;
    mpz_t lpb;
    mpz_init(lpb);
    mpz_ui_pow_ui(lpb, 2, (unsigned long) lpb_bit[k]);

    int size_line = 1024;
    char line [size_line];

    uint64_t fbb_tmp;
    ret = fscanf(file, "fbb:%" PRIu64 "\n", &fbb_tmp);
    ASSERT_ALWAYS(ret == 1);
    ASSERT(fbb[k] <= fbb_tmp);

    unsigned int lpb_bit_tmp;
    ret = fscanf(file, "lpb:%u\n", &lpb_bit_tmp);
    ASSERT_ALWAYS(ret == 1);
    ASSERT(lpb_bit[k] <= lpb_bit_tmp);

    int n;
    ret = fscanf(file, "deg:%d\n", &n);
    ASSERT_ALWAYS(ret == 1);

    mpz_poly f_tmp;
    mpz_poly_init(f_tmp, n);
    if (fgets(line, size_line, file) == NULL) {
      return;
    }
    ASSERT_ALWAYS(parse_line_mpz_poly(f_tmp, line, n) == n + 1);
    ASSERT(mpz_poly_cmp(f->pols[k], f_tmp) == 0);
    mpz_poly_clear(f_tmp);

    unsigned int t_tmp;
    ret = fscanf(file, "t:%u\n", &t_tmp);
    ASSERT_ALWAYS(ret == 1);
    ASSERT(t == t_tmp);

    uint64_t max_number_element_1, max_number_element_u, max_number_element_pr;
    ret = fscanf(file, "%" PRIu64 ":%" PRIu64 ":%" PRIu64 "\n",
        &max_number_element_1, &max_number_element_u, &max_number_element_pr);
    ASSERT_ALWAYS(ret == 3);

    factor_base_init(fb[k], max_number_element_1, max_number_element_u,
        max_number_element_pr);

    uint64_t number_element_1 = 0, number_element_u = 0, number_element_pr = 0;

    ideal_1_t ideal_1;
    ideal_1_init(ideal_1);
    for (uint64_t i = 0; i < max_number_element_1; i++) {
      if (fgets(line, size_line, file) == NULL) {
        return;
      }
      parse_ideal_1(ideal_1, line, t);
      ideal_1->log = ideal_1->log / log2_base[k];
      if (ideal_1->ideal->r <= fbb[k]) {
        factor_base_set_ideal_1(fb[k], number_element_1, ideal_1, t);
        number_element_1++;
      }
    }
    ideal_1_clear(ideal_1, t);

    for (uint64_t i = 0; i < max_number_element_u; i++) {
      if (fgets(line, size_line, file) == NULL) {
        return;
      }
      ideal_u_t ideal_u;
      ideal_u_init(ideal_u);
      parse_ideal_u(ideal_u, line, t);
      mpz_t limit;
      mpz_init(limit);
      mpz_set_d(limit, pow((double)ideal_u->ideal->r,
            (double)ideal_u->ideal->h->deg));
      ideal_u->log = ideal_u->log / log2_base[k];
      if (mpz_cmp(lpb, limit) >= 0 && ideal_u->ideal->r <= fbb[k]) {
        //TODO: Problem here.
        factor_base_set_ideal_u(fb[k], number_element_u, ideal_u, t);
        number_element_u++;
      }
      mpz_clear(limit);
      ideal_u_clear(ideal_u, t);
    }

    if (max_number_element_pr != 0) {
      ideal_pr_t ideal_pr;
      ideal_pr_init(ideal_pr);
      for (uint64_t i = 0; i < max_number_element_pr; i++) {
        if (fgets(line, size_line, file) == NULL) {
          return;
        }
        parse_ideal_pr(ideal_pr, line, t);
        ideal_pr->log = ideal_pr->log / log2_base[k];
        if (ideal_pr->ideal->r <= fbb[k]) {
          factor_base_set_ideal_pr(fb[k], number_element_pr, ideal_pr->ideal->r,
              t);
          number_element_pr++;
        }
      }
      ideal_pr_clear(ideal_pr, t);
    }

    factor_base_realloc(fb[k], number_element_1, number_element_u,
        number_element_pr);

    mpz_clear(lpb);

    ret = fscanf(file, "----------------------------------------\n");
    ASSERT_ALWAYS(ret == 0);
  }
}

#ifdef MAIN
/*
 * Write ideal in a file.
 */
void write_ideal(FILE * file, ideal_srcptr ideal)
{
  fprintf(file, "%d:", ideal->h->deg);
  fprintf(file, "%" PRIu64 ":", ideal->r);
  gmp_fprintf(file, "%Zd", ideal->h->coeff[0]);
  for (int i = 1; i <= ideal->h->deg; i++) {
    gmp_fprintf(file, ",%Zd", ideal->h->coeff[i]);
  }
  fprintf(file, ":");
}

/*
 * Write ideal_1 in a file.
 */
void write_ideal_1(FILE * file, ideal_1_srcptr ideal, unsigned int t) {
  write_ideal(file, ideal->ideal);
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd,", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd:", ideal->Tr[t - 2]);
  fprintf(file, "%f\n", ideal->log);
}

/*
 * Write ideal_u in a file.
 */
void write_ideal_u(FILE * file, ideal_u_srcptr ideal, unsigned int t) {
  write_ideal(file, ideal->ideal);
  for (int row = 0; row < ideal->ideal->h->deg - 1; row++) {
    for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
      gmp_fprintf(file, "%Zd,", ideal->Tr[row][col]);
    }
    gmp_fprintf(file, "%Zd,", ideal->Tr[row]
                [t - (unsigned int)ideal->ideal->h->deg - 1]);
  }
  for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
    gmp_fprintf(file, "%Zd,", ideal->Tr[ideal->ideal->h->deg - 1][col]);
  }
  gmp_fprintf(file, "%Zd:", ideal->Tr[ideal->ideal->h->deg - 1]
              [t - (unsigned int)ideal->ideal->h->deg - 1]);
  fprintf(file, "%f\n", ideal->log);
}

/*
 * Write ideal_pr in a file.
 */
void write_ideal_pr(FILE * file, ideal_pr_srcptr ideal, unsigned int t)
{
  write_ideal(file, ideal->ideal);
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd,", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd:", ideal->Tr[t - 2]);
  fprintf(file, "%f\n", ideal->log);
}

/*
 * Write factor base in a file.
 *
 * file: the file.
 * fb: the factor base.
 * f: polynomial that defines the number field.
 * fbb: factor base bound.
 * lpb: large prime bound.
 * t: dimension of the lattice.
 */
void export_factor_base(FILE * file, factor_base_t * fb, cado_poly_srcptr f,
    uint64_t * fbb, unsigned int * lpb, unsigned int t)
{
  for (unsigned int k = 0; k < (unsigned int)f->nb_polys; k++) {
    fprintf(file, "fbb:%" PRIu64 "\n", fbb[k]);
    fprintf(file, "lpb:%u\n", lpb[k]);
    fprintf(file, "deg:%d\n", f->pols[k]->deg);
    fprintf(file, "f:");
    for (int i = 0; i < f->pols[k]->deg; i++) {
      gmp_fprintf(file, "%Zd,", f->pols[k]->coeff[i]);
    }
    gmp_fprintf(file, "%Zd\n", f->pols[k]->coeff[f->pols[k]->deg]);
    fprintf(file, "t:%u\n", t);
    fprintf(file, "%" PRIu64 ":%" PRIu64 ":%" PRIu64 "\n",
        fb[k]->number_element_1, fb[k]->number_element_u,
        fb[k]->number_element_pr);

    for (uint64_t i = 0; i < fb[k]->number_element_1; i++) {
      write_ideal_1(file, fb[k]->factor_base_1[i], t);
    }

    for (uint64_t i = 0; i < fb[k]->number_element_u; i++) {
      write_ideal_u(file, fb[k]->factor_base_u[i], t);
    }

    for (uint64_t i = 0; i < fb[k]->number_element_pr; i++) {
      write_ideal_pr(file, fb[k]->factor_base_pr[i], t);
    }
    fprintf(file, "----------------------------------------\n");
  }
}

void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "t", "dimension of the lattice");
  param_list_decl_usage(pl, "poly", "path to the polynomial file");
  param_list_decl_usage(pl, "fbb", "factor base bounds");
  param_list_decl_usage(pl, "lpb", "large prime bounds");
  param_list_decl_usage(pl, "out", "path to file where factor bases are \
      stored");
}

/*
 * Initialise all the parameters. Frequently read from command line.
 *
 * f: array of mpz_poly that define the number fields.
 * fbb: factor base bounds.
 * fb: factor bases.
 * t: dimension of the lattice.
 * lpb: large prime bounds.
 * V: number of number fields.
 * p: characteristic of the finite field.
 * n: extension of the finite field.
 */
void initialise_parameters(int argc, char * argv[], cado_poly_ptr f,
    uint64_t ** fbb, factor_base_t ** fb, unsigned int * t, unsigned int ** lpb,
    char ** file_r)
{
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  FILE * fpl;
  char * argv0 = argv[0];

  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

    /* Could also be a file */
    if ((fpl = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, fpl, 0);
      fclose(fpl);
      argv++,argc--;
      continue;
    }

    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit (EXIT_FAILURE);
  }

  cado_poly_init(f);
  unsigned int size_path = 1024;
  char path [size_path];
  param_list_parse_string(pl, "poly", path, size_path);
  cado_poly_read(f, path);
  ASSERT(f->nb_polys >= 2);
  unsigned int V = (unsigned int) f->nb_polys;

  * fbb = malloc(sizeof(uint64_t) * (V));
  * fb = malloc(sizeof(factor_base_t) * (V));
  * lpb = malloc(sizeof(unsigned int) * (V));

  param_list_parse_uint64_list(pl, "fbb", * fbb, (size_t) V, ",");
  param_list_parse_uint_list(pl, "lpb", * lpb, (size_t) V, ",");

  for (unsigned int i = 0; i < V; i++) {
    factor_base_init((*fb)[i], (*fbb)[i], (*fbb)[i], (*fbb)[i]);
  }

  /*for (unsigned int i = 0; i < * V; i++) {*/
    /*ASSERT(mpz_cmp_ui((*lpb)[i], (*fbb)[i]) >= 0);*/
  /*}*/

  param_list_parse_uint(pl, "t", t);
  ASSERT(* t > 2);

  size_t size = 1024;
  * file_r = malloc(sizeof(char) * size);
  param_list_parse_string(pl, "out", * file_r, size);

  param_list_clear(pl);
}

int main(int argc, char ** argv)
{
  cado_poly f;
  uint64_t * fbb;
  unsigned int t;
  unsigned int * lpb;
  factor_base_t * fb;
  char * file_r;

  initialise_parameters(argc, argv, f, &fbb, &fb, &t, &lpb, &file_r);

#ifndef TIME_MAKEFB
  double sec = seconds();
#endif // TIME_MAKEFB

  makefb(fb, f, fbb, t, lpb);

  //5 because name of the file is p,n,V.
  FILE * file;
  file = fopen (file_r, "w+");

#ifndef PRINT_INFO
  printf("# Write the factor bases in the file %s.\n", file_r);
#endif // PRINT_INFO

  export_factor_base(file, fb, f, fbb, lpb, t);

  fclose(file);
#ifndef TIME_MAKEFB
  printf("# Time to build makefb: %fs.\n", seconds() - sec);
#endif // TIME_MAKEFB

  for (unsigned int i = 0; i < (unsigned int)f->nb_polys; i++) {
    factor_base_clear(fb[i], t);
  }

  cado_poly_clear(f);
  free(fbb);
  free(lpb);
  free(fb);
  free(file_r);

  return 0;
}
#endif // MAIN
