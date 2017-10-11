#include "cado.h"
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include "utils.h"
#include "polynomials.h"
#include "sieving_bound.h"
#include "polyselect/auxiliary.h"

static void mpz_poly_coeff_fprintf(FILE * file, mpz_poly_srcptr f)
{
  for (int i = 0; i < f->deg; i++) {
    gmp_fprintf(file, "%Zd,", f->coeff[i]);
  }
  gmp_fprintf(file, "%Zd\n", f->coeff[f->deg]);
}

void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "p", "prime number");
  param_list_decl_usage(pl, "n", "degree of the extension");
  param_list_decl_usage(pl, "coeff", "range of coefficients of g");
  param_list_decl_usage(pl, "nb_times", "number of tested polynomials");
  param_list_decl_usage(pl, "weight",
      "weight for alphas of when sum alphas");
  param_list_decl_usage(pl, "type", "0 if classical, 1 if special-q");
  param_list_decl_usage(pl, "q", "a typical q of a special-q (for type 1)");
  param_list_decl_usage(pl, "t", "dimension of the lattice (for type 1)");
  param_list_decl_usage(pl, "h", "polynomial to build f0 and f1 (for type 1)");
  param_list_decl_usage(pl, "out", "path to the output file");
  param_list_decl_usage(pl, "c_tol", "tolerance around c (for type 1)");
  param_list_decl_usage(pl, "gal", "galois of order 2");
}

void initialise_parameters(int argc, char * argv[], mpz_ptr p, mpz_poly_ptr h,
    int ** coeff, unsigned int * q, unsigned int * t,
    unsigned int * nb_times, double ** weight,
    unsigned int * n, unsigned int * type, FILE ** outstd, int * c_tol,
    unsigned int * h_set, unsigned int * gal)
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

  * type = 1;
  param_list_parse_uint(pl, "type", type);
  ASSERT(* type < 2);

  param_list_parse_mpz(pl, "p", p);
  param_list_parse_uint(pl, "n", n);

  * nb_times = 100;
  param_list_parse_uint(pl, "nb_times", nb_times);
  ASSERT(* nb_times > 0);

  * weight = (double *) malloc(2 * sizeof(double));
  (*weight)[0] = 0.5;
  (*weight)[1] = 0.5;
  param_list_parse_double_and_double(pl, "weight", * weight, ",");
  ASSERT((* weight)[0] >= 0.0);
  ASSERT((* weight)[1] >= 0.0);
  ASSERT((* weight)[0] + (* weight)[1] == 1.0);

  * coeff = (int *) malloc(2 * sizeof(int));
  param_list_parse_int_and_int(pl, "coeff", * coeff, ",");
  ASSERT((*coeff)[0] < (*coeff)[1]);
  ASSERT((*coeff)[0] < 0);
  ASSERT((*coeff)[1] > 0);

  if (* type == 1) {
    param_list_parse_mpz_poly(pl, "h", h);
    if (h->deg < (int) * n) {
      ASSERT(h->deg == 0);

      * h_set = 0;
      mpz_poly_realloc(h, (int)* n + 1);
      h->deg = * n;
      ASSERT(h->alloc == (int)* n + 1);
    } else {
      * h_set = 1;
      ASSERT(* n % ((unsigned int )h->deg) == 0);
    }

    param_list_parse_uint(pl, "q", q);

    param_list_parse_uint(pl, "t", t);
    ASSERT(* t > 2);

    * c_tol = 20;
    param_list_parse_int(pl, "c_tol", c_tol);
    ASSERT(* c_tol > 0);

    * gal = 0;
    param_list_parse_uint(pl, "gal", gal);
    ASSERT(* gal == 0 || * gal == 2);

    if (* h_set == 1 && * gal == 2) {
      ASSERT(mpz_poly_is_reciprocal(h));
    }
  }

  unsigned int size_path = 1024;
  char path [size_path];
  path[0] = '\0';
  param_list_parse_string(pl, "out", path, size_path);
  if (path[0] == '\0') {
    * outstd = stdout;
  } else {
    * outstd = fopen(path, "w");
  }


  param_list_clear(pl);
}

int main(int argc, char * argv[])
{
  mpz_t p;
  mpz_init(p);
  unsigned int n;
  mpz_poly h;
  mpz_poly_init(h, -1);
  unsigned int q;
  int * coeff;
  unsigned int t;
  unsigned int nb_times;
  double * weight;
  unsigned int type;
  FILE * outstd;
  int c_tol;
  unsigned int h_set;
  unsigned int gal;

  initialise_parameters(argc, argv, p, h, &coeff, &q, &t, &nb_times,
      &weight, &n, &type, &outstd, &c_tol, &h_set, &gal);

  mpz_poly f0;
  mpz_poly_init(f0, h->deg);
  mpz_poly f1;
  mpz_poly_init(f1, h->deg);
  mpz_poly g;
  mpz_poly_init(g, h->deg);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, time(NULL));
  srand(time(NULL));

  mpz_t a;
  mpz_init(a);
  mpz_t b;
  mpz_init(b);
  mpz_t c;
  mpz_init(c);

  double epsilon = 0.0;

  gmp_fprintf(outstd, "# p^n: %Zd^%u\n", p, n);
  gmp_fprintf(outstd, "n: %Zd\n", p);
  //TODO: check if it is needed or not for cado_poly.
  fprintf(outstd, "skew: 1.000\n");

  if (type == 0) {
    function_classical(f0, f1, p, n, coeff, nb_times, weight, state);
  } else if (type == 1) {
    epsilon = function_special_q(f0, f1, g, a, b, c, p, h, coeff, q, t,
        nb_times, weight, c_tol, state, h_set, gal);
  }

#ifndef NDEBUG
  mpz_t res;
  mpz_init(res);
  mpz_t tmp;
  mpz_init(tmp);
  mpz_poly_resultant(res, f0, f1);
  mpz_pow_ui(tmp, p, n);
  mpz_mod(tmp, res, tmp);
  ASSERT(mpz_cmp_ui(tmp, 0) == 0);
  gmp_fprintf(outstd, "# res(f0, f1) = %Zd\n", res);
  mpz_pow_ui(tmp, p, n);
  mpz_divexact(tmp, res, tmp);
  gmp_fprintf(outstd, "# Excess: %Zd\n", tmp);
  mpz_clear(res);
  mpz_clear(tmp);
#endif // NDEBUG

  fprintf(outstd, "# f0 = ");
  mpz_poly_fprintf(outstd, f0);
  fprintf(outstd, "poly0: ");
  mpz_poly_coeff_fprintf(outstd, f0);
  fprintf(outstd, "# alpha0: %f\n", get_alpha(f0, ALPHA_BOUND));
  fprintf(outstd, "# f1 = ");
  mpz_poly_fprintf(outstd, f1);
  fprintf(outstd, "poly1: ");
  mpz_poly_coeff_fprintf(outstd, f1);
  fprintf(outstd, "# alpha1: %f\n", get_alpha(f1, ALPHA_BOUND));

  if (type == 1) {
    fprintf(outstd, "# t = %u\n", t);
    fprintf(outstd, "# log2(q) = %u\n", q);
    fprintf(outstd, "# epsilon = %f\n", epsilon);
    gmp_fprintf(outstd, "# a = %Zd; b = %Zd; c = %Zd\n", a, b, c);
    fprintf(outstd, "# g = ");
    mpz_poly_fprintf(outstd, g);
    fprintf(outstd, "# h = ");
    mpz_poly_fprintf(outstd, h);
  }

  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_poly_clear(f0);
  mpz_poly_clear(f1);

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);
  mpz_clear(p);

  gmp_randclear(state);

  fclose(outstd);

  return 0;
}
