/* This standalone program checks matrix products produced by lingen
   (as stored in the cp/ subdirectory).

   Example: check -p ... -dim 13 -k 3 pi.3.2.127 pi.3.127.251 pi.2.2.251

   Optional arguments: [-seed xxx] [-v]
*/

#include "cado.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gmp.h>
#include <cassert>
#include <sys/types.h>
#include <unistd.h>
#include "macros.h"
#include "params.h"
#include "cxx_mpz.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

// #define FORMAT 0 /* matrices of polynomials */
#define FORMAT 1 /* polynomials of matrices */

/* define WARNING to get warning for non-zero padding coefficients */
// #define WARNING

struct {
    unsigned int m,n;
} bw_parameters;

cxx_mpz p;                /* prime modulus */
unsigned long lingen_p; /* number of limbs per coefficient */
unsigned long dim = 0;  /* matrix dimension */
int k = 1;              /* matrix is cut in k x k submatrices */
gmp_randstate_t state;
int verbose = 0;
unsigned long seed;

typedef struct
{
  unsigned long n;   /* matrix dimension */
  unsigned long k;   /* matrix is cut into k x k submatrices */
  unsigned long deg; /* degree of coefficients */
  mpz_t **coeff;     /* coeff[i][j] is P(x=v) mod p */
} matrix_struct;
typedef matrix_struct matrix[1];

void
init_matrix (matrix M, unsigned long n, unsigned long k, unsigned long deg)
{
  M->n = n;
  M->k = k;
  M->deg = deg;
  M->coeff = (mpz_t **) malloc (n * sizeof (mpz_t*));
  for (unsigned long i = 0; i < n; i++)
    {
      M->coeff[i] = (mpz_t *) malloc (n * sizeof (mpz_t));
      for (unsigned long j = 0; j < n; j++)
        mpz_init (M->coeff[i][j]);
    }
}

/* initializes a vector of n numbers */
mpz_t *
init_vector (unsigned long n)
{
  mpz_t *u;
  u = (mpz_t *) malloc (n * sizeof (mpz_t));
  for (unsigned long i = 0; i < n; i++)
    mpz_init (u[i]);
  return u;
}

/* return a vector of n random numbers mod p */
void
random_vector (mpz_t *u, unsigned long n)
{
  for (unsigned long i = 0; i < n; i++)
    mpz_urandomm (u[i], state, p);
  // mpz_set_ui (u[i], i == 0);
}

void
clear_vector (mpz_t *u, unsigned long n)
{
  for (unsigned long i = 0; i < n; i++)
    mpz_clear (u[i]);
  free (u);
}

void
clear_matrix (matrix M)
{
  unsigned long n = M->n;
  for (unsigned long i = 0; i < n; i++)
    {
      for (unsigned long j = 0; j < n; j++)
        mpz_clear (M->coeff[i][j]);
      free (M->coeff[i]);
    }
  free (M->coeff);
}

void
read_coeff (FILE *fp, mpz_t c)
{
  mp_size_t n;
  size_t size;
  mpz_realloc2 (c, lingen_p * mp_bits_per_limb);
  size = fread (c->_mp_d, sizeof (mp_limb_t), lingen_p, fp);
  ASSERT_ALWAYS (size == lingen_p);
  n = lingen_p;
  while (n > 0 && c->_mp_d[n-1] == 0)
    n--;
  c->_mp_size = n;
}

/* return non-zero if next lingen_p limbs read are zero */
int
read_zero (FILE *fp)
{
  mp_limb_t buf[1];
  size_t size;
  for (unsigned long i = 0; i < lingen_p; i++)
    {
      size = fread (buf, sizeof (mp_limb_t), 1, fp);
      ASSERT_ALWAYS (size == 1);
    }
  return true;
}

unsigned long
read_degree (const char *s)
{
  char t[1024];
  strcpy (t, s);
  sprintf (t + strlen (t), ".aux");
  FILE *fp;
  fp = fopen (t, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "Error, unable to read file %s\n", t);
      exit (1);
    }
  int format;
  int ret = fscanf (fp, "format %d\n", &format);
  ASSERT_ALWAYS (ret == 1);
  unsigned long ncoeff;
  ret = fscanf (fp, "%lu\n", &ncoeff);
  ASSERT_ALWAYS (ret == 1);
  fclose (fp);
  return ncoeff - 1;
}

#if FORMAT == 0
void
read_submatrix (matrix M, FILE *fp, unsigned long starti, unsigned long ni,
                unsigned long startj, unsigned long nj, unsigned long m,
                mpz_t x)
{
  unsigned long deg = M->deg;
  unsigned long i, j, k;
  mpz_t *Buf;
  Buf = init_vector (M->deg + 1);
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      {
        if (i < ni && j < nj)
          {
            for (k = 0; k <= deg; k++)
              read_coeff (fp, Buf[k]);
            eval_pol (M->coeff[i][j], Buf, x, M->deg);
          }
        else
          {
            int ok = 1;
            for (k = 0; k <= deg; k++)
              ok = ok &&  read_zero (fp);
#ifdef WARNING
            if (ok == 0)
              fprintf (stderr, "Warning, padding coefficient i=%lu,j=%lu is not zero\n", i, j);
#endif
          }
      }
  clear_vector (Buf);
}
#else
void
read_submatrix (matrix M, FILE *fp, unsigned long starti, unsigned long ni,
                unsigned long startj, unsigned long nj, unsigned long m,
                mpz_t x)
{
  unsigned long deg = M->deg;
  unsigned long i, j, k;
  char *T;
  mpz_t x_power_k, tmp;
  T = (char *) malloc (m * m);
  memset (T, 0, m * m);
  mpz_init (tmp);
  mpz_init_set_ui (x_power_k, 1);
  for (k = 0; k <= deg; k++)
    {
      /* invariant: x_power_k = x^k mod p */
      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
          {
            if (i < ni && j < nj)
              {
                // read_coeff (fp, M->coeff[starti + i][startj + j][k]);
                if (k == 0)
                  read_coeff (fp, M->coeff[starti + i][startj + j]);
                else
                  {
                    read_coeff (fp, tmp);
                    mpz_mul (tmp, tmp, x_power_k);
                    mpz_add (tmp, tmp, M->coeff[starti + i][startj + j]);
                    mpz_mod (M->coeff[starti + i][startj + j], tmp, p);
                  }
              }
            else
              {
#ifndef WARNING
                read_zero (fp);
#else
                int ok = read_zero (fp);
                if (ok == 0 && T[m+i+j] == 0)
                  {
                    fprintf (stderr, "Warning, padding coefficient %lu,%lu is not zero\n", i, j);
                    T[m+i+j] = 1;
                  }
#endif
              }
          }
      mpz_mul (x_power_k, x_power_k, x);
      mpz_mod (x_power_k, x_power_k, p);
    }
  mpz_clear (x_power_k);
  mpz_clear (tmp);
  free (T);
}
#endif

/* read a matrix of dimension n, divided into kxk submatrices */
void
read_matrix (matrix M, const char *s, unsigned long n, unsigned long k, mpz_t x)
{
  unsigned long remi = n; /* it remains remi rows to read */
  unsigned long starti = 0;
  unsigned long deg = read_degree (s);
  init_matrix (M, n, k, deg);
  unsigned long m = (n + k - 1) / k; /* size of submatrices */
  for (unsigned long i = 0; i < k; i++)
    {
      /* it remains remi rows to read, for (k-i) submatrices */
      unsigned long ni = (remi + (k - i) - 1) / (k - i); /* ceil (remi/(k-i)) */
      unsigned long remj = n; /* it remains remi columns to read */
      unsigned long startj = 0;
      for (unsigned long j = 0; j < k; j++)
        {
          unsigned long nj = (remj + (k - j) - 1) / (k - j);
          FILE *fp;
          char t[1024];
          int nij = M->k * i + j;
          strcpy (t, s);
          if (k > 1) {
              sprintf (t + strlen (t), ".%d.data", nij);
          } else {
              sprintf (t + strlen (t), ".single.data");
          }
          if (verbose) {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
	    printf ("Reading %lux%lu matrix %s\n", ni, nj, t);
          }
          fp = fopen (t, "r");
	  if (fp == NULL)
	    {
	      fprintf (stderr, "Error, unable to read file %s\n", t);
	      exit (1);
	    }
          read_submatrix (M, fp, starti, ni, startj, nj, m, x);
          // assert (feof (fp));
          fclose (fp);
          startj += nj;
          remj -= nj;
        }
      starti += ni;
      remi -= ni;
    }
}

/* w <- P(x=v) mod p */
void
eval_pol (mpz_t w, mpz_t *P, mpz_t v, unsigned long deg)
{
  mpz_set (w, P[deg]);
  for (unsigned long i = deg; i-- > 0;)
    {
      mpz_mul (w, w, v);
      mpz_add (w, w, P[i]);
      mpz_mod (w, w, p);
    }
}

/* w <- v*M evaluated at x and modulo p */
void
mul_left (mpz_t *w, mpz_t *v, matrix M)
{
  unsigned long n = M->n;
  mpz_t tmp;
  mpz_init (tmp);
  for (unsigned long j = 0; j < n; j++)
    {
      mpz_set_ui (w[j], 0);
      for (unsigned long i = 0; i < n; i++)
        {
          /* w[j] += v[i]*M[i,j] */
          mpz_mul (tmp, v[i], M->coeff[i][j]);
          mpz_add (w[j], w[j], tmp);
        }
      mpz_mod (w[j], w[j], p);
    }
  mpz_clear (tmp);
}

/* w <- M*v evaluated */
void
mul_right (mpz_t *w, matrix M, mpz_t *v)
{
  unsigned long n = M->n;
  mpz_t tmp;
  mpz_init (tmp);
  for (unsigned long i = 0; i < n; i++)
    {
      mpz_set_ui (w[i], 0);
      for (unsigned long j = 0; j < n; j++)
        {
          /* w[i] += M[i,j]*v[j] */
          mpz_mul (tmp, M->coeff[i][j], v[j]);
          mpz_add (w[i], w[i], tmp);
        }
      mpz_mod (w[i], w[i], p);
    }
  mpz_clear (tmp);
}

void
scalar_product (mpz_t res, mpz_t *u, mpz_t *v, unsigned long n)
{
  mpz_set_ui (res, 0);
  for (unsigned long i = 0; i < n; i++)
    mpz_addmul (res, u[i], v[i]);
  mpz_mod (res, res, p);
}

void
print_vector (mpz_t *u, unsigned long n)
{
  for (unsigned long i = 0; i < n; i++)
    gmp_printf ("%Zd ", u[i]);
  printf ("\n");
}

/* print a 0 for zero coefficients, otherwise 1 */
void
print_matrix (matrix M)
{
  unsigned long n = M->n;
  for (unsigned long i = 0; i < n; i++)
    {
      for (unsigned long j = 0; j < n; j++)
        if (mpz_cmp_ui (M->coeff[i][j], 0) == 0)
          printf ("0 ");
        else
          printf ("1 ");
      printf ("\n");
    }
}

void declare_usage(cxx_param_list & pl)
{
  param_list_usage_header(pl,
      "Usage: lingen_verify_checkpoints [options] -- [list of file names]\n"
      "\n"
      "The list of file names can have one of the following formats:\n"
      " - pi0 pi1 pi2 : check that pi0*pi1==pi2\n"
      " - E pi : check that E*pi=O(x^length(E))\n"
      "Options are as follows.\n"
      );
  param_list_decl_usage(pl, "prime", "characteristic of the base field");
  param_list_decl_usage(pl, "mpi", "mpi geometry");
  param_list_decl_usage(pl, "seed", "random seed");
  param_list_decl_usage(pl, "m", "block Wiedemann parameter m");
  param_list_decl_usage(pl, "n", "block Wiedemann parameter n");
  param_list_decl_usage(pl, "v", "More verbose output");
}

void lookup_parameters(cxx_param_list & pl)
{
  param_list_lookup_string(pl, "prime");
  param_list_lookup_string(pl, "mpi");
  param_list_lookup_string(pl, "seed");
}

int do_check_pi(const char * pi_left_filename, const char * pi_right_filename, const char * pi_filename)
{
  int ret;
  matrix Mab, Mbc, Mac;
  unsigned long n = dim;

  mpz_t *u = init_vector (n);
  random_vector (u, n);
  mpz_t *v = init_vector (n);
  random_vector (v, n);
  mpz_t x;
  mpz_init (x);
  mpz_urandomm (x, state, p); /* random variable */
  printf ("Using seed %lu\n", seed);
  if (verbose)
    gmp_printf ("x=%Zd\n", x);

  mpz_t *u_times_piab, *pibc_times_v, *piac_times_v;
#ifdef HAVE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef HAVE_OPENMP
#pragma omp section
#endif
    {
      u_times_piab = init_vector (n);
      read_matrix (Mab, pi_left_filename, dim, k, x);
      mul_left (u_times_piab, u, Mab);
      clear_matrix (Mab);
    }

#ifdef HAVE_OPENMP
#pragma omp section
#endif
    {
      pibc_times_v = init_vector (n);
      read_matrix (Mbc, pi_right_filename, dim, k, x);
      mul_right (pibc_times_v, Mbc, v);
      clear_matrix (Mbc);
    }

#ifdef HAVE_OPENMP
#pragma omp section
#endif
    {
      piac_times_v = init_vector (n);
      read_matrix (Mac, pi_filename, dim, k, x);
      mul_right (piac_times_v, Mac, v);
      clear_matrix (Mac);
    }
  }

  mpz_t res_left, res_right;
  mpz_init (res_left);
  mpz_init (res_right);

  scalar_product (res_left, u_times_piab, pibc_times_v, n);
  scalar_product (res_right, u, piac_times_v, n);
  if (mpz_cmp (res_left, res_right) != 0)
  {
    fprintf (stderr, "Error, results differ\n");
    gmp_printf ("res_left  = %Zd\n", res_left);
    gmp_printf ("res_right = %Zd\n", res_right);
    ret = 0;
  }
  else
  {
    printf ("Check is ok\n");
    ret = 1;
  }

  mpz_clear (res_left);
  mpz_clear (res_right);
  mpz_clear (x);
  clear_vector (u, n);
  clear_vector (v, n);
  clear_vector (u_times_piab, n);
  clear_vector (pibc_times_v, n);
  clear_vector (piac_times_v, n);

  return ret;
}

int
main (int argc, char *argv[])
{
  cxx_param_list pl;

  seed = getpid ();

  declare_usage(pl);

  const char * argv0 = argv[0];

  param_list_configure_switch(pl, "-v", &verbose);
  for(argc--, argv++ ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    if (strcmp(argv[0], "--") == 0) {
      argc--,argv++;
      break;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  if (!param_list_parse_mpz(pl, "prime", (mpz_ptr) p)) {
    fprintf(stderr, "Missing parameter: prime\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }
  int mpi_dims[2]={0,0};
  if (param_list_parse_int_and_int(pl, "mpi", mpi_dims, "x")) {
    ASSERT_ALWAYS(mpi_dims[0] == mpi_dims[0]);
    k = mpi_dims[0];
  }
  param_list_parse_ulong(pl, "seed", &seed);
  if (!param_list_parse_uint(pl, "m", &bw_parameters.m)) {
    fprintf(stderr, "Missing parameter: m\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }
  if (!param_list_parse_uint(pl, "n", &bw_parameters.m)) {
    fprintf(stderr, "Missing parameter: n\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }
  if (param_list_warn_unused(pl)) {
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  gmp_randinit_default (state);
  gmp_randseed_ui (state, seed);

  lingen_p = mpz_size (p);
  if (verbose)
    printf ("Number of limbs: %zu\n", lingen_p);

  int ret = 0;

  if (argc == 3) {
    ret = do_check_pi(argv[0], argv[1], argv[2]);
  } else if (argc == 2) {
    fprintf(stderr, "not implemented\n");
    ret = 1;
  } else {
    fprintf(stderr, "bad argument list\n");
    exit(EXIT_FAILURE);
  }

  gmp_randclear (state);

  return ret ? EXIT_SUCCESS : EXIT_FAILURE;
}
/* vim :sw=2 sta et: */
