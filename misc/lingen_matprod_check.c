/* This standalone program checks matrix products produced by lingen
   (as stored in the cp/ subdirectory).

   Example: check -p ... -dim 13 -k 3 pi.3.2.127 pi.3.127.251 pi.2.2.251

   Optional arguments: [-seed xxx] [-v]
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

// #define FORMAT 0 /* matrices of polynomials */
#define FORMAT 1 /* polynomials of matrices */

/* define WARNING to get warning for non-zero padding coefficients */
// #define WARNING

mpz_t p;                /* prime modulus */
unsigned long lingen_p; /* number of limbs per coefficient */
unsigned long dim = 0;  /* matrix dimension */
int k = 0;              /* matrix is cut in k x k submatrices */
gmp_randstate_t state;
int verbose = 0;

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
  M->coeff = malloc (n * sizeof (mpz_t*));
  for (unsigned long i = 0; i < n; i++)
    {
      M->coeff[i] = malloc (n * sizeof (mpz_t));
      for (unsigned long j = 0; j < n; j++)
        mpz_init (M->coeff[i][j]);
    }
}

/* initializes a vector of n numbers */
mpz_t *
init_vector (unsigned long n)
{
  mpz_t *u;
  u = malloc (n * sizeof (mpz_t));
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
  assert (size == lingen_p);
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
  int ok;
  for (unsigned long i = 0; i < lingen_p; i++)
    {
      size = fread (buf, sizeof (mp_limb_t), 1, fp);
      assert (size == 1);
      ok = ok && buf == 0;
    }
  return ok;
}

unsigned long
read_degree (char *s)
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
  assert (ret == 1);
  unsigned long ncoeff;
  ret = fscanf (fp, "%lu\n", &ncoeff);
  assert (ret == 1);
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
  T = malloc (m * m);
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
read_matrix (matrix M, char *s, unsigned long n, unsigned long k, mpz_t x)
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
          sprintf (t + strlen (t), ".%d.data", nij);
          if (verbose)
#pragma omp critical
	    printf ("Reading %lux%lu matrix %s\n", ni, nj, t);
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

int
main (int argc, char *argv[])
{
  matrix Mab, Mbc, Mac;
  unsigned long seed = getpid ();

  mpz_init (p);

  while (argc >= 2 && argv[1][0] == '-')
    {
      if (argc >= 3 && strcmp (argv[1], "-p") == 0)
        {
          mpz_set_str (p, argv[2], 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-dim") == 0)
        {
          dim = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-k") == 0)
        {
          k = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-seed") == 0)
        {
          seed = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 2 && strcmp (argv[1], "-v") == 0)
        {
          verbose ++;
          argc -= 1;
          argv += 1;
        }
    }
  
  if (mpz_cmp_ui (p, 0) == 0)
    {
      fprintf (stderr, "Error, missing -p argument\n");
      exit (1);
    }

  if (dim == 0)
    {
      fprintf (stderr, "Error, missing -dim argument\n");
      exit (1);
    }

  if (k == 0)
    {
      fprintf (stderr, "Error, missing -k argument\n");
      exit (1);
    }

  gmp_randinit_default (state);
  gmp_randseed_ui (state, seed);

  lingen_p = mpz_size (p);
  if (verbose)
    printf ("Number of limbs: %zu\n", lingen_p);

  assert (argc == 4);

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
#pragma omp parallel sections
  {
    #pragma omp section
    {
      u_times_piab = init_vector (n);
      read_matrix (Mab, argv[1], dim, k, x);
      mul_left (u_times_piab, u, Mab);
      clear_matrix (Mab);
    }

    #pragma omp section
    {
      pibc_times_v = init_vector (n);
      read_matrix (Mbc, argv[2], dim, k, x);
      mul_right (pibc_times_v, Mbc, v);
      clear_matrix (Mbc);
    }

    #pragma omp section
    {
      piac_times_v = init_vector (n);
      read_matrix (Mac, argv[3], dim, k, x);
      mul_right (piac_times_v, Mac, v);
      clear_matrix (Mac);
    }
  }

  mpz_t res_left, res_right;
  mpz_init (res_left);
  mpz_init (res_right);

  scalar_product (res_left, u_times_piab, pibc_times_v, n);
  scalar_product (res_right, u, piac_times_v, n);
  int ret;
  if (mpz_cmp (res_left, res_right) != 0)
      {
        fprintf (stderr, "Error, results differ\n");
        gmp_printf ("res_left  = %Zd\n", res_left);
        gmp_printf ("res_right = %Zd\n", res_right);
        ret = 1;
      }
  else
    {
      printf ("Check is ok\n");
      ret = 0; /* return 0 if ok */
    }

  mpz_clear (res_left);
  mpz_clear (res_right);
  mpz_clear (x);
  clear_vector (u, n);
  clear_vector (v, n);
  clear_vector (u_times_piab, n);
  clear_vector (pibc_times_v, n);
  clear_vector (piac_times_v, n);
  mpz_clear (p);
  gmp_randclear (state);

  return ret;
}
