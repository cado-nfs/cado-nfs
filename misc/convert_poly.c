/* Convert between different polynomial formats
   (CADO, Franke-Kleinjung, GGNFS, CWI).
   First version written by Paul Zimmermann, 20030804.
   Rewritten to allow different input/output formats, 20071220.

   Call with --help or any invalid option for usage information.
 */

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define MAX_DEGREE 6

/************************** input routines ***********************************/

static int
readline (char *s)
{
  int c;

  while ((c = getchar ()) != '\n' && c != EOF)
    *s++ = c;
  *s = '\0';
  return c;
}

static void
read_franke (mpz_t N, mpz_t *X, mpz_t *Y, mpz_t M)
{
  int c, i;

  if (mpz_inp_str (N, stdin, 0) == 0)
    {
      fprintf (stderr, "Error while reading N\n");
      exit (1);
    }
  c = getchar ();
  if (c != '\n')
    {
      fprintf (stderr, "Error: end of line expected after N\n");
      exit (1);
    }

  /* read comment lines */
  while ((c = getchar ()) == '#')
    {
      while ((c = getchar ()) != '\n');
    }

  /* now read X lines */
  if (c != 'X')
    {
      fprintf (stderr, "Error: Xnnn expected\n");
      exit (1);
    }
  do
    {
      if (scanf ("%d ", &i) != 1)
        {
          fprintf (stderr, "Error: Xnnn expected\n");
          exit (1);
        }
      if (i > MAX_DEGREE)
        {
          fprintf (stderr, "Error: too large degree X%d expected\n", i);
          exit (1);
        }
      mpz_inp_str (X[i], stdin, 0);
      if (getchar () != '\n')
        {
          fprintf (stderr, "Error: end of line expected after X%d\n", i);
          exit (1);
        }
      c = getchar ();
    }
  while (c == 'X');

  /* now read Y lines */
  if (c != 'Y')
    {
      fprintf (stderr, "Error: Ynnn expected\n");
      exit (1);
    }
  do
    {
      if (scanf ("%d ", &i) != 1)
        {
          fprintf (stderr, "Error: Ynnn expected\n");
          exit (1);
        }
      if (i > MAX_DEGREE)
        {
          fprintf (stderr, "Error: too large degree Y%d expected\n", i);
          exit (1);
        }
      mpz_inp_str (Y[i], stdin, 0);
      if (getchar () != '\n')
        {
          fprintf (stderr, "Error: end of line expected after Y%d\n", i);
          exit (1);
        }
      c = getchar ();
    }
  while (c == 'Y');

  if (c != 'M' || getchar () != ' ')
    {
      fprintf (stderr, "Error: M expected\n");
      exit (1);
    }
  mpz_inp_str (M, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after M\n");
      exit (1);
    }
}

#define MAX_BUF 4096

static void
read_ggnfs (mpz_t N, mpz_t *X, mpz_t *Y, mpz_t M)
{
  int i, ret;
  char s[MAX_BUF]; /* input buffer */

  while (feof (stdin) == 0)
    {
      ret = readline (s);
      if (ret == EOF)
        break;
      if (strlen (s) + 1 >= MAX_BUF)
        {
          fprintf (stderr, "Error, buffer overflow\n");
          exit (1);
        }
      if (strncmp (s, "n:", 2) == 0) /* n: input number */
        {
          if (mpz_set_str (N, s + 2, 0) != 0)
            {
              fprintf (stderr, "Error while reading N");
              exit (1);
            }
        }
      else if (sscanf (s, "c%d:", &i) == 1) /* ci: coeff of degree i */
        {
          if (i > MAX_DEGREE)
            {
              fprintf (stderr, "Error, too large degree %d\n", i);
              exit (1);
            }
          if (mpz_set_str (X[i], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading X[%d]", i);
              exit (1);
            }
        }
      else if (strncmp (s, "Y1:", 3) == 0)
        {
          if (mpz_set_str (Y[1], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading Y1");
              exit (1);
            }
        }
      else if (strncmp (s, "Y0:", 3) == 0)
        {
          if (mpz_set_str (Y[0], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading Y0");
              exit (1);
            }
        }
      else if (strncmp (s, "# M ", 4) == 0 || strncmp (s, "m:", 2) == 0)
        {
          if (mpz_set_str (M, s + 2 + (s[0] == '#'), 0) != 0)
            {
              fprintf (stderr, "Error while reading M or m:\n");
              fprintf (stderr, "%s\n", s);
              exit (1);
            }
        }
    }
}

static void
read_msieve (mpz_t N, mpz_t *X, mpz_t *Y)
{
  int i, ret;
  char s[MAX_BUF]; /* input buffer */

  while (feof (stdin) == 0)
    {
      ret = readline (s);
      if (ret == EOF)
        break;
      if (strlen (s) + 1 >= MAX_BUF)
        {
          fprintf (stderr, "Error, buffer overflow\n");
          exit (1);
        }
      if (strncmp (s, "N ", 2) == 0) /* N <input number> */
        {
          if (mpz_set_str (N, s + 2, 0) != 0)
            {
              fprintf (stderr, "Error while reading N");
              exit (1);
            }
        }
      else if (sscanf (s, "A%d ", &i) == 1) /* Ai <coeff of degree i> */
        {
          if (i > MAX_DEGREE)
            {
              fprintf (stderr, "Error, too large degree %d\n", i);
              exit (1);
            }
          if (mpz_set_str (X[i], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading X[%d]", i);
              exit (1);
            }
        }
      else if (strncmp (s, "R1 ", 3) == 0)
        {
          if (mpz_set_str (Y[1], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading Y1");
              exit (1);
            }
        }
      else if (strncmp (s, "R0 ", 3) == 0)
        {
          if (mpz_set_str (Y[0], s + 3, 0) != 0)
            {
              fprintf (stderr, "Error while reading Y0");
              exit (1);
            }
        }
    }
}

static void
read_cwi (mpz_t N, mpz_t *X, mpz_t *Y, mpz_t M)
{
  int npoly; /* number of polynomials */
  int degree, i;

  if (mpz_inp_str (N, stdin, 0) == 0)
    {
      fprintf (stderr, "Error, can't read number to factor N\n");
      exit (1);
    }
  if (mpz_inp_str (M, stdin, 0) == 0)
    {
      fprintf (stderr, "Error, can't read common root M\n");
      exit (1);
    }
  if (scanf ("%d\n", &npoly) != 1)
    {
      fprintf (stderr, "Error, can't read number of polynomials\n");
      exit (1);
    }
  if (npoly != 2)
    {
      fprintf (stderr, "Error, case of >=3 polynomials not yet implemented\n");
      exit (1);
    }
  if (scanf ("%d\n", &degree) != 1)
    {
      fprintf (stderr, "Error, can't read degree of 1st polynomial\n");
      exit (1);
    }
  if (degree != 1)
    {
      fprintf (stderr, "Error, first polynomial must be linear\n");
      exit (1);
    }
  if (mpz_inp_str (Y[0], stdin, 0) == 0)
    {
      fprintf (stderr, "Error, can't read coefficient Y0\n");
      exit (1);
    }
  if (mpz_inp_str (Y[1], stdin, 0) == 0)
    {
      fprintf (stderr, "Error, can't read coefficient Y1\n");
      exit (1);
    }
  if (scanf ("%d\n", &degree) != 1)
    {
      fprintf (stderr, "Error, can't read degree of 2nd polynomial\n");
      exit (1);
    }
  if (degree > MAX_DEGREE)
    {
      fprintf (stderr, "Error, too large degree\n");
      exit (1);
    }
  for (i = 0; i <= degree; i++)
    {
      if (mpz_inp_str (X[i], stdin, 0) == 0)
        {
          fprintf (stderr, "Error, can't read coefficient X[%d]\n", i);
          exit (1);
        }
    }
}

/************************** output routines **********************************/

static void
out_cwi (mpz_t N, mpz_t *X, int deg, mpz_t *Y, int degY, mpz_t M)
{
  int i;

  if (degY != 1)
    {
      fprintf (stderr, "Error, only degree 1 implemented for CWI format output\n");
      exit (1);
    }

  mpz_out_str (stdout, 10, N);
  printf ("\n");
  mpz_out_str (stdout, 10, M);
  printf ("\n\n2\n\n1\n");
  mpz_out_str (stdout, 10, Y[0]);
  printf (" ");
  mpz_out_str (stdout, 10, Y[1]);
  printf ("\n\n%d\n", deg);
  for (i = 0; i <= deg; i++)
    {
      mpz_out_str (stdout, 10, X[i]);
      printf (" ");
    }
  printf ("\n\n10000000\n");
}

static void
out_franke (mpz_t N, mpz_t *X, int deg, mpz_t *Y, int degY, mpz_t M)
{
  int i;

  if (degY != 1)
    {
      fprintf (stderr, "Error, only degree 1 implemented for FK format output\n");
      exit (1);
    }

  mpz_out_str (stdout, 10, N);  printf ("\n");
  for (i = deg; i >= 0; i--)
    {
      printf ("X%d ", i);
      mpz_out_str (stdout, 10, X[i]);
      printf ("\n");
    }
  printf ("Y1 "); mpz_out_str (stdout, 10, Y[1]); printf ("\n");
  printf ("Y0 "); mpz_out_str (stdout, 10, Y[0]); printf ("\n");
  if (mpz_cmp_ui (M, 0) != 0)
    gmp_printf ("M %Zd\n", M);
  printf ("# parameters for the linear and algebraic polynomials:\n");
  printf ("0 2000000 2.7 27 54\n");
  printf ("0 4000000 2.7 27 54\n");
}

static void
out_cado (mpz_t N, mpz_t *X, int degX, mpz_t *Y, int degY, mpz_t M)
{
  int i;

  printf ("n: ");
  mpz_out_str (stdout, 10, N);
  printf ("\n");
  for (i = degX; i >= 0; i--)
    {
      printf ("c%d: ", i);
      mpz_out_str (stdout, 10, X[i]);
      printf ("\n");
    }
  for (i = degY; i >= 0; i--)
    {
      printf ("Y%d: ", i);
      mpz_out_str (stdout, 10, Y[i]);
      printf ("\n");
    }
  printf ("m: ");
  mpz_out_str (stdout, 10, M);
  printf ("\n");
}

static void
out_msieve (mpz_t N, mpz_t *X, int deg, mpz_t *Y, int degY)
{
  int i;

  if (degY != 1)
    {
      fprintf (stderr, "Error, only degree 1 implemented for msieve format output\n");
      exit (1);
    }

  gmp_printf ("N %Zd\n", N);
  gmp_printf ("R0 %Zd\n", Y[0]);
  gmp_printf ("R1 %Zd\n", Y[1]);
  for (i = 0; i <= deg; i++)
    gmp_printf ("A%d %Zd\n", i, X[i]);
}

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s [options] < file1 > file2\n", argv0);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "        -if xxx - specify input format\n");
  fprintf (stderr, "        -of yyy - specify output format\n");
  fprintf (stderr, "Available formats:\n");
  fprintf (stderr, "        cado   - CADO-NFS format (default)\n");
  fprintf (stderr, "        fk     - Franke-Kleinjung format\n");
  fprintf (stderr, "        ggnfs  - GGNFS format\n");
  fprintf (stderr, "        cwi    - CWI format\n");
  fprintf (stderr, "        msieve - msieve format\n");
}

#define FORMAT_CADO   0
#define FORMAT_FK     1
#define FORMAT_GGNFS  2
#define FORMAT_CWI    3
#define FORMAT_MSIEVE 4

int
main (int argc, char *argv[])
{
  mpz_t N, *X, *Y, M;
  int iformat = FORMAT_CADO;
  int oformat = FORMAT_CADO;
  int i, degX = 0, degY = 0;
  
  while (argc > 1)
    {
      if (argc >= 3 && strcmp (argv[1], "-if") == 0)
        {
          if (strcmp (argv[2], "cado") == 0)
            iformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            iformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            iformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            iformat = FORMAT_CWI;
          else if (strcmp (argv[2], "msieve") == 0)
            iformat = FORMAT_MSIEVE;
          else
            {
              fprintf (stderr, "Error, invalid format: %s\n", argv[2]);
              exit (1);
            }
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-of") == 0)
        {
          if (strcmp (argv[2], "cado") == 0)
            oformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            oformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            oformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            oformat = FORMAT_CWI;
          else if (strcmp (argv[2], "msieve") == 0)
            oformat = FORMAT_MSIEVE;
          else
            {
              fprintf (stderr, "Error, invalid format: %s\n", argv[2]);
              exit (1);
            }
          argc -= 2;
          argv += 2;
        }
      else
        {
          usage (argv[0]);
          exit (1);
        }
    }

  mpz_init (N);
  X = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
  Y = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
  for (i = 0; i <= MAX_DEGREE; i++)
    {
      mpz_init (X[i]);
      mpz_init (Y[i]);
    }
  mpz_init (M);

  if (iformat == FORMAT_FK)
    read_franke (N, X, Y, M);
  else if (iformat == FORMAT_GGNFS || iformat == FORMAT_CADO)
    read_ggnfs (N, X, Y, M);
  else if (iformat == FORMAT_CWI)
    read_cwi (N, X, Y, M);
  else if (iformat == FORMAT_MSIEVE)
    read_msieve (N, X, Y);
  else
    {
      fprintf (stderr, "Input format not yet implemented\n");
      exit (1);
    }

  for (degX = MAX_DEGREE; degX > 0 && mpz_cmp_ui (X[degX], 0) == 0; degX --);
  for (degY = MAX_DEGREE; degY > 0 && mpz_cmp_ui (Y[degY], 0) == 0; degY --);

  if (mpz_cmp_ui (M, 0) == 0 && degY == 1)
    {
      mpz_t t;
      /* M = -Y0/Y1 mod N */
      mpz_invert (M, Y[1], N);
      mpz_neg (M, M);
      mpz_mul (M, M, Y[0]);
      mpz_mod (M, M, N);
      /* check M is also a root of the algebraic polynomial mod N */
      mpz_init_set (t, X[degX]);
      for (i = degX - 1; i >= 0; i --)
        {
          mpz_mul (t, t, M);
          mpz_add (t, t, X[i]);
          mpz_mod (t, t, N);
        }
      if (mpz_cmp_ui (t, 0) != 0)
        {
          fprintf (stderr, "Polynomials have no common root mod N\n");
          exit (1);
        }
      mpz_clear (t);
    }

  if (oformat == FORMAT_CWI)
    out_cwi (N, X, degX, Y, degY, M);
  else if (oformat == FORMAT_FK)
    out_franke (N, X, degX, Y, degY, M);
  else if (oformat == FORMAT_CADO)
    out_cado (N, X, degX, Y, degY, M);
  else if (oformat == FORMAT_MSIEVE)
    out_msieve (N, X, degX, Y, degY);
  else
    {
      fprintf (stderr, "Output format not yet implemented\n");
      exit (1);
    }

  mpz_clear (M);
  for (i = 0; i <= MAX_DEGREE; i++)
    {
      mpz_clear (X[i]);
      mpz_clear (Y[i]);
    }
  mpz_clear (N);
  free (X);
  free (Y);

  return 0;
}
