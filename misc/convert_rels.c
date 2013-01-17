/* Convert relations from format X to format Y,
   where X, Y are in {CADO, Franke/Kleinjung, GGNFS}.
   
   Run with no argument to get usage.

   TODO: msieve uses the same format as CADO-NFS, except for free relations,
   which should be converted to single lines like:

   p,0:

*/

#define FORMAT_CADO 0  /* CADO format (default output format) */
#define FORMAT_FK   1  /* Franke-Kleinjung's format */
#define FORMAT_GGNFS 2 /* GGNFS format (default input format) */
#define FORMAT_CWI   3 /* CWI format */

static const char *format_names[4] = {
  "CADO", "Franke-Kleinjung", "GGNFS", "CWI"
};

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> /* for UINT32_MAX */
#include <string.h>
#include <assert.h>
#include "../utils/gzip.h"
#include "gmp.h"

#define MAX_PRIMES 255 /* maximal number of factor base primes */
#define MAX_LPRIMES 3  /* maximal number of large primes */
#define DEGF_MAX 6
#define REPS 10

int line_number;

uint32_t
get_uint32 (FILE *fp)
{
  uint32_t w;
  fread (&w, sizeof (uint32_t), 1, fp);
  return w;
}

int32_t
get_int32 (FILE *fp)
{
  int32_t w;
  fread (&w, sizeof (int32_t), 1, fp);
  return w;
}

typedef struct
{
  int64_t a;
  uint32_t b;
  unsigned int rfb_entries;
  unsigned int afb_entries;
  unsigned int sp_entries;
  unsigned int num_lrp;
  unsigned int num_lap;
  unsigned long rprimes[MAX_PRIMES]; /* index i of prime p(i) */
  unsigned long rexp[MAX_PRIMES];
  unsigned long aprimes[MAX_PRIMES];
  unsigned long aexp[MAX_PRIMES];
  unsigned long sprimes[MAX_PRIMES];
  unsigned long sexp[MAX_PRIMES];
  unsigned long large_rprimes[MAX_LPRIMES];
  unsigned long large_aprimes[MAX_LPRIMES];
  int32_t qcb1, qcb2; /* quadratic characters */
  int end_of_set;
} relation_t;

/* prints a large prime, given its low and high 32-bit parts.
   Currently not used. */
void
print_large_prime (uint32_t l, uint32_t h)
{
  mpz_t P;

  mpz_init (P);
  mpz_set_ui (P, h);
  mpz_mul_2exp (P, P, 32);
  mpz_add_ui (P, P, l);
  mpz_out_str (stdout, 16, P);
  mpz_clear (P);
}

/* Return 1 if the largest prime in rel is < 2^lpb or if lpb == 0, 
   return 0 otherwise */
int
checksize_relation_cado (relation_t *rel, int32_t *rfb, int32_t *afb, int lpb)
{
  unsigned int i;
  unsigned long p;

  if (lpb == 0)
    return 1;

  /* rational side */
  for (i = 0; i < rel->rfb_entries; i++)
    {
      assert(rel->rexp[i] > 0);
      if (rfb == NULL)
        p = rel->rprimes[i];
      else
        p = rfb[rel->rprimes[i]];
      if (p >> lpb)
        return 0;
    }
  for (i = 0; i < rel->num_lrp; i++)
    {
      p = rel->large_rprimes[i];
      if (p >> lpb)
        return 0;
    }

  /* algebraic side */
  for (i = 0; i < rel->afb_entries; i++)
    {
      assert (rel->aexp[i] > 0);
      if (afb == NULL)
        p = rel->aprimes[i];
      else
        p = afb[rel->aprimes[i]];
      if (p >> lpb)
        return 0;
    }
  for (i = 0; i < rel->sp_entries; i++)
    {
      if (rel->sexp[i] > 0)
        {
          p = rel->sprimes[i];
          if (p >> lpb)
            return 0;
        }
    }
  for (i = 0; i < rel->num_lap; i++)
    {
      p = rel->large_aprimes[i];
      if (p >> lpb)
        return 0;
    }

  return 1;
}

/* output a relation in CADO format:
   a,b:p1,p2,...,pn:q1,q2,...,qn
   where a and p are printed in decimal, p1,...,qn are printed in hexadecimal,
   p1,p2,...,pn are the primes on the rational side.
*/
void
print_relation_cado (relation_t *rel, int32_t *rfb, int32_t *afb)
{
  unsigned int i, j;
  int need_comma;

  printf ("%ld,%u:", rel->a, rel->b);
  need_comma = 0;

  /* rational side */
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      {
	if (need_comma++) putchar (',');
        if (rfb == NULL) /* primes are given directly */
          printf ("%lx", rel->rprimes[i]);
        else /* primes are given by their index */
          printf ("%x", rfb[rel->rprimes[i]]);
        fflush (stdout);
      }
  for (i = 0; i < rel->num_lrp; i++)
    {
      if (need_comma++) putchar (',');
      printf ("%lx", rel->large_rprimes[i]);
    }

  putchar (':'); 
  need_comma = 0;

  /* algebraic side */
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      {
	if (need_comma++) putchar (',');
        if (afb == NULL) /* primes are given directly */
          printf ("%lx", rel->aprimes[i]);
        else /* primes are given by their index */
          printf ("%x", afb[rel->aprimes[i]]);
      }
  for (i = 0; i < rel->sp_entries; i++)
    {
      if (rel->sexp[i] > 0)
	for (j = 0; j < rel->sexp[i]; j++)
	  {
	    if (need_comma++) putchar (',');
	    printf ("%lx", rel->sprimes[i]);
	  }
    }
  for (i = 0; i < rel->num_lap; i++)
    {
      if (need_comma++) putchar (',');
      printf ("%lx", rel->large_aprimes[i]);
    }
  printf ("\n");
}

/* output a relation in Franke/Kleinjung format:
   - the first line is "F 0 X 5 1" where the two last integers are the degrees
   - then it comes by groups of 3 consecutive lines (W gives a and b in
     hexadecimal form, X gives the algebraic factors, Y the linear factors.
   W 2335e65b 1d9
   X 96f8323 1413403 1369 16B5 1A0C77 347 55D 3 5 65 71 9D F42439
   Y 444B A6CF C827 2C8627 B 25 209 527 B 2 2 2
*/
void
print_relation_fk (relation_t *rel, int32_t *rfb, int32_t *afb, int iformat)
{
  unsigned int i, j;

  if (rel->a > 0)
    printf ("W %lx %x\n", rel->a, rel->b);
  else
    printf ("W -%lx %x\n", -rel->a, rel->b);
  /* algebraic side */
  printf ("X");
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      printf (" %x", (iformat == FORMAT_GGNFS) ? afb[rel->aprimes[i]]
              : (int32_t) rel->aprimes[i]);
  for (i = 0; i < rel->sp_entries; i++)
    for (j = 0; j < rel->sexp[i]; j++)
      printf (" %lx", rel->sprimes[i]);
  for (i = 0; i < rel->num_lap; i++)
    printf (" %lx", rel->large_aprimes[i]);
  printf ("\n");
  /* rational side */
  printf ("Y");
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      printf (" %x", (iformat == FORMAT_GGNFS) ? rfb[rel->rprimes[i]]
              : (int32_t) rel->rprimes[i]);
  for (i = 0; i < rel->num_lrp; i++)
    printf (" %lx", rel->large_rprimes[i]);
  printf ("\n");
}


/* Map 0, 1, ..., 61 to '0',...,'9','A',...,'Z','a',...'z' */
static unsigned char
cwi_idx_to_char (unsigned char i)
{
  if (i < 10)
    return i + '0';
  if (i < 36)
    return i - 10 + 'A';
  if (i < 62)
    return i - 36 + 'a';
  abort();
}


/* Map '0',...,'9','A',...,'Z','a',...'z' to 0, 1, ..., 61 */
static unsigned char
cwi_char_to_idx (unsigned char i)
{
  assert (i >= '0');
  if (i <= '9')
    return i - '0';
  assert (i >= 'A');
  if (i <= 'Z')
    return i - 'A' + 10;
  assert (i >= 'a');
  if (i <= 'z')
    return i - 'a' + 36;
  abort();
}


static void 
print_primes_cwi (unsigned long *primes, unsigned long *exp,
                  int32_t *fb, const unsigned int np)
{
  unsigned long i,j;

  for (i = 0; i < np; i++)
    for (j = 0; j < (exp != NULL ? exp[i] : 1); j++)
      {
        if (fb == NULL) /* primes are given directly */
          printf (" %lu", primes[i]);
        else /* primes are given by their index */
          printf (" %u", fb[primes[i]]);
      }
}

void
print_relation_cwi (relation_t *rel, int32_t *rfb, int32_t *afb, 
                    const int alg_first)
{
  unsigned long i, nr, na;
  
  nr = 0;
  na = 0;
  for (i = 0; i < rel->rfb_entries; i++)
    nr += rel->rexp[i];
  for (i = 0; i < rel->afb_entries; i++)
    na += rel->aexp[i];
  for (i = 0; i < rel->sp_entries; i++)
    na += rel->sexp[i];
  
  printf ("01%c%c %ld %u", cwi_idx_to_char (alg_first ? na : nr), 
          cwi_idx_to_char (alg_first ? nr : na), rel->a, rel->b);

  for (i = 0; i < 2; i++)
    {
      /* if alg_first == 0, do rational primes for side == 0 */
      if (i == (alg_first ? 1 : 0))
        {
          /* rational side */
          print_primes_cwi (rel->rprimes, rel->rexp, rfb, rel->rfb_entries);
          print_primes_cwi (rel->large_rprimes, NULL, NULL, rel->num_lrp);
        }
      else
        {
          /* algebraic side and special primes */
          print_primes_cwi (rel->aprimes, rel->aexp, afb, rel->afb_entries);
          print_primes_cwi (rel->sprimes, rel->sexp, NULL, rel->sp_entries);
          print_primes_cwi (rel->large_aprimes, NULL, NULL, rel->num_lap);
        }
    }
  printf ("%c\n", rel->end_of_set ? ';' : ':');

  fflush (stdout);
}

/* norm <- |f(a,b)| = |f[d]*a^d + ... + f[0]*b^d| */
void
eval_algebraic (mpz_t norm, mpz_t *f, int d, long a, unsigned long b)
{
  mpz_t B; /* powers of b */
  mpz_init_set_ui (B, 1);
  mpz_set (norm, f[d]);
  while (d-- > 0)
    {
      mpz_mul_si (norm, norm, a);
      mpz_mul_ui (B, B, b);
      mpz_addmul (norm, f[d], B);
    }
  mpz_clear (B);
  mpz_abs (norm, norm);
}


/* Add the prime p to *primes with exponents in *exp. 
   If p is already present in *primes, its exponent is increased. 
   i is the number of entries in *primes so far, the new number 
   of entries (after adding p) is the return value. */
static void 
add_prime (unsigned long *primes, unsigned long *exps, unsigned int *n, 
           const unsigned long p)
{
  unsigned int i;
  for (i = 0; i < *n; i++)
    if (primes[i] == p)
      {
        exps[i] ++;
        return;
      }
  
  assert (*n < MAX_PRIMES);
  primes[*n] = p;
  exps[*n] = 1;
  (*n)++;
}

static void
trialdiv (unsigned long *primes, unsigned long *exps, unsigned int *n, 
          mpz_t N)
{
  unsigned long k = 2; /* prime factors to try */
  
  while (mpz_cmp_ui (N, 1) != 0)
    {
      if (mpz_probab_prime_p (N, REPS))
        {
          if (!mpz_fits_ulong_p (N))
            {
              gmp_fprintf (stderr, "Error, cofactor %Zd does not fit in an unsigned long\n", N);
              abort ();
            }
          add_prime (primes, exps, n, mpz_get_ui (N));
          mpz_set_ui (N, 1);
        }
      else /* search for a prime factor of N */
        {
          while (mpz_divisible_ui_p (N, k) == 0)
            k = (k + 1) | 1;
          do
            {
              add_prime (primes, exps, n, k);
              mpz_divexact_ui (N, N, k);
            }
          while (mpz_divisible_ui_p (N, k));
        }
    }
}

/* Read one relation in CADO format from fp, and put it in rel.
   Return 1 if relation is valid, 0 if invalid, -1 if end of file. */
int
read_relation_cado (FILE *fp, relation_t *rel)
{
  int c;
  unsigned long p;

  if (fscanf (fp, "%ld,%u:", &(rel->a), &(rel->b)) != 2)
    {
      if (feof (fp))
        return 0;
      if ((c = getc (fp)) == '#')
        {
          /* Discard comments? Copy to output? For now, discard */
          do
              c = getc (fp);
          while (c != '\n');
          return 0;
        }
      ungetc (c, fp);
      
      fprintf (stderr, "Error, invalid relation (could not read a,b): ");
      do
        {
          c = getc (fp);
          fputc (c, stderr);
        }
      while (c != '\n');
      exit (1);
    }

  rel->rfb_entries = 0; /* number of rational primes */
  do
    if (fscanf (fp, "%lx", &p) == 1)  /* new rational prime */
      add_prime (rel->rprimes, rel->rexp, &rel->rfb_entries, p);
  while ((c = getc (fp)) == ',');
  rel->num_lrp = 0;

  if (c != ':')
    {
      fprintf (stderr, "Error, invalid relation (expected \":\"): ");
      exit(1);
    }

  rel->afb_entries = 0; /* number of algebraic primes */
  /* fscanf() skips leading whitespace so we need to test for a 
     newline here which can occur if there are no algebraic primes */
  if ((c = getc (fp)) != '\n')
    {
      ungetc (c, fp);
      do
        if (fscanf (fp, "%lx", &p) == 1) /* new algebraic prime */
          add_prime (rel->aprimes, rel->aexp, &rel->afb_entries, p);
      while ((c = getc (fp)) == ',');
    }
    
  if (c != '\n')
    {
      fprintf (stderr, "Error, invalid relation (expected newline): ");
      exit(1);
    }

  rel->sp_entries = 0;
  rel->num_lap = 0;
  rel->end_of_set = 1;

  return 1;
}

int 
fk_read_line (char *lp, const int maxlen, FILE *fp, char *file)
{
  int i, skip;

  do {
    if (fgets (lp, maxlen, fp) == NULL)
      return 0; /* Possibly EOF */

    line_number ++;
    
    i = strlen (lp); /* fgets() always puts a '\0' so this is safe */
    if (i == maxlen - 1 && lp[i] != '\n')
      {
	fprintf (stderr, "Error, input line too long\n");
	exit(1);
      }
    skip = 0;
    if (i == 0 || lp[i - 1] != '\n')
      {
	fprintf (stderr, "Skipping incomplete line %d in file %s:\n",
                 line_number, file);
        fprintf (stderr, "   %s\n", lp);
	skip = 1;
      }
  } while (skip || lp[0] == '#' || strncmp (lp, "F 0 X", 5) == 0);
  /* Lines beginning with '#' are comments and are skipped over. Each file
     begins with "F 0 X", but the input may be several files concatenated,
     so we allow (and ignore) "F 0 X" anywhere */
  return 1;
}

int 
fk_read_primes (char **lp, unsigned long *exponent, unsigned long *primes)
{
  unsigned int i;

  i = 0; /* number of primes */
  while ((*lp)[0] != '\n')
    {
      char *nlp;
      unsigned long p;
      
      p = strtoul (*lp, &nlp, 16);
      if (nlp == *lp)
	{
	  fprintf (stderr, "Error, could not parse prime\n");
	  exit (1);
	}
      *lp = nlp;
      
      add_prime (primes, exponent, &i, p);
    }
  
  return i;
}

int
read_relation_fk (FILE *fp, relation_t *rel, char *file)
{
  char line[512];
  char *lp;

  /* Read the "W" line with the a and b values */
  if (fk_read_line (line, 512, fp, file) == 0)
    return -1; /* Signal EOF */

  if (line[0] != 'W' || line[1] != ' ')
    {
      fprintf (stderr, "Error, no W line at start of relation\n");
      exit (1);
    }
  lp = line + 2;

  /* Read a and b values */
  rel->a = strtol (lp, &lp, 16);
  rel->b = strtoul (lp, &lp, 16);
  if (lp[0] != '\n')
    {
      fprintf (stderr, "Error, could not read a and/or b value\n");
      exit (1);
    }
  
  /* Read the "X" line, which has the algebraic primes */
  if (fk_read_line (line, 512, fp, file) != 1)
    {
      fprintf (stderr, "Error, incomplete relation at end of file\n");
      exit (1);
    }
  if (line[0] != 'X' || (line[1] != ' ' && line[1] != '\n'))
    {
      fprintf (stderr, "Error, no X line after W line\n");
      exit (1);
    }
  lp = line + (line[1] == ' ' ? 2 : 1);

  rel->afb_entries = fk_read_primes (&lp, rel->aexp, rel->aprimes);
  rel->num_lap = 0;

  /* Read the "Y" line, which has the rational primes */
  if (fk_read_line (line, 512, fp, file) != 1)
    {
      fprintf (stderr, "Error, incomplete relations at end of file\n");
      exit (1);
    }
  if (line[0] != 'Y' || (line[1] != ' ' && line[1] != '\n'))
    {
      fprintf (stderr, "Error, no Y line after X line\n");
      exit (1);
    }
  lp = line + (line[1] == ' ' ? 2 : 1);

  rel->rfb_entries = fk_read_primes (&lp, rel->rexp, rel->rprimes);
  rel->num_lrp = 0;
  rel->sp_entries = 0;
  rel->end_of_set = 1;

  return 1;
}

/* Read one relation in CWI format from fp, and put it in rel.
   Return 1 if relation is valid, 0 if invalid, -1 if end of file. */
int
read_relation_cwi (FILE *fp, relation_t *rel, const int alg_first)
{
  int ret, side;
  unsigned long p;
  unsigned int i; /* number of given primes */
  char c, flag[4];

  ret = fscanf (fp, "%4s %ld %u", flag, &(rel->a), &(rel->b));
  if (ret != 3)
    {
      if (feof (fp))
        return 0;
      fprintf (stderr, "Error, invalid relation: ");
      do
        {
          c = getc (fp);
          fputc (c, stderr);
        }
      while (c != '\n');
      exit (1);
    }

  /* check flag=01xy */
  if (flag[0] != '0' || flag[1] != '1')
    {
      fprintf (stderr, "Error, flag differs from 01xy: %.4s\n", flag);
      exit (1);
    }

  for (side = 0; side < 2; side++)
    {
      /* if alg_first == 0, do rational primes for side == 0 */
      int rat_now = side == (alg_first ? 1 : 0);
      unsigned int *fb_entries = 
        rat_now ? &rel->rfb_entries : &rel->afb_entries;
      unsigned long *primes = rat_now ? rel->rprimes : rel->aprimes;
      unsigned long *exps = rat_now ? rel->rexp : rel->aexp;
      char *side_name = rat_now ? "rational" : "algebraic";
      unsigned int n;

      n = cwi_char_to_idx(flag[2 + side]);
      *fb_entries = 0;
      for (i = 0; i < n; i++)
        {
          if (fscanf (fp, " %ld", &p) != 1)
            {
              fprintf (stderr, "Error, can't read next %s prime\n", side_name);
              exit (1);
            }
          add_prime (primes, exps, fb_entries, p);
        }
    }

  rel->num_lrp = 0;
  rel->sp_entries = 0;
  rel->num_lap = 0;

  ret = fscanf (fp, "%c\n", &c);
  if (ret != 1 || (c != ';' && c != ':'))
    {
      fprintf (stderr, "Error, invalid relation for a=%ld b=%u\n",
               rel->a, rel->b);
      exit (1);
    }

  rel->end_of_set = (c == ';');

  return 1;
}

/* Read one relation in GGNFS format (adapted from function dataConvertToRel
   in GGNFS, file rels.c).
   Return value: 1 if relation is ok, 0 if invalid, -1 if end of file.
*/
int
read_relation_ggnfs (FILE *fp, relation_t *rel, mpz_t norm, mpz_t *f,
                         int degf, int32_t *afb)
{
  uint32_t size_field; /* 32 bits */
  unsigned int i = 0; /* gcc 4.1.2 warns that i may be used uninitialized,
                         but this does not seem correct */
  unsigned int k;
  unsigned sp_entries; /* Number of sprimes according to file */
  int32_t p;

  size_field = get_uint32 (fp);
  rel->a = (int32_t) get_int32 (fp);
  rel->b = (int32_t) get_uint32 (fp);
  eval_algebraic (norm, f, degf, rel->a, rel->b);
  rel->rfb_entries = size_field >> 24;
  for (i = 0; i < rel->rfb_entries; i++)
    {
      rel->rprimes[i] = get_uint32 (fp);
      rel->rexp[i] = get_uint32 (fp);
    }
  rel->afb_entries = (size_field >> 16) & 255;
  for (i = 0; i < rel->afb_entries; i++)
    {
      rel->aprimes[i] = get_uint32 (fp);
      p = afb[rel->aprimes[i]];
      rel->aexp[i] = get_uint32 (fp);
      for (k = 0; k < rel->aexp[i]; k++)
        {
          if (!mpz_divisible_ui_p (norm, p))
            {
              static unsigned int count = 0;
              if (count++ < 10)
                fprintf (stderr, "Warning, f(%ld,%u) not divisible by factor base prime %d\n", rel->a, rel->b, p);
            }
          mpz_divexact_ui (norm, norm, p);
        }
    }
  sp_entries = (size_field >> 8) & 255;
  for (i = 0; i < sp_entries; i++)
    {
      /* since GGNFS format gives the *index* of the special primes
         and not their value, we have no easy mean to recover them.
         The workaround is to factor the norm residue after all factor
         base and large primes have been divided out. We thus put the
         rel->sp_entries to 0 for now. */
      get_uint32 (fp); /* reads the special prime index (not used) */
      get_int32 (fp); /* reads the exponent field (not used) */
    }
  rel->qcb1 = get_uint32 (fp);
  rel->qcb2 = get_uint32 (fp);
  rel->num_lrp = (size_field >> 6) & 3;
  for (i = 0; i < rel->num_lrp; i++)
    rel->large_rprimes[i] = (unsigned long) get_uint32 (fp);
  rel->num_lap = (size_field >> 4) & 3;
  for (i = 0; i < rel->num_lap; i++)
    {
      unsigned long p;
      p = (unsigned long) get_uint32 (fp);
      get_uint32 (fp); /* corresponding root, not used here */
      rel->large_aprimes[i] = p;
      if (!mpz_divisible_ui_p (norm, p))
        fprintf (stderr, "Warning, f(%ld,%u) not divisible by large prime %lu\n",
                 rel->a, rel->b, p);
      else
        mpz_divexact_ui (norm, norm, p);
    }

  rel->sp_entries = 0;
  trialdiv (rel->sprimes, rel->sexp, &rel->sp_entries, norm);
  if (rel->sp_entries != sp_entries)
    {
      /* Trial division found a different number of primes than the 
         number of sprimes according to the file. Discard this relation. */
      fprintf (stderr, "Error: %u special primes according to file but "
               "trial division found %u\n", sp_entries, rel->sp_entries);
      return 0;
    }

  rel->end_of_set = 1;
  return 1; /* valid relation */
}

/* Read relations from file fp, with input format is 'iformat'. 
   and prints them with given output format 'oformat'.

   degf is the degree of the algebraic polynomial (needed for GGNFS
   input formats).
   
   For GGNFS input: the rational factor base is into rfb, the algebraic factor
   base is into afb, and the coefficients of the algebraic polynomial are into
   f[0]...f[degf].
*/
unsigned long
convert_relations (char *rels, int32_t *rfb, int32_t *afb, mpz_t *f, int degf,
                   int iformat, int oformat, int lpb, int verbose, int multi, 
                   const int in_alg_first, const int out_alg_first)
{
  FILE *fp;
  relation_t rel[1];
  uint32_t relsInFile;
  unsigned long outputRels = 0;
  unsigned int j;
  mpz_t norm;
  int ok;
  int compressed;

  fp = fopen_compressed_r (rels, &compressed, NULL);
  if (fp == NULL)
    {
      fprintf (stderr, "Error, unable to open relation file %s\n", rels);
      exit (1);
    }

  mpz_init (norm);

  /* the first 4 bytes give the number of relations in the file */
  rewind (fp);
  line_number = 0;

  if (iformat == FORMAT_GGNFS)
    {
      fread (&relsInFile, sizeof (int32_t), 1, fp);
      if (verbose)
        fprintf (stderr, "File %s: header says %d relations.\n",
                 rels, relsInFile);
    }
  else
    relsInFile = UINT32_MAX;

  for (j = 0; feof (fp) == 0 && j < relsInFile; j++)
    {
      if (iformat == FORMAT_GGNFS)
        ok = read_relation_ggnfs (fp, rel, norm, f, degf, afb);
      else if (iformat == FORMAT_CADO)
        ok = read_relation_cado (fp, rel);
      else if (iformat == FORMAT_FK)
        ok = read_relation_fk (fp, rel, rels);
      else if (iformat == FORMAT_CWI)
        ok = read_relation_cwi (fp, rel, in_alg_first);
      else
        {
          fprintf (stderr, "Error, unknown format %d\n", iformat);
          exit (1);
        }

      if (ok == -1) /* end of file */
        break;

      if (ok != 0 && checksize_relation_cado (rel, rfb, afb, lpb))
        {
          if (!multi)
            {
              unsigned long i;
              for (i = 0; i < rel->rfb_entries; i++)
                rel->rexp[i] = 1;
              for (i = 0; i < rel->afb_entries; i++)
                rel->aexp[i] = 1;
              for (i = 0; i < rel->sp_entries; i++)
                rel->sexp[i] = 1;
            }
            
          switch (oformat)
            {
            case FORMAT_CADO:
              print_relation_cado (rel, rfb, afb);
              break;
            case FORMAT_FK:
              print_relation_fk (rel, rfb, afb, iformat);
              break;
            case FORMAT_CWI:
              print_relation_cwi (rel, rfb, afb, out_alg_first);
              break;
            default:
              fprintf (stderr, "Error, unknown format %d\n", oformat);
              exit (1);
            }
          fflush (stdout);
          outputRels ++;
        }
    }

  if (iformat == FORMAT_GGNFS)
    {
      if (j != relsInFile)
        {
          fprintf (stderr, "Warning: number of relations read differs from expectation\n");
          exit (1);
        }
    }
  else
    relsInFile = j;

  mpz_clear (norm);

  if (compressed)
    pclose (fp);
  else
    fclose (fp);

  return outputRels;
}

#define BF_DELIMITER 0x0a
#define END_POLY "END_POLY"
#define END_HEADER "END_HEADER"

/* read fp until the string s is found,
   the next char read by getc (fp) is the one following s. */
void
get_string (char *s, FILE *fp, int verbose)
{
  int l, i, c;

  l = strlen (s);
  i = 0;
  while (i < l)
    {
      c = getc (fp);
      if (verbose)
	putc (c, stderr);
      if (c == s[i])
	i ++;
      else
	i = 0;
    }
}

/* read factor base from file fp, puts rational factor base into *rfb,
   (size into *rfb_size), algebraic factor base into *afb (size into
   *afb_size), and put coefficients from algebraic polynomial in f[0]..f[d].
   Returns degree d.
*/
int
read_fb (FILE *fp, int32_t **rfb, int32_t *rfb_size, int32_t **afb, 
         int32_t *afb_size, int verbose, mpz_t *f)
{
  int c, i, degf;
  long RFBsize, AFBsize;

  /* read polynomial */
  get_string ("Y1: ", fp, verbose);
  mpz_inp_str (f[0], fp, 10);
  if (verbose)
    mpz_out_str (stderr, 10, f[0]);
  c = getc (fp); /* newline */
  fscanf (fp, "c%d: ", &degf);
  if (degf > DEGF_MAX)
    {
      fprintf (stderr, "Error, too large degree\n");
      exit (1);
    }
  if (verbose)
    fprintf (stderr, "%cc%d: ", c, degf);
  mpz_inp_str (f[degf], fp, 10);
  if (verbose)
    mpz_out_str (stderr, 10, f[degf]);
  for (i = degf - 1; i >= 0; i--)
    {
      getc (fp); /* newline */
      fscanf (fp, "c%d: ", &c);
      if (c != i)
	{
	  fprintf (stderr, "Error, missing coefficient of degree %d\n", i);
	  exit (1);
	}
      if (verbose)
	fprintf (stderr, "\nc%d: ", i);
      mpz_inp_str (f[i], fp, 10);
      if (verbose)
	mpz_out_str (stderr, 10, f[i]);
    }
  get_string (END_POLY, fp, verbose);
  c = getc (fp);
  if (c != BF_DELIMITER)
    {
      fprintf (stderr, "Error, delimited expected\n");
      exit (1);
    }
  if (verbose)
    putc ('\n', stderr);

  /* read header */
  get_string ("RFBsize: ", fp, verbose);
  fscanf (fp, "%ld", &RFBsize);
  if (verbose)
    fprintf (stderr, "%ld", RFBsize);
  get_string ("AFBsize: ", fp, verbose);
  fscanf (fp, "%ld", &AFBsize);
  if (verbose)
    fprintf (stderr, "%ld", AFBsize);
  get_string (END_HEADER, fp, verbose);
  c = getc (fp);
  if (c != BF_DELIMITER)
    {
      fprintf (stderr, "Error, delimited expected\n");
      exit (1);
    }
  if (verbose)
    fprintf (stderr, "\nRFBsize=%ld AFBsize=%ld\n", RFBsize, AFBsize);
  *rfb_size = RFBsize;
  *rfb = (int32_t*) malloc (RFBsize * sizeof(int32_t));
  *afb_size = AFBsize;
  *afb = (int32_t*) malloc (AFBsize * sizeof(int32_t));
  for (i = 0; i < RFBsize; i++)
    {
      (*rfb)[i] = get_uint32 (fp);
      get_uint32 (fp); /* don't need the root mod g */
    }
  for (i = 0; i < AFBsize; i++)
    {
      (*afb)[i] = get_uint32 (fp);
      get_uint32 (fp); /* don't need the root mod f */
    }
  return degf;
}

static void
usage (char *s)
{
  fprintf (stderr, "Usage: %s [options] f0 f1 ...\n", s);
  fprintf (stderr, "       where f0 f1 ... are input files with relations\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -v       - verbose\n");
  fprintf (stderr, "       -if xxx  - specifies input format (default ggnfs)\n");
  fprintf (stderr, "       -of yyy  - specifies output format (default cado)\n");
  fprintf (stderr, "                  where xxx, yyy are in {cado, fk, cwi, cwi2, ggnfs}\n");
  fprintf (stderr, "                  cwi has rational side first, cwi2 algebraic first\n");
  fprintf (stderr, "       -fb file - factor base (needed for ggnfs input)\n");
  fprintf (stderr, "       -deg d   - algebraic degree (for fk output)\n");
  fprintf (stderr, "       -lpb l   - discard relations with primes >= 2^l\n");
  fprintf (stderr, "       -nomulti - print repeated primes only once\n");
}

int
main (int argc, char *argv[])
{
  FILE *fp;
  unsigned long num_rels = 0;
  int32_t *rfb = NULL, *afb = NULL;
  int32_t rfb_size, afb_size;
  int verbose = 0, i, degf = 0;
  int iformat = FORMAT_GGNFS; /* default input format */
  int oformat = FORMAT_CADO; /* default output format */
  mpz_t f[DEGF_MAX + 1];
  char *FB = NULL; /* name of factor base file */
  char *RELS; /* name of binary file to read */
  int num_files = 0;
  char *program_name = argv[0];
  int lpb = 0; /* 0 means no bound */
  int multi = 1, in_alg_first = 0, out_alg_first = 0;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-if") == 0)
	{
          if (strcmp (argv[2], "cado") == 0)
            iformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            iformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            iformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            iformat = FORMAT_CWI;
          else if (strcmp (argv[2], "cwi2") == 0)
            {
              iformat = FORMAT_CWI;
              in_alg_first = 1;
            }
          else
            {
              fprintf (stderr, "Unknown format: %s\n", argv[2]);
              exit (1);
            }
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-of") == 0)
	{
          if (strcmp (argv[2], "cado") == 0)
            oformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            oformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            oformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            oformat = FORMAT_CWI;
          else if (strcmp (argv[2], "cwi2") == 0)
            {
              oformat = FORMAT_CWI;
              out_alg_first = 1;
            }
          else
            {
              fprintf (stderr, "Unknown format: %s\n", argv[2]);
              exit (1);
            }
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-deg") == 0)
        {
          degf = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-v") == 0)
	{
	  verbose = 1;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-fb") == 0)
	{
	  FB = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-lpb") == 0)
	{
	  lpb = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (strcmp (argv[1], "-nomulti") == 0)
	{
	  multi = 0;
	  argv ++;
	  argc --;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (1);
	}
    }

  /* for ggnfs input, we need the factor base,
     and for Franke-Kleinjung output, we need the algebraic degree */
  if ((iformat == FORMAT_GGNFS && FB == NULL) || argc <= 1)
    {
      usage (program_name);
      exit (1);
    }

  fprintf (stderr, "Convert relations from %s to %s format\n",
           format_names[iformat], format_names[oformat]);
  if (lpb != 0)
    fprintf (stderr, "(keeping relations with primes < 2^%d only)\n", lpb);

  if (oformat == FORMAT_GGNFS)
    {
      fprintf (stderr, "Error, GGNFS output format not yet implemented\n");
      exit (1);
    }

  for (i = 0; i <= DEGF_MAX; i++)
    mpz_init (f[i]);

  /* check we can open the factor base for GGNFS input format */
  if (iformat == FORMAT_GGNFS)
    {
      fp = fopen (FB, "rb");
      if (fp == NULL)
        {
          fprintf (stderr, "Error, unable to open factor base file %s\n", FB);
          exit (1);
        }
      degf = read_fb (fp, &rfb, &rfb_size, &afb, &afb_size, verbose, f);
      fclose (fp);
    }
  
  while (argc > 1)
    {
      RELS = argv[1];
      argv ++;
      argc --;
      /* if output format is Franke-Kleinjung, print degree if provided */
      if (oformat == FORMAT_FK && num_files == 0 && degf != 0)
	printf ("F 0 X %d 1\n", degf);
      num_files ++;
      num_rels += convert_relations (RELS, rfb, afb, f, degf, iformat, oformat,
                                     lpb, verbose, multi, in_alg_first, 
                                     out_alg_first);
    }

  fprintf (stderr, "Output %lu relations from %d file(s).\n", num_rels,
	   num_files);

  free (rfb);
  free (afb);
  for (i = 0; i <= DEGF_MAX; i++)
    mpz_clear (f[i]);

  return 0;
}
