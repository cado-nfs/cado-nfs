/* Convert relations from format X to format Y,
   where X, Y are in {CADO, Franke/Kleinjung, GGNFS}.
   
   Run with no argument to get usage.

   TODO: msieve uses the same format as CADO-NFS, except for free relations,
   which should be converted to single lines like:

   p,0:


  
XXX: When converting from Franke-Kleinjung format to CADO format the following
command line is faster than the current code:  
  zgrep -vh -e "^F 0 X" -e "^#" $FILES | awk -f convert_rels.awk | tr " " "," | gzip -c --fast > rels.out.gz

  where:
    - convert_rels.awk contains
       BEGIN { RS = "" ; FS = "[WYX] |\n" }
       {
         printf("%s:%s:%s\n", $2 ,$6 ,$4);
       }
    - $FILES is a list of relations files (in .gz format in this example) where
      the three lines containing W Y and X should be separate by a blank line
*/

#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint> /* for UINT32_MAX */
#include <cinttypes>
#include <cstring>

#include <vector>

#include <pthread.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/syscall.h>

#include <gmp.h>

#include "cado_poly.h"
#include "getprime.h"
#include "gzip.h"
#include "macros.h"
#include "mpz_poly.h"
#include "renumber.hpp"
#include "timing.h"
#include "typedefs.h"

#define MAX_PRIMES 255 /* maximal number of factor base primes */
#define MAX_LPRIMES 3  /* maximal number of large primes */
#define DEGF_MAX 6
#define REPS 10


#define MAX_THREADS 128
#define MAX_LINE_LEN 8192
#define END_OF_REL_MARKER "\xFF"

#define INTERLOCKED_EXCHANGE(x, y) \
    __sync_bool_compare_and_swap(&x, x, y)

#define INTERLOCKED_INC(x) \
        __sync_fetch_and_add(&x, 1)

#define gettid() syscall(__NR_gettid)

#define FORMAT_CADO    0 /* CADO format (default output format) */
#define FORMAT_FK      1 /* Franke-Kleinjung's format */
#define FORMAT_GGNFS   2 /* GGNFS format (default input format) */
#define FORMAT_CWI     3 /* CWI format */
#define FORMAT_INDEXED 4 /* CADO index of prime format */
static const char *format_names[] = {
  "CADO", "Franke-Kleinjung", "GGNFS", "CWI", "CADO Renumbered"
};


static int line_number = 0;
static int verbose = 0;
static int skip_freerel_primes = 0;
static int num_threads = 1;

static uint64_t relations_written = 0;

/* used for counting time in different processes */
static timingstats_dict_t stats = {0};

struct worker {
    pid_t pid = 0;
    int fd[2] = {0,};

    char readbuf[16 * 1024] = {0,};
    int bufsize = 0;
    int readpos = 0;
};

struct thread_rel_args {
    worker* workers;
    int num_workers;
    FILE* fp;
    int lines_per_relation;
};

struct relation
{
  int64_t a;
  uint64_t b;
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
};

typedef int (*read_relation_func)(FILE* fp, relation * rel, struct relation_data* data);
typedef void (*print_relation_func)(FILE* fp, relation * rel, struct relation_data* data);

struct relation_data {
    mpz_t* f;
    int degf;

    int32_t *rfb;
    int32_t *afb;

    const int in_alg_first;
    const int out_alg_first;

    cado_poly_ptr cpoly;
    const renumber_t * renumber;

    read_relation_func read_relation;
    print_relation_func print_relation;

    const char * relation_filename;

    int iformat;
    int oformat;

};

static uint32_t
get_uint32 (FILE *fp)
{
  uint32_t w;
  size_t const ok = fread (&w, sizeof (uint32_t), 1, fp);
  ASSERT_ALWAYS (ok != 0);
  return w;
}

static int32_t
get_int32 (FILE *fp)
{
  int32_t w;
  size_t const ok = fread (&w, sizeof (int32_t), 1, fp);
  ASSERT_ALWAYS (ok != 0);
  return w;
}

static int get_format_lines(int iformat) {
    switch (iformat) {
    case FORMAT_CADO:
    case FORMAT_CWI:
    case FORMAT_GGNFS:
    case FORMAT_INDEXED:
        return 1;
    case FORMAT_FK:
        return 3;
    default:
        fprintf(stderr, "Unknown format!\n");
        exit(1);
    }
}

/* prints a large prime, given its low and high 32-bit parts.
   Currently not used. */
#if 0
static void
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
#endif

/* Return 1 if the largest prime in rel is < 2^lpb or if lpb == 0, 
   return 0 otherwise */
static int
checksize_relation_cado (relation  *rel, int32_t *rfb, int32_t *afb, int lpb)
{
  unsigned int i;
  unsigned long p;

  if (lpb == 0)
    return 1;

  /* rational side */
  for (i = 0; i < rel->rfb_entries; i++)
    {
      ASSERT(rel->rexp[i] > 0);
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
      ASSERT (rel->aexp[i] > 0);
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
static void
print_relation_cado (FILE* fp, relation  *rel, relation_data * data)
{

  unsigned int i, j;
  int need_comma;
  int32_t *rfb = data->rfb;
  int32_t *afb = data->afb;

  fprintf (fp, "%" PRId64 ",%" PRIu64 ":", rel->a, rel->b);

  if (rel->b == 0 && skip_freerel_primes) {
      fprintf(fp, "\n");
      fflush (fp);
      return;
  }
  need_comma = 0;

  /* rational side */
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      {
	if (need_comma++) fputc (',', fp);
        if (rfb == NULL) /* primes are given directly */
          fprintf (fp, "%lx", rel->rprimes[i]);
        else /* primes are given by their index */
          fprintf (fp, "%x", rfb[rel->rprimes[i]]);
        fflush (fp);
      }
  for (i = 0; i < rel->num_lrp; i++)
    {
      if (need_comma++) putchar (',');
      fprintf(fp, "%lx", rel->large_rprimes[i]);
    }

  fputc (':', fp);
  need_comma = 0;

  /* algebraic side */
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      {
	if (need_comma++) fputc (',', fp);
        if (afb == NULL) /* primes are given directly */
          fprintf(fp, "%lx", rel->aprimes[i]);
        else /* primes are given by their index */
          fprintf(fp, "%x", afb[rel->aprimes[i]]);
      }
  for (i = 0; i < rel->sp_entries; i++)
    {
      if (rel->sexp[i] > 0)
	for (j = 0; j < rel->sexp[i]; j++)
	  {
	    if (need_comma++) fputc (',', fp);
	    fprintf(fp, "%lx", rel->sprimes[i]);
	  }
    }
  for (i = 0; i < rel->num_lap; i++)
    {
      if (need_comma++) fputc (',', fp);
      fprintf(fp, "%lx", rel->large_aprimes[i]);
    }
  fprintf(fp, "\n");
}

/* output a relation in Franke/Kleinjung format:
   - the first line is "F 0 X 5 1" where the two last integers are the degrees
   - then it comes by groups of 3 consecutive lines (W gives a and b in
     hexadecimal form, X gives the algebraic factors, Y the linear factors.
   W 2335e65b 1d9
   X 96f8323 1413403 1369 16B5 1A0C77 347 55D 3 5 65 71 9D F42439
   Y 444B A6CF C827 2C8627 B 25 209 527 B 2 2 2
*/
static void
print_relation_fk (FILE* fp, relation  *rel, relation_data * data)
{
  unsigned int i, j;
  int32_t *rfb = data->rfb;
  int32_t *afb = data->afb;
  int const iformat = data->iformat;

  if (rel->a > 0)
    fprintf(fp, "W %" PRIx64 " %" PRIx64 "\n", rel->a, rel->b);
  else
    fprintf(fp, "W %" PRIx64 " %" PRIx64 "\n", -rel->a, rel->b);
  /* algebraic side */
  fprintf(fp, "X");
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      fprintf(fp, " %x", (iformat == FORMAT_GGNFS) ? afb[rel->aprimes[i]]
              : (int32_t) rel->aprimes[i]);
  for (i = 0; i < rel->sp_entries; i++)
    for (j = 0; j < rel->sexp[i]; j++)
      fprintf(fp, " %lx", rel->sprimes[i]);
  for (i = 0; i < rel->num_lap; i++)
    fprintf(fp, " %lx", rel->large_aprimes[i]);
  fprintf(fp, "\n");
  /* rational side */
  fprintf(fp, "Y");
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      fprintf(fp, " %x", (iformat == FORMAT_GGNFS) ? rfb[rel->rprimes[i]]
              : (int32_t) rel->rprimes[i]);
  for (i = 0; i < rel->num_lrp; i++)
    fprintf(fp, " %lx", rel->large_rprimes[i]);
  fprintf(fp, "\n");
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
  ASSERT (i >= '0');
  if (i <= '9')
    return i - '0';
  ASSERT (i >= 'A');
  if (i <= 'Z')
    return i - 'A' + 10;
  ASSERT (i >= 'a');
  if (i <= 'z')
    return i - 'a' + 36;
  abort();
}


static void 
print_primes_cwi (FILE* fp, unsigned long *primes, unsigned long *exp,
                  int32_t *fb, const unsigned int np)
{
  unsigned long i,j;

  for (i = 0; i < np; i++)
    for (j = 0; j < (exp != NULL ? exp[i] : 1); j++)
      {
        if (fb == NULL) /* primes are given directly */
          fprintf (fp, " %lu", primes[i]);
        else /* primes are given by their index */
          fprintf (fp, " %u", fb[primes[i]]);
      }
}

static void
print_relation_cwi (FILE* fp, relation  *rel, relation_data * data)
{
  unsigned long i, nr, na;
  int32_t *rfb = data->rfb;
  int32_t *afb = data->afb;
  const int alg_first = data->out_alg_first;
  
  nr = 0;
  na = 0;
  for (i = 0; i < rel->rfb_entries; i++)
    nr += rel->rexp[i];
  for (i = 0; i < rel->afb_entries; i++)
    na += rel->aexp[i];
  for (i = 0; i < rel->sp_entries; i++)
    na += rel->sexp[i];
  
  fprintf (fp, "01%c%c %" PRId64 " %" PRIu64,
              cwi_idx_to_char (alg_first ? na : nr),
              cwi_idx_to_char (alg_first ? nr : na),
              rel->a,
              rel->b);

  for (i = 0; i < 2; i++)
    {
      /* if alg_first == 0, do rational primes for side == 0 */
      if (i == (alg_first ? 1 : 0))
        {
          /* rational side */
          print_primes_cwi (fp, rel->rprimes, rel->rexp, rfb, rel->rfb_entries);
          print_primes_cwi (fp, rel->large_rprimes, NULL, NULL, rel->num_lrp);
        }
      else
        {
          /* algebraic side and special primes */
          print_primes_cwi (fp, rel->aprimes, rel->aexp, afb, rel->afb_entries);
          print_primes_cwi (fp, rel->sprimes, rel->sexp, NULL, rel->sp_entries);
          print_primes_cwi (fp, rel->large_aprimes, NULL, NULL, rel->num_lap);
        }
    }
  fprintf (fp, "%c\n", rel->end_of_set ? ';' : ':');

  fflush (fp);
}

/* norm <- |f(a,b)| = |f[d]*a^d + ... + f[0]*b^d| */
static void
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
  
  ASSERT (*n < MAX_PRIMES);
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
static int
read_relation_cado (FILE *fp, relation  *rel, relation_data * data MAYBE_UNUSED)
{
  int c;
  unsigned long p;

  if (fscanf (fp, "%" SCNd64 ",%" SCNu64 ":", &(rel->a), &(rel->b)) != 2)
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

static int 
fk_read_line (const char * lp, const int maxlen, FILE *fp, const char *file)
{
  int i, skip;

  do {
    if (fgets ((char*) lp, maxlen, fp) == NULL)
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

static int 
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

static int
read_relation_fk (FILE *fp, relation  *rel, relation_data * data)
{
  char line[512];
  char *lp;
  const char * file = data->relation_filename;

  /* Read the "W" line with the a and b values */
  do
  {
    if (fk_read_line (line, 512, fp, file) == 0)
      return -1; /* Signal EOF */
  } while (line[0] == '\n');

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
static int
read_relation_cwi (FILE *fp, relation  *rel, relation_data * data)
{
  int ret, side;
  unsigned long p;
  unsigned int i; /* number of given primes */
  char c, flag[5];
  const int alg_first = data->in_alg_first;

  ret = fscanf (fp, "%4s %" SCNd64 " %" SCNu64 "", flag, &(rel->a),
          &(rel->b));
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
      int const rat_now = side == (alg_first ? 1 : 0);
      unsigned int *fb_entries = 
        rat_now ? &rel->rfb_entries : &rel->afb_entries;
      unsigned long *primes = rat_now ? rel->rprimes : rel->aprimes;
      unsigned long *exps = rat_now ? rel->rexp : rel->aexp;
      const char *side_name = rat_now ? "rational" : "algebraic";
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
      fprintf (stderr, "Error, invalid relation for a=%" PRId64 
                       " b=%" PRIu64 "\n", rel->a, rel->b);
      exit (1);
    }

  rel->end_of_set = (c == ';');

  return 1;
}

/* Read one relation in GGNFS format (adapted from function dataConvertToRel
   in GGNFS, file rels.c).
   Return value: 1 if relation is ok, 0 if invalid, -1 if end of file.
*/
static int
read_relation_ggnfs (FILE *fp, relation  *rel, relation_data * data)
{
  uint32_t size_field; /* 32 bits */
  unsigned int i = 0; /* gcc 4.1.2 warns that i may be used uninitialized,
                         but this does not seem correct */
  unsigned int k;
  unsigned sp_entries; /* Number of sprimes according to file */
  int32_t p;
  mpz_t *f = data->f;
  int const degf = data->degf;
  int32_t *afb = data->afb;
  mpz_t norm;

  mpz_init (norm);

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
                fprintf (stderr, "Warning, f(%" PRId64 ",%" PRIu64 ") not "
                    "divisible by factor base prime %d\n", rel->a, rel->b, p);
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
        fprintf (stderr, "Warning, f(%" PRId64 ",%" PRIu64 ") not divisible "
                         "by large prime %lu\n", rel->a, rel->b, p);
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

      mpz_clear (norm);
      return 0;
    }

  rel->end_of_set = 1;

  mpz_clear (norm);
  return 1; /* valid relation */
}

/**
 * Bluntly stolen from check_rels, should be refactored.
 * renumbered relations do not contain all the prime factors for the relation, we use this function to add the missing primes.
 */
static int fix_relation(relation  *rel, cado_poly_ptr cpoly, unsigned int * lpb)
{
  mpz_t norm[2];
  int err = 0;
  unsigned long lpb_max[2];

  for(int side = 0 ; side < 2 ; side++) {
      mpz_init (norm[side]);
      lpb_max[side] = 1UL << lpb[side];
  }

  /* compute the norm on alg and rat sides */
  for(int side = 0 ; side < cpoly->nb_polys ; side++)
  {
      ASSERT_ALWAYS(side < 2);
      mpz_poly_ptr ps = cpoly->pols[side];
      mpz_poly_homogeneous_eval_siui (norm[side], ps, rel->a, rel->b);
      mpz_abs (norm[side], norm[side]);
  }

  unsigned int* entries = &rel->rfb_entries;
  unsigned long* primes = rel->rprimes;
  unsigned long* exp = rel->rexp;

  for (int side = 0 ; side < 2;
          side++,
          entries = &rel->afb_entries,
          primes = rel->aprimes,
          exp = rel->aexp) {

      /* check for correctness of the factorization of the norms */
      for(weight_t i = 0; i < *entries; i++)
      {
        p_r_values_t const p = primes[i];
        exponent_t const e = exp[i];
        ASSERT_ALWAYS(p != 0); /* could reveal a problem in parsing */
        ASSERT_ALWAYS(e > 0); /* non positive exponent is not possible */
        for (int j = 0; j < e; j++)
        {
          if (!mpz_divisible_ui_p (norm[side], p))
          {
            err |= -1;
          }
          else
            mpz_divexact_ui (norm[side], norm[side], p);
        }
      }

      if (mpz_cmp_ui (norm[side], 1) != 0)
      {
          /* complete at least for primes up to 10000 (because of GGNFS and Msieve
           * that skip these primes) */
          unsigned long const max_p = MAX(lpb_max[side], 10000);
          prime_info pi;
          prime_info_init (pi);
          for (unsigned long p = 2; mpz_cmp_ui (norm[side], 1) != 0 && p < max_p ;
               p = getprime_mt (pi))
          {
              while (mpz_divisible_ui_p (norm[side], p))
              {
                mpz_divexact_ui (norm[side], norm[side], p);
                add_prime(primes, exp, entries, p);
              }
          }
          prime_info_clear (pi);

          if (mpz_cmp_ui (norm[side], 1) != 0) {
              err |= -2;
          }
    }
  }

  for(int side = 0 ; side < 2 ; side++)
    mpz_clear(norm[side]);

  return err;
}


/* Read one relation in CADO format from fp, and put it in rel.
   Return 1 if relation is valid, 0 if invalid, -1 if end of file. */
static int
read_relation_renumbered (FILE *fp, relation  *rel, relation_data * data)
{
  int c;
  cado_poly_ptr cpoly = data->cpoly;
  const renumber_t * renumber_table = data->renumber;

  if (fscanf (fp, "%" SCNx64 ",%" SCNx64 ":", &(rel->a), &(rel->b)) != 2)
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

  rel->num_lrp = 0;
  rel->rfb_entries = 0; /* number of rational primes */
  rel->afb_entries = 0; /* number of algebraic primes */

  rel->sp_entries = 0;
  rel->num_lap = 0;
  rel->end_of_set = 1;

  do {
    unsigned long i;
    renumber_t::p_r_side x;
    if (fscanf (fp, "%lx", &i) == 1) {
        x = renumber_table->p_r_from_index(i);
        if (x.side == 0) {
            add_prime (rel->rprimes, rel->rexp, &rel->rfb_entries, x.p);
        }
        else if (x.side == 1) {
            add_prime (rel->aprimes, rel->aexp, &rel->afb_entries, x.p);
        }
        else {
            fprintf (stderr, "Got unexpected side: %d", x.side);
            exit(1);
        }
    }
  } while ((c = getc (fp)) == ',');

  if (c != '\n')
  {
      fprintf (stderr, "Error, invalid relation (expected newline)\n");
      exit(1);
  }

  if (rel->b == 0 && skip_freerel_primes) {
      // Avoid fix_relation when skip_freerel_primes is set.
      return 1;
  }

  /* Add missing primes - Not they will be added out-of-order */
  std::vector<unsigned int> lpb;
  lpb.reserve(renumber_table->get_nb_polys());
  for(int side = 0 ; side < renumber_table->get_nb_polys() ; side++)
      lpb.push_back(renumber_table->get_lpb(side));
  fix_relation(rel, cpoly, lpb.data());

  return 1;
}


/* Read relations from file fp, with input format is 'iformat'. 
   and prints them with given output format 'oformat'.

   degf is the degree of the algebraic polynomial (needed for GGNFS
   input formats).
   
   For GGNFS input: the rational factor base is into rfb, the algebraic factor
   base is into afb, and the coefficients of the algebraic polynomial are into
   f[0]...f[degf].
*/
static unsigned long
convert_relations (FILE* infp, FILE* outfp, int lpb, int multi, unsigned int rels_in_file, relation_data * data)
{
  relation  rel[1];
  unsigned long output_rels = 0;
  unsigned int j;
  int ok;

  line_number = 0;

  for (j = 0; feof (infp) == 0 && j < rels_in_file; j++)
    {
      ok = data->read_relation(infp, rel, data);

      if (ok == -1) /* end of file */
        break;

      if (ok != 0 && checksize_relation_cado (rel, data->rfb, data->afb, lpb))
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

          data->print_relation(outfp, rel, data);
          output_rels++;
        }

      // Mark end of relation
      if (num_threads > 1)
          fputs(END_OF_REL_MARKER, outfp);

      fflush (outfp);
    }

  fclose(infp);

  return output_rels;
}

#define BF_DELIMITER 0x0a
#define END_POLY "END_POLY"
#define END_HEADER "END_HEADER"

/* read fp until the string s is found,
   the next char read by getc (fp) is the one following s. */
static void
get_string (const char *s, FILE *fp, int verbose)
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
static int
read_fb (FILE *fp, int32_t **rfb, int32_t *rfb_size, int32_t **afb, 
         int32_t *afb_size, int verbose, mpz_t *f)
{
  int c, i, degf;
  long RFBsize, AFBsize;
  int retscanf MAYBE_UNUSED;

  /* read polynomial */
  get_string ("Y1: ", fp, verbose);
  mpz_inp_str (f[0], fp, 10);
  if (verbose)
    mpz_out_str (stderr, 10, f[0]);
  c = getc (fp); /* newline */
  retscanf = fscanf (fp, "c%d: ", &degf);
  ASSERT (retscanf == 1);
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
      retscanf = fscanf (fp, "c%d: ", &c);
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
  retscanf = fscanf (fp, "%ld", &RFBsize);
  if (verbose)
    fprintf (stderr, "%ld", RFBsize);
  get_string ("AFBsize: ", fp, verbose);
  retscanf = fscanf (fp, "%ld", &AFBsize);
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
usage (const char *s)
{
  fprintf (stderr, "Usage: %s [options] f0 f1 ...\n", s);
  fprintf (stderr, "       where f0 f1 ... are input files with relations\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -v       - verbose\n");
  fprintf (stderr, "       -if xxx  - specifies input format (default ggnfs)\n");
  fprintf (stderr, "       -of yyy  - specifies output format (default cado)\n");
  fprintf (stderr, "                  where xxx, yyy are in {cado, fk, cwi, cwi2, ggnfs, renumbered}\n");
  fprintf (stderr, "                  cwi has rational side first, cwi2 algebraic first\n");
  fprintf (stderr, "       -fb file - factor base (needed for ggnfs input)\n");
  fprintf (stderr, "       -deg d   - algebraic degree (for fk output)\n");
  fprintf (stderr, "       -lpb l   - discard relations with primes >= 2^l\n");
  fprintf (stderr, "       -nomulti - print repeated primes only once\n");
  fprintf (stderr, "       -skip-freerel-primes - do not print primes for free relations (for msieve)\n");

  fprintf (stderr, "       -poly file - the polynomial file (for -if renumbered)\n");
  fprintf (stderr, "       -renumber file - the renumber file (for -if renumbered)\n");
  fprintf (stderr, "       -out file - output file (or stdout)\n");
  fprintf (stderr, "       -t NUM - Number of threads to run\n");

}

static void* read_rels(void* _args) {
    thread_rel_args * args = (thread_rel_args *) _args;
    worker * workers = args->workers;
    int const num_workers = args->num_workers;
    FILE* fp = args->fp;

    char line[MAX_LINE_LEN];
    int len = 0;
    uint32_t i = 0;

    if (verbose) fprintf(stderr, "# read_rels() started\n");

    for (i = 0; fgets(line, sizeof(line) / sizeof(line[0]), fp) != NULL; i++) {
        len = strlen(line);
        if (write(workers[i % num_workers].fd[0], line, len) != len) {
            perror("write");
            exit(1);
        }
    }

    fprintf(stderr, "Finished reading %u lines\n", i);

    // Notify workers with EOF
    for (i = 0; i < (uint32_t)num_workers; i++) {
        shutdown(workers[i % num_workers].fd[0], SHUT_WR);
    }

    return NULL;
}


static void* read_ggnfs_rels(void* _args) {
    thread_rel_args * args = (thread_rel_args *) _args;
    worker * workers = args->workers;
    int const num_workers = args->num_workers;
    FILE* fp = args->fp;

    uint32_t i = 0;
    size_t ok = 0;
    size_t written = 0;

    if (verbose) fprintf(stderr, "# read_ggnfs_rels() started\n");

    for (i = 0; !feof(fp); i++) {
        int const fd = workers[i % num_workers].fd[0];
        uint32_t size_field;
        uint32_t rfb_entries, afb_entries, sp_entries, num_lrp, num_lap;

        char buf[4096];
        size_field = get_uint32(fp);
        written = write(fd, &size_field, sizeof(size_field));
        ASSERT_ALWAYS (written == sizeof(size_field));

        ok = fread(buf, sizeof(int32_t), 2, fp);
        ASSERT_ALWAYS (ok == 2);
        written = write(fd, buf, sizeof(int32_t) * 2);
        ASSERT_ALWAYS (written == sizeof(int32_t) * 2);

        rfb_entries = size_field >> 24;
        ok = fread(buf, sizeof(uint32_t) * 2, rfb_entries, fp);
        ASSERT_ALWAYS (ok == rfb_entries);
        written = write(fd, buf, sizeof(uint32_t) * 2 * rfb_entries);
        ASSERT_ALWAYS (written == sizeof(uint32_t) * 2 * rfb_entries);

        afb_entries = (size_field >> 16) & 255;
        ok = fread(buf, sizeof(uint32_t) * 2, afb_entries, fp);
        ASSERT_ALWAYS (ok == afb_entries);
        written = write(fd, buf, sizeof(uint32_t) * 2 * afb_entries);
        ASSERT_ALWAYS (written == sizeof(uint32_t) * 2 * afb_entries);

        sp_entries = (size_field >> 8) & 255;
        ok = fread(buf, sizeof(uint32_t) * 2, sp_entries, fp);
        ASSERT_ALWAYS (ok == sp_entries);
        written = write(fd, buf, sizeof(uint32_t) * 2 * sp_entries);
        ASSERT_ALWAYS (written == sizeof(uint32_t) * 2 * sp_entries);

        // qcb1, qcb2
        ok = fread(buf, sizeof(int32_t), 2, fp);
        ASSERT_ALWAYS (ok == 2);
        written = write(fd, buf, sizeof(int32_t) * 2);
        ASSERT_ALWAYS (written == sizeof(int32_t) * 2);

        num_lrp = (size_field >> 6) & 3;
        ok = fread(buf, sizeof(uint32_t), num_lrp, fp);
        ASSERT_ALWAYS (ok == num_lrp);
        written = write(fd, buf, sizeof(uint32_t) * num_lrp);
        ASSERT_ALWAYS (written == sizeof(uint32_t) * num_lrp);

        num_lap = (size_field >> 4) & 3;
        ok = fread(buf, sizeof(uint32_t) * 2, num_lap, fp);
        ASSERT_ALWAYS (ok == num_lap);
        written = write(fd, buf, sizeof(uint32_t) * 2 * num_lap);
        ASSERT_ALWAYS (written == sizeof(uint32_t) * 2 * num_lap);
    }

    fprintf(stderr, "Finished reading %u relations\n", i);

    // Notify workers with EOF
    for (i = 0; i < (uint32_t)num_workers; i++) {
        shutdown(workers[i % num_workers].fd[0], SHUT_WR);
    }

    return NULL;
}

#define READ_EOF     -1
#define READ_PARTIAL  0
#define READ_COMPLETE 1

static int read_until(worker * worker, char** buf, int* bufsize, const char* needle) {
    int ret = 0;
    if (worker->readpos >= worker->bufsize) {
        ret = read(worker->fd[0], worker->readbuf, sizeof(worker->readbuf));

        worker->readpos = 0;
        worker->bufsize = ret;

        if (ret == 0) {
            return READ_EOF;
        }
    }

    *buf = &worker->readbuf[worker->readpos];
    *bufsize = worker->bufsize - worker->readpos;

    const char* ptr = (const char *) memmem(&worker->readbuf[worker->readpos], *bufsize, needle, strlen(needle));
    if (ptr == NULL) {
        // Consume all buffer
        worker->readpos = worker->bufsize;
        return READ_PARTIAL;
    }

    // Set the buffer size to end at ptr
    *bufsize = ptr - &worker->readbuf[worker->readpos];

    // Advance the position
    worker->readpos += (ptr - &worker->readbuf[worker->readpos] + 1);

    return READ_COMPLETE;
}


static void* write_rels(void* _args) {
    thread_rel_args * args = (thread_rel_args *) _args;
    worker * workers = args->workers;
    int const num_workers = args->num_workers;
    FILE* fp MAYBE_UNUSED = args->fp;

    if (verbose) fprintf(stderr, "# write_rels() started\n");
    char* buf = NULL;
    int bufsize = 0;
    int ret = READ_EOF;

    for (;;) {
        int any_alive = 0;
        for (int i = 0; i < num_workers; i++) {
            do {
                ret = read_until(&workers[i % num_workers],
                                                   &buf,
                                                   &bufsize,
                                                   END_OF_REL_MARKER);
                if (bufsize <= 0) {
                    // Empty line
                    continue;
                }
                if (fwrite(buf, 1, bufsize, fp) != (size_t)bufsize) {
                    perror("fwrite");
                    exit(1);
                }
            } while (ret == READ_PARTIAL);

            any_alive |= (ret != READ_EOF);
            INTERLOCKED_INC(relations_written);
        }

        if (!any_alive) {
            break;
        }
    }

    if (verbose) fprintf(stderr, "# Finished processing %lu relations\n", relations_written);
    return NULL;
}

int main(int argc, char const * argv[])
{
  FILE *fp;
  int32_t *rfb = NULL, *afb = NULL;
  int32_t rfb_size, afb_size;
  int i, j, degf = 0;
  int iformat = FORMAT_GGNFS; /* default input format */
  int oformat = FORMAT_CADO; /* default output format */
  mpz_t f[DEGF_MAX + 1];
  int num_files = 0;
  const char * program_name = argv[0];
  int lpb = 0; /* 0 means no bound */
  int multi = 1, in_alg_first = 0, out_alg_first = 0;

  uint32_t rels_in_file = UINT32_MAX;

  worker  workers[MAX_THREADS];

  cado_poly cpoly;
  renumber_t renumber_table;

  read_relation_func read_relation = NULL;
  print_relation_func print_relation = NULL;

  const char * fbfile = NULL; /* name of factor base file */
  const char * relsfile = NULL; /* name of binary file to read */
  const char * polyfile = NULL; /* The polynomial file */
  const char * renumberfile = NULL; /* The renumbered file */
  const char * outfile = NULL;

  FILE* out_fp = stdout;

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
          else if (strcmp (argv[2], "renumbered") == 0)
            iformat = FORMAT_INDEXED;
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
      fbfile = argv[2];
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
      else if (argc > 2 && strcmp (argv[1], "-skip-freerel-primes") == 0)
    {
      skip_freerel_primes = 1;
      argv ++;
      argc --;
    }
      else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
    {
      polyfile = argv[2];
      argv += 2;
      argc -= 2;
    }
      else if (argc > 2 && strcmp (argv[1], "-renumber") == 0)
    {
      renumberfile = argv[2];
      argv += 2;
      argc -= 2;
    }
      else if (argc > 2 && strcmp (argv[1], "-t") == 0)
    {
      num_threads = atoi (argv[2]);
      if (num_threads <= 1) {
          num_threads = 1;
      }
      argv += 2;
      argc -= 2;
    }
      else if (argc > 2 && strcmp (argv[1], "-out") == 0)
    {
      outfile = argv[2];
      argv += 2;
      argc -= 2;
    }
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (1);
	}
    }

  /* for ggnfs input, we need the factor base,
     and for Franke-Kleinjung output, we need the algebraic degree */
  if ((iformat == FORMAT_GGNFS && fbfile == NULL) || argc <= 1)
    {
      usage (program_name);
      exit (1);
    }

  if (outfile != NULL) {
      out_fp = fopen_maybe_compressed(outfile, "w");
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


  printf ("# Verbose output: %s\n", verbose ? "yes" : "no");
  if (verbose) {

      printf ("# multi: %s\n", multi ? "yes" : "no");
      printf ("# skip freerel primes: %s\n", skip_freerel_primes ? "yes" : "no");
      printf ("# number of workers: %d\n", num_threads);
      printf ("# input format: %s\n", format_names[iformat]);
      printf ("# output format: %s\n", format_names[oformat]);
  }

  timingstats_dict_init(stats);

  /* check we can open the factor base for GGNFS input format */
  if (iformat == FORMAT_GGNFS)
    {
      fp = fopen (fbfile, "rb");
      if (fp == NULL)
        {
          fprintf (stderr, "Error, unable to open factor base file %s\n", fbfile);
          exit (1);
        }
      degf = read_fb (fp, &rfb, &rfb_size, &afb, &afb_size, verbose, f);
      fclose (fp);
    }

  if (iformat == FORMAT_INDEXED) {
      if (renumberfile == NULL) {
          fprintf(stderr, "Missing -renumber for -if renumbered!\n");
          exit(1);
      }
      if (polyfile == NULL) {
          fprintf(stderr, "Missing -poly for -if renumbered!\n");
          exit(1);
      }

      cado_poly_init(cpoly);
      if (!cado_poly_read (cpoly, polyfile))
      {
          fprintf (stderr, "Error reading polynomial file\n");
          exit(1);
      }

      renumber_table = renumber_t(cpoly);
      renumber_table.read_from_file(renumberfile, 0);
  }

  switch (iformat) {
  case FORMAT_GGNFS:
      read_relation = read_relation_ggnfs;
      break;
  case FORMAT_CADO:
      read_relation = read_relation_cado;
      break;
  case FORMAT_FK:
      read_relation = read_relation_fk;
      break;
  case FORMAT_CWI:
      read_relation = read_relation_cwi;
      break;
  case FORMAT_INDEXED:
      read_relation = read_relation_renumbered;
      break;
  default:
      fprintf (stderr, "Error, unknown format %d\n", iformat);
      exit (1);
  }
  switch (oformat)
  {
   case FORMAT_CADO:
     print_relation = print_relation_cado;
     break;
   case FORMAT_FK:
     print_relation = print_relation_fk;
     break;
   case FORMAT_CWI:
     print_relation = print_relation_cwi;
     break;
   default:
     fprintf (stderr, "Error, unknown format %d\n", oformat);
     exit (1);
  }


  relation_data  data = {
          .f = f,
          .degf = degf,

          .rfb = rfb,
          .afb = afb,

          .in_alg_first = in_alg_first,
          .out_alg_first = out_alg_first,

          .cpoly = cpoly,
          .renumber = & renumber_table,

          .read_relation = read_relation,
          .print_relation = print_relation,

          .relation_filename = relsfile,

          .iformat = iformat,
          .oformat = oformat
  };

  while (argc > 1)
    {
      relsfile = argv[1];
      argv ++;
      argc --;
      /* if output format is Franke-Kleinjung, print degree if provided */
      if (oformat == FORMAT_FK && num_files == 0 && degf != 0)
          fprintf (out_fp, "F 0 X %d 1\n", degf);

      num_files++;

      if (verbose) fprintf(stderr, "# Opening file: %s\n", relsfile);

      fp = fopen_maybe_compressed (relsfile, "r");
      if (fp == NULL)
      {
          fprintf (stderr, "Error, unable to open relation file %s\n", relsfile);
          exit (1);
      }

      if (iformat == FORMAT_GGNFS)
      {
          size_t retfread MAYBE_UNUSED;
          retfread = fread (&rels_in_file, sizeof (int32_t), 1, fp);
          ASSERT (retfread != 0);
          if (verbose)
            fprintf (stderr, "File %s: header says %d relations.\n",
                     relsfile, rels_in_file);

          /* the first 4 bytes give the number of relations in the file */
          rewind (fp);
      }

      if (num_threads <= 1) {
          unsigned long const num_rels = convert_relations (fp, out_fp, lpb, multi, rels_in_file, &data);
          if (verbose) fprintf(stderr, "# Converted %lu relations from file: %s!\n", num_rels, relsfile);
          continue;
      }

      if (verbose) fprintf(stderr, "# Spawning %d workers\n", num_threads);

      int* fd = workers[i].fd;

      for (i = 0; i < num_threads; i++) {
          fd = workers[i].fd;

          if ( socketpair( AF_UNIX, SOCK_STREAM, 0, fd ) < 0 ) {
            perror("Cannot open socketpair");
            exit(1);
          }

          /**
           * Setting a big buffer is very important. otherwise a deadlock might occur.
           * Consider the following:
           * read_rels is stuck on write to worker 0 (since it has not read stdin in awhile)
           * worker 0 is stuck on write to write_rels
           * - At this point other workers do do not get additional rels.
           * write_rels is stuck on read from a different worker.
           *
           * To avoid this situation the buffer must cover the time to iterate over all of the workers
           * from read_rels and from write_rels
           */
          int sendbuf = 1024 * 1024;
          setsockopt(fd[0], SOL_SOCKET, SO_SNDBUF, &sendbuf, sizeof(sendbuf));
          setsockopt(fd[1], SOL_SOCKET, SO_SNDBUF, &sendbuf, sizeof(sendbuf));

          int const ret = fork();
          if (ret == -1) {
              perror("fork()");
              exit (1);
          } else if (ret == 0) {
              if (fd[1] != STDIN_FILENO) { /*Redirect standard input to socketpair*/
                if (dup2(fd[1], STDIN_FILENO) != STDIN_FILENO) {
                  perror("Cannot dup2 stdin");
                  exit(1);
                }
              }
              if (fd[1] != STDOUT_FILENO) { /*Redirect standard output to socketpair*/
                if (dup2(fd[1], STDOUT_FILENO) != STDOUT_FILENO) {
                  perror("Cannot dup2 stdout");
                  exit(1);
                }
              }
              close(fd[0]);
              close(fd[1]);
              fd[0] = -1;
              fd[1] = -1;
              break;
          }
          else {
              // parent
              workers[i].pid = ret;
              if (verbose) fprintf(stderr, "# Worker-%d (tid=%d) fd=%d\n", i, ret, fd[0]);
              close(fd[1]);
              fd[1] = -1;
              fd = NULL;
          }
      }

      if (fd == NULL) {
          pthread_t reader, writer;
          thread_rel_args  reader_args = {
                  .workers = workers,
                  .num_workers = num_threads,
                  .fp = fp,
                  .lines_per_relation = get_format_lines(iformat)
          };
          thread_rel_args  writer_args = {
                .workers = workers,
                .num_workers = num_threads,
                .fp = out_fp,
                .lines_per_relation = 0 //unused for writing
          };
          void* retval;
          if (pthread_create(&reader, NULL, iformat == FORMAT_GGNFS ? read_ggnfs_rels : read_rels, &reader_args) != 0) {
              perror("thread_create(reader)");
              exit(1);
          }
          if (pthread_create(&writer, NULL, write_rels, &writer_args) != 0) {
              perror("thread_create(writer)");
              exit(1);
          }

          // Wait for workers to finish
          for (j = 0; j < num_threads; j++) {
              int status = 0;
              waitpid(workers[j].pid, &status, 0);
              if (verbose) fprintf(stderr, "# waitpid(%d) = %d\n", workers[j].pid, status);

          }

          pthread_join(reader, &retval);
          pthread_join(writer, &retval);
      }
      else {
          char filename[4096] = {0};

          snprintf(filename, sizeof(filename) - 1, "%s (Worker-%d)", relsfile, i);
          data.relation_filename = filename;
          if (verbose) fprintf(stderr, "# Worker-%d started tid=%ld\n", i, gettid());

          unsigned long const num_rels = convert_relations (stdin, stdout, lpb, multi, rels_in_file, &data);

          // Worker finished successfully
          if (verbose) fprintf(stderr, "# Worker-%d finished converting %lu relations!\n", i, num_rels);
          shutdown(STDIN_FILENO, SHUT_WR);
          close(STDIN_FILENO);

          exit(0);
      }
  }

  fprintf (stderr, "Output %lu relations from %d file(s).\n", relations_written,
	   num_files);

  free (rfb);
  free (afb);
  for (i = 0; i <= DEGF_MAX; i++)
    mpz_clear (f[i]);

  if (outfile != NULL) {
      fclose_maybe_compressed(out_fp, outfile);
  }

  timingstats_dict_add_myprocess(stats, "main");
  timingstats_dict_disp(stats);
  timingstats_dict_clear(stats);

  return 0;
}
