/* This program is an experimental implementation of out-of-core singleton
   removal for CADO-NFS.

   Input: a set of relation files (might be compressed), each one on a single
          line in the "filelist" argument. Each file name can be in
          subdirectories like "foo/bar/toto.rels.gz".
   Output: each file with singletons removed will be written in the "new"
          directory (relative to where the program was started). This directory
          must already exist. In addition, if file names contain subdirectories
          those subdirectories must already exist. For example the output for
          the input file "foo/bar/toto.rels.gz" will be stored in
          "new/foo/bar/toto.rels.gz".

   Parameters:
   -M nnn : gives the size nnn of the hash table used; nnn should be 2*p where
            p is a prime. In addition for efficiency nnn should be at least
            two times larger than the maximal number of primes ideals, which
            is about 2^lpba/log(2^lpba) + 2^lpbr/log(2^lpbr).
   -minpa nnn : considers only algebraic ideals >= nnn
   -minpr nnn : considers only rational ideals >= nnn
   -t nnn     : uses nnn threads

   Memory usage: this program uses only 2 bits per entry in the hash table.
   Here are possible values of -M for different values of lpb[ar], and the
   corresponding memory usage.

   lpba,lpbr         -M nnn       memory usage
     27,27           28686706         7.2Mb
     28,28           55324382        13.8Mb
     29,29          106833262        26.7Mb
     30,30          206544298        51.6Mb
     31,31          399763178        99.9Mb
     32,32          774541058         194Mb
     33,33         1502140186         376Mb
     34,34         2915919074         729Mb
     35,35         5665214234         1.4Gb
     36,36        11015694286         2.8Gb
     37,37        21435945614         5.4Gb
     38,38        41743683574        10.4Gb
     39,39        81346665278        20.3Gb
     40,40       158625997294        39.7Gb

   Each run of the program performs two passes on the relations:
   - pass 1 reads all relations and counts the weight of each ideal
   - pass 2 discards relations with ideals of weight 1

   Since the program uses a hash function which has collisions (the larger the
   -M option, the less collisions), some singletons might not be detected.
   However the program is intended to be run again on the "new" data, with a
   *different* value of the -M option, to split ideals which had collisions in
   previous runs. After say 10 runs of the program with different -M values,
   only very few singletons should remain.
*/

#define _GNU_SOURCE     /* FIXME. Include directly cado.h instead */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>

#define MAXNAME 1024
#define MAX_THREADS 128

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* mutual exclusion lock */
long wct0;

uint64_t M = 0, minpa = 0, minpr = 0, Hsize, *H, nrels = 0, remains = 0;

/* thread structure */
typedef struct
{
  char g[MAXNAME];
  int thread;
  unsigned long nlines;
  unsigned long nideal;
  unsigned long nsingl;
  unsigned long ndupli;
  unsigned long nthreads;
} __tab_struct;
typedef __tab_struct tab_t[1];

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "gzip -dc %s", "gzip -c --fast > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    /* These two have to be present */
    { "", NULL, NULL },
    { NULL, NULL, NULL },
};

long
cputime ()
{
  return (long) seconds();
}

long
realtime ()
{
  struct timeval tv[1];
  /* struct timeval {
     time_t      tv_sec;     seconds
     suseconds_t tv_usec;    microseconds
     }; */

  gettimeofday (tv, NULL);
  return tv->tv_sec;
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

/* The pipe capacity is 2^16 by default, we can increase it, but it does not
   seem to make a difference, thus we don't change it by default (patch from
   Alain Filbois). */
#define PIPE_CAPACITY 1UL << 18

FILE*
fopen_compressed_r (const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_in) { /* this suffix has an associated deflate command */
            char * command;
            int ret = asprintf(&command, r->pfmt_in, name);
            assert(ret >= 0);
            f = popen(command, "r");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            /* remove the 'xxx' to change the pipe capacity */
#ifdef F_SETPIPE_SZ
            fcntl (fileno (f), F_SETPIPE_SZ, PIPE_CAPACITY);
#endif
            return f;
        } else {
            f = fopen(name, "r");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE*
fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_out) { /* suffix has an associated compression command */
            char * command;
            int ret = asprintf(&command, r->pfmt_out, name);
            assert (ret >= 0);
            f = popen(command, "w");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            /* remove the 'xxx' to change the pipe capacity */
#ifdef F_SETPIPE_SZ
            /* increase the pipe capacity (2^16 by default), thanks to
               Alain Filbois */
            fcntl (fileno (f), F_SETPIPE_SZ, PIPE_CAPACITY);
#endif
            return f;
        } else {
            f = fopen(name, "w");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE * gzip_open(const char * name, const char * mode)
{
    if (strcmp(mode, "r") == 0) {
        return fopen_compressed_r (name, NULL, NULL);
    } else if (strcmp(mode, "w") == 0) {
        return fopen_compressed_w (name, NULL, NULL);
    } else {
        abort();
    }
}

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

void gzip_close(FILE * f, const char * name)
{
    const char * e = name + strlen(name);
    const struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        const char * ne = e - MIN(strlen(name), strlen(r->suffix));
        if (strcmp(r->suffix, ne) != 0)
            continue;
        if (r->pfmt_out) {
            pclose(f);
            return;
        } else {
            fclose(f);
            return;
        }
    }
}

static uint64_t
index_rat (uint64_t p)
{
  return (p + 1) % M;
}

static uint64_t
index_alg (uint64_t p, uint64_t r)
{
  return (p + (M - 2) * r) % M; /* always odd for odd p and even M */
}

#if 0
static void
insert (uint64_t h)
{
  uint64_t r, s, v;

  assert (h < M);
  r = h >> 5;        /* h div 32 */
  s = (h & 31) << 1; /* 2*(h mod 32) */
  v = (H[r] >> s) & 3;
  if (v < 3)
    H[r] += (uint64_t) 1 << s;
}
#else
static void
insert (uint64_t h)
{
  uint64_t r, s, v;
  const uint64_t mask[32] = {0x3, 0xc, 0x30, 0xc0, 0x300, 0xc00, 0x3000, 0xc000, 0x30000, 0xc0000, 0x300000, 0xc00000, 0x3000000, 0xc000000, 0x30000000, 0xc0000000, 0x300000000, 0xc00000000, 0x3000000000, 0xc000000000, 0x30000000000, 0xc0000000000, 0x300000000000, 0xc00000000000, 0x3000000000000, 0xc000000000000, 0x30000000000000, 0xc0000000000000, 0x300000000000000, 0xc00000000000000, 0x3000000000000000, 0xc000000000000000};

  assert (h < M);
  r = h >> 5;        /* h div 32 */
  s = h & 31       ; /* h mod 32 */
  v = H[r] & mask[s];
  if (v != mask[s])
    H[r] += (uint64_t) 1 << (s << 1); /* add 1 to bit 2s */
}
#endif

static void
delete (uint64_t h)
{
  uint64_t r, s, v;

  assert (h < M);
  r = h >> 5;        /* h div 32 */
  s = (h & 31) << 1; /* 2*(h mod 32) */
  v = (H[r] >> s) & 3;
  if (v == 0)
    fprintf (stderr, "Warning: deleting ideal with weight 0 at h=%lu\n", h);
  else if (v == 2) /* no need to update for v=1 (singleton) */
    H[r] -= (uint64_t) 1 << s;
}

static int
is_singleton (uint64_t h)
{
  uint64_t r, s, v;

  h = h % M;
  r = h >> 5;
  s = (h & 31) << 1;
  v = (H[r] >> s) & 3;
  if (v == 0)
    fprintf (stderr, "Warning: unexpected zero counter for h=%lu\n", h);
  return v == 1;
}

/* return a*b mod p, assuming 0 <= a, b, p < 2^40 */
uint64_t
mulmod (uint64_t a, uint64_t b, uint64_t p)
{
  uint64_t a1, a0, b1, b0, r;

  a1 = a >> 20;
  a0 = a & 1048575;
  b1 = b >> 20;
  b0 = b & 1048575;
  r = (a1 * b1) % p; /* r < 2^40 */
  r = ((r << 20) + a1 * b0 + a0 * b1) % p;
  r = ((r << 20) + a0 * b0) % p;
  return r;
}

/* return a/b mod p */
static unsigned long
root (long a0, long b0, long p)
{
  long u, w, q, r, t, a, b;

  a = a0 % p;
  b = b0 % p;
  assert (b != 0);
  if (a < 0)
    a += p;

  u = 1;
  w = 0;
  a = b;
  b = p;
  /* invariant: a = u*b0 mod p, b = w*b0 mod p */
  while (b != 0)
    {
      q = a / b;
      r = a % b;
      a = b;
      b = r;
      t = u;
      u = w;
      w = (t - q * w) % p;
    }
  assert (a == 1);
  /* compute a0*u % p */
  if (u < 0)
    u += p;
  assert (u > 0);
  if (a0 >= 0)
    u = mulmod (a0, u, p);
  else
    u = p - mulmod (-a0, u, p);
  assert (u > 0);
  r = ((long) mulmod (u, b0, p) - a0) % p;
  assert (r == 0);
  return u;
}

/* return the number of lines read */
static unsigned long
foo (char *g)
{
  FILE *f;
  char *ret, s[1024], *t, *endptr;
  unsigned long line = 0, b, r;
  long a;
  uint64_t p, pfree;

  f = gzip_open (g, "r");

  while (1)
    {
      ret = fgets (s, 1024, f);
      if (ret == 0)
        break;
      line ++;
      t = s;
      a = 0;
      a = strtol (t, &endptr, 10);
      assert (a != 0);
      assert (endptr > t);
      t = endptr;
      assert (*t == ',');
      t ++;
      b = strtoul (t, &endptr, 10);
      assert (endptr > t);
      t = endptr;
      assert (*t == ':');
      t ++;
      pfree = 0;
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (b == 0)
            pfree = p;
          if (p >= minpr)
            insert (index_rat (p));
          if (*t++ != ',')
            break; /* end of rational primes */
        }
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (pfree >= minpa)
            insert (index_alg (pfree, p));
          else if (p >= minpa)
            {
              r = root (a, b, p);
              insert (index_alg (p, r));
            }
          if (*t++ != ',')
            break; /* end of algebraic primes */
        }
    }

  gzip_close (f, g);
  return line;
}

static void
bar (char *g)
{
  FILE *f, *newf;
  char *ret, s[1024], *t, *endptr, newg[1024];
  unsigned long line = 0, b, r, output = 0;
  long a;
  uint64_t p, pfree;
  int keep, nr, na;
  uint64_t rp[16], ap[16], ar[16];

  f = gzip_open (g, "r");
  strcpy (newg, "new/");
  strcpy (newg + 4, g);
  newf = gzip_open (newg, "w");
  assert (newf != NULL);

  while (1)
    {
      ret = fgets (s, 1024, f);
      if (ret == 0)
        break;
      line ++;
      t = s;
      a = 0;
      a = strtol (t, &endptr, 10);
      assert (a != 0);
      assert (endptr > t);
      t = endptr;
      assert (*t == ',');
      t ++;
      b = strtoul (t, &endptr, 10);
      assert (endptr > t);
      t = endptr;
      assert (*t == ':');
      t ++;
      pfree = 0;
      keep = 1;
      nr = na = 0;
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (b == 0)
            pfree = p;
          if (p >= minpr)
            {
              rp[nr++] = p;
              if (is_singleton (index_rat (p)))
                keep = 0;
            }
          if (*t++ != ',')
            break; /* end of rational primes */
        }
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (pfree >= minpa)
            {
              ap[na] = pfree;
              ar[na++] = p;
              if (is_singleton (index_alg (pfree, p)))
                keep = 0;
            }
          else if (p >= minpa)
            {
              r = root (a, b, p);
              ap[na] = p;
              ar[na++] = r;
              if (is_singleton (index_alg (p, r)))
                keep = 0;
            }
          if (*t++ != ',')
            break; /* end of algebraic primes */
        }
      if (keep)
        {
          fprintf (newf, "%s", s);
          output ++;
        }
      else /* update 2-bit counters of removed ideals */
        {
          while (nr-- > 0)
            delete (index_rat (rp[nr]));
          while (na-- > 0)
            delete (index_alg (ap[na], ar[na]));
        }
    }

  gzip_close (f, g);
  gzip_close (newf, newg);

  fprintf (stderr, "   new/%s done: remains %lu rels out of %lu\n",
           g, output, line);
  fflush (stderr);
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
}

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  fprintf (stderr, "1: Thread %d deals with %s\n", tab[0]->thread, tab[0]->g);
  fflush (stderr);
  tab[0]->nlines = foo (tab[0]->g);
  fprintf (stderr, "   %s done (%lu relations)\n", tab[0]->g, tab[0]->nlines);
  fflush (stderr);
  return NULL;
}

void*
one_thread2 (void* args)
{
  tab_t *tab = (tab_t*) args;

  fprintf (stderr, "2: Thread %d deals with %s\n", tab[0]->thread, tab[0]->g);
  fflush (stderr);
  bar (tab[0]->g);
  return NULL;
}

void*
one_thread3 (void* args)
{
  tab_t *tab = (tab_t*) args;
  uint64_t h, r, wt[4];

  wt[1] = 0;
  wt[2] = 0;
  wt[3] = 0;
  for (r = tab[0]->thread; r < Hsize; r += tab[0]->nthreads)
    {  /* deal with 64 bits, i.e., 32 entries at a time */
      h = H[r];
      while (h != 0)
        {
          /* we could save a factor of 2 with wt[h & 15]++ and a table wt[] of
             16 entries, dealing with 2 counters at a time */
          wt[h & 3] ++;
          h >>= 2;
        }
    }

  tab[0]->nideal = wt[1] + wt[2] + wt[3];
  tab[0]->nsingl = wt[1];
  tab[0]->ndupli = wt[2];

  return NULL;
}

/* given a hash table of n entries, of which ek entries are already filled,
   return an estimate of the number k of insertions in the table.
   We assume that e(k+1) = e(k) + 1 - e(k)/n with e(0) = 1,
   which gives e(k) = n - n*(1-1/n)^k,
   thus k = log(1-e(k)/n)/log(1-1/n) */
double
estimate_ideals (double n, double ek)
{
  return -n * log (1.0 - ek / n);
}

/* estimate the number of prime ideals below a */
static double
est_excess (double a)
{
  return a / log (a);
}

#if 0
/* mono-thread version */
static void
stat (void)
{
  uint64_t nideal = 0, nsingl = 0, h, r;
  long st = cputime (), rt = realtime ();

  /* if a word is 0, skip 32 values of h */
  for (r = 0; r < Hsize; r++)
    {  /* deal with 64 bits, i.e., 32 entries at a time */
      h = H[r];
      while (h != 0)
        {
          nideal += h & 2;
          nsingl += h & 1;
          h >>= 2;
        }
    }
  nideal >>= 1;     /* number of entries with value 3 */
  nsingl -= nideal; /* number of entries with value 1 */
  nideal += nsingl; /* number of entries with value 1 or 3 */
  fprintf (stderr, "stat() took %lds (cpu), %lds (real)\n",
           cputime () - st, realtime () - rt);
  fprintf (stderr, "Total %lu ideals (%.2f%%), %lu singletons\n",
           nideal, 100.0 * (double) nideal / (double) M, nsingl);
  fflush (stderr);
}
#endif

/* multi-thread version */
static void
stat_mt (int nthreads)
{
  int i;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  unsigned long nideal = 0, nsingl = 0, ndupli = 0;
  long st = cputime (), rt = realtime ();
  double est_ideals;

  T = malloc (nthreads * sizeof (tab_t));
  for (i = 0; i < nthreads; i++)
    {
      T[i]->thread = i;
      T[i]->nthreads = nthreads;
      pthread_create (&tid[i], NULL, one_thread3, (void *) (T + i));
    }
  for (i = 0; i < nthreads; i++)
    {
      pthread_join (tid[i], NULL);
      nideal += T[i]->nideal;
      nsingl += T[i]->nsingl;
      ndupli += T[i]->ndupli;
    }
  free (T);
  fprintf (stderr, "stat_mt() took %lds (cpu), %lds (real)\n",
           cputime () - st, realtime () - rt);
  est_ideals = estimate_ideals ((double) M, (double) nideal);
  fprintf (stderr, "Read %lu rels (%lu/s), %lu entries (%.2f%%), est. %1.0f ideals, %lu singletons, %lu of weight 2\n",
           nrels, nrels / (rt - wct0 + (rt == wct0)),
           nideal, 100.0 * (double) nideal / (double) M,
           est_ideals, nsingl, ndupli);
  fprintf (stderr, "Estimated excess %1.0f, need %1.0f\n",
           (double) nrels - est_ideals,
           est_excess ((double) minpr) + est_excess ((double) minpa));
  fflush (stderr);
}

/* read all files and fills the hash table */
static void
doit (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], *ret;
  int i = 0, j;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  long st, rt;

  f = fopen (filelist, "r");
  if (f == NULL)
    {
      fprintf (stderr, "Error, cannot open %s\n", filelist);
      exit (1);
    }
  T = malloc (nthreads * sizeof (tab_t));
  while (!feof (f))
    {
      ret = fgets (g, MAXNAME, f);
      if (ret == NULL)
        break;
      g[strlen(g) - 1] = '\0';
      strcpy (T[i]->g, g);
      T[i]->thread = i;
      i ++;
      if (i == nthreads)
        {
          st = cputime ();
          rt = realtime ();
          for (i = 0; i < nthreads; i++)
            pthread_create (&tid[i], NULL, one_thread, (void *) (T + i));
          for (i = 0; i < nthreads; i++)
            {
              pthread_join (tid[i], NULL);
              nrels += T[i]->nlines;
            }
          fprintf (stderr, "Batch took %lds (cpu), %lds (real)\n",
                   cputime () - st, realtime () - rt);
          i = 0;
          stat_mt (nthreads);
        }
    }
  if (i > 0)
    {
      for (j = 0; j < i; j++)
        pthread_create (&tid[j], NULL, one_thread, (void *) (T + j));
      for (j = 0; j < i; j++)
        {
          pthread_join (tid[j], NULL);
          nrels += T[j]->nlines;
        }
      stat_mt (nthreads);
    }
  fclose (f);
  free (T);
}

/* second pass: outputs remaining relations */
static void
doit2 (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], *ret;
  int i = 0, j;
  tab_t *T;
  pthread_t tid[MAX_THREADS];

  f = fopen (filelist, "r");
  T = malloc (nthreads * sizeof (tab_t));
  while (!feof (f))
    {
      ret = fgets (g, MAXNAME, f);
      if (ret == NULL)
        break;
      g[strlen(g) - 1] = '\0';
      strcpy (T[i]->g, g);
      T[i]->thread = i;
      i ++;
      if (i == nthreads)
        {
          for (i = 0; i < nthreads; i++)
            pthread_create (&tid[i], NULL, one_thread2, (void *) (T + i));
          for (i = 0; i < nthreads; i++)
            pthread_join (tid[i], NULL);
          i = 0;
          fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
                   remains, nrels, 100.0 * (double) remains / (double) nrels);
          fflush (stderr);
        }
    }
  if (i > 0)
    {
      for (j = 0; j < i; j++)
        pthread_create (&tid[j], NULL, one_thread2, (void *) (T + j));
      for (j = 0; j < i; j++)
        pthread_join (tid[j], NULL);
      fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
               remains, nrels, 100.0 * (double) remains / (double) nrels);
      fflush (stderr);
    }
  fclose (f);
  free (T);
}

int
isprime (uint64_t n)
{
  uint64_t p;

  if ((n % 2) == 0)
    return 0;
  for (p = 3; p * p <= n; p += 3)
    if ((n % p) == 0)
      return 0;
  return 1;
}

int
main (int argc, char *argv[])
{
  int nthreads = 1, i;
  char *filelist = NULL;

  /* print command-line arguments */
  fprintf (stderr, "%s", argv[0]);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");
  fflush (stderr);

  fprintf (stderr, "Using PIPE_CAPACITY = %lu\n", PIPE_CAPACITY);
  fflush (stderr);

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-M") == 0)
        {
          M = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-minpa") == 0)
        {
          minpa = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-minpr") == 0)
        {
          minpr = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-t") == 0)
        {
          nthreads = atoi (argv[2]);
          assert (nthreads <= MAX_THREADS);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-filelist") == 0)
        {
          filelist = argv[2];
          argc -= 2;
          argv += 2;
        }
      else if (argv[1][0] == '-')
        {
          fprintf (stderr, "Unknown option: %s\n", argv[1]);
          exit (1);
        }
    }

  if (M == 0)
    {
      fprintf (stderr, "Missing option -M\n");
      exit (1);
    }

  if (minpa == 0)
    {
      fprintf (stderr, "Missing option -minpa\n");
      exit (1);
    }

  if (minpr == 0)
    {
      fprintf (stderr, "Missing option -minpr\n");
      exit (1);
    }

  if (filelist == NULL)
    {
      fprintf (stderr, "Missing option -filelist\n");
      exit (1);
    }

  if (argc != 1)
    {
      fprintf (stderr, "Unknown command-line argument: %s\n", argv[1]);
      exit (1);
    }

  /* check that M = 2p where p is prime */
  assert ((M % 2) == 0);
  assert (isprime (M/2));

  Hsize = 1 + (M - 1) / 32;
  H = (uint64_t*) malloc (Hsize * sizeof (uint64_t));
  memset (H, 0, Hsize * sizeof (uint64_t));

  wct0 = realtime ();

  doit (nthreads, filelist);
  nrels = 0;
  doit2 (nthreads, filelist);

  free (H);

  return 0;
}
