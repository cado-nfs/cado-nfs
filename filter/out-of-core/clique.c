/* This program is an experimental implementation of out-of-core clique
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
   -remove nnn : the number of "cliques" to be removed
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
#include <linux/limits.h>
#include <semaphore.h>
#include <libgen.h>
#include <limits.h> /* for CHAR_BIT */

#define NEW_DIR "new/"
#define MAX_FILES_PER_THREAD 1
#define MAX_STOCK_RAT_PER_LINE 16
#define MAX_STOCK_ALG_PER_LINE 24
#define MAXNAME 1024
#define MAX_THREADS 128
#define PR_TYPE uint32_t
#define MAXLINE 1024

#define LIKELY(X)    __builtin_expect((X), 1)
#define UNLIKELY(X)  __builtin_expect((X), 0)

/* thread structure */
typedef struct
{
  FILE *fin, *fout;
  unsigned long nlines;
  unsigned long nideal;
  unsigned long nsingl;
  unsigned long ndupli;
  char g[ARG_MAX];
  uint8_t busy; /* 0 = free; 1 = work done; 2 = busy */
  uint8_t thread, nthreads;
} __tab_struct;
typedef __tab_struct tab_t[1];

typedef union
{
  double weight;
  int64_t pointer; /* h+1 if points to H2[h] */
} tab2_t;

/* Mix of hashtable, (p,r) array, and weight, to optimize Lx caches :
   200 bytes with LN2BPRHWP = 4.
   So, when (p,r) is found, h has ~80% to be in the same L0 cache line
   and 99.9% to be in the same 16K page (so, less TLB miss and so on).
   CAREFUL: __prhwp_struct must be compact and h type must encode exactly
   BPRHWP entries. The procedures stat and one_thread3 hardcode the
   lenght of h (uint64_t). So don't change LN2BPRHWP.
*/
#define LN2BPRHWP 4
#define BPRHWP (1U<<LN2BPRHWP)
typedef struct
{
  PR_TYPE  pr[BPRHWP]; /* BPRHWP (p,r) of 32 bits */
  uint64_t h;          /* entries of 4 bits : WARNING: hardcoded type! */
  tab2_t   wp[BPRHWP]; /* BPRHWP weight/pointeur of 64 bits */
} __prhwp_struct /* __attribute__ ((__packed__)) */ ;
__prhwp_struct *prhwp;

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "antebuffer 24 %s|gzip -dc", "gzip -c1>%s", },
    { ".bz2", "antebuffer 24 %s|bzip2 -dc", "bzip2 -c1>%s", },
    { ".lzma", "lzma -dc %s", "lzma -c0>%s", },
    /* These two have to be present */
    { "", "antebuffer 24 %s", NULL },
    { NULL, NULL, NULL },
};

#ifndef HAVE_SYNC_FETCH
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* mutual exclusion lock */
#endif
long wct0;

uint64_t M = 0, minpa = 0, minpr = 0, Hsize, nrels, remains;
PR_TYPE *maxpr;
uint64_t target = 0;
double threshold_weight = 999.0;
sem_t sem_pt;


static const unsigned char ugly[256] = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

void
create_directories (char *filelist)
{
  static char odirg[MAXNAME];
  FILE *f;
  char g[MAXNAME], *pdirg, mkdirg[MAXNAME<<1];
  size_t lg;
  int ret;

  f = fopen (filelist, "r");
  while (!feof (f)) {
    if (fgets (g, MAXNAME, f) && (lg = strlen(g))) {
      g[lg-1] = 0;
      pdirg = dirname(g);
      if (strcmp(odirg,pdirg)) {
	strcpy(odirg,pdirg);
	strcpy(mkdirg, "mkdir -p ");
	strcat(mkdirg, NEW_DIR);
	strcat(mkdirg, pdirg);
	ret = system (mkdirg);
        if (ret == -1)
          {
            fprintf (stderr, "An error occurred during system() call\n");
            exit (1);
          }
      }
    }
  }
  fclose(f);
}

double
compute_delay_stat ()
{
  /*
  uint64_t h;
  uint32_t i;
  */
  return (60 + (Hsize>>24));
  /*
  for (h = (Hsize>>24), i = 0; h; h >>=1, i++);
  return (i) ? (double) (i<<6) : 32.;
  */
}

double
realtime ()
{
  struct timeval tv[1];
  /* struct timeval {
     time_t      tv_sec;     seconds
     suseconds_t tv_usec;    microseconds
     }; */

  gettimeofday (tv, NULL);
  return tv->tv_sec + tv->tv_usec / 1000000.0 ;
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

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
    if (strcmp(r->suffix, ne))
      continue;
    if (r->pfmt_out)
      pclose(f);
    else
      fclose(f);
    return;
  }
}

static inline uint64_t
index_hash (uint64_t pr)
{
  PR_TYPE *prh, *prm, mpr;

  mpr = (PR_TYPE) pr;
  pr %= M;
  prm = prhwp[pr >> LN2BPRHWP].pr;
  prh = prm + (pr & (BPRHWP-1));
  
  /* Fastest & most common procedure: hash[pr] == pr or hash[pr] == 0 */
  if (*prh == mpr)
    return pr;
  if (!(*prh)) {
    *prh = mpr;
    return pr;
  }
  
  prh++;
  pr++;
  prm += BPRHWP;
  if (prm > maxpr) prm = maxpr;
  for (;;) {
    if (UNLIKELY(prh == prm)) {
      if (LIKELY(prm != maxpr)) {
	prm = (PR_TYPE *) ((void *) prm + sizeof(*prhwp));
	prh = prm - BPRHWP;
      } else {
	pr = 0;
	prh = prhwp[0].pr;
	prm = prh + BPRHWP;
      }
    }
    if (*prh == mpr)
      return pr;
    if (!(*prh)) {
      *prh = mpr;
      return pr;
    }
    prh++;
    pr++;
  }
  return pr;
}

#define INDEX_RAT(P) (index_hash ((P) + 1)) /* even for even M */
static inline uint64_t
index_rat (uint64_t p)
{
  return INDEX_RAT(p);
}

#define INDEX_ALG(P,R) (index_hash ((P) + (M - 2) * (R)))  /* always odd for odd p and even M */
static inline uint64_t
index_alg (uint64_t p, uint64_t r)
{
  return INDEX_ALG(p,r);
}

#define WEIGHT(P) ((prhwp[(P) >> LN2BPRHWP].h >> (((P)&(BPRHWP-1)) << 2)) & (BPRHWP-1))
static inline uint8_t
weight (uint64_t h)
{
  return (uint8_t) WEIGHT(h);
}

static void
insert (uint64_t h)
{
  uint8_t *hh;
  hh = ((uint8_t *) &(prhwp[h >> LN2BPRHWP].h)) + ((h & (BPRHWP-1)) >> 1);
  if (h & 1) {
    if ((*hh & 0xf0) != 0xf0) *hh += 0x10;
  } else
    if ((*hh & 0x0f) != 0x0f) *hh += 0x01;
}

static void
delete (uint64_t h)
{
  uint8_t *hh, v;
  static int count = 0;

  hh = ((uint8_t *) &(prhwp[h >> LN2BPRHWP].h)) + ((h & (BPRHWP-1)) >> 1);
  if (h & 1) {
    v = *hh & 0xf0;
    if (v) {
      if (v != 0xf0) *hh -= 0x10;
    } else
      if (count++ < 10) fprintf (stderr, "Warning: deleting weight-0 ideal 0 at h=%lu\n", h);
  } else {
    v = *hh & 0x0f;
    if (v) {
      if (v != 0x0f) *hh -= 0x01;
    } else
      if (count++ < 10) fprintf (stderr, "Warning: deleting weight-0 ideal 0 at h=%lu\n", h);
  }
}

static double
clique_weight (uint64_t h)
{
  tab2_t *pw;
  
  for (;;) {
    pw = prhwp[h >> LN2BPRHWP].wp + (h & (BPRHWP-1));
    if (pw->pointer <= 0) return -(pw->weight);
    h = pw->pointer - 1;
  }
}

static double
propagate_weight (uint64_t h)
{
  tab2_t *pw;

  pw = prhwp[h >> LN2BPRHWP].wp + (h & (BPRHWP-1));
  if (pw->pointer > 0)
    pw->weight = propagate_weight (pw->pointer - 1);
  return pw->weight;
}

/* adds the weight w to the connected component of H2[h] */
static uint64_t
insert2 (uint64_t h, double w)
{
  tab2_t *pw;

  for (;;) {
    pw = prhwp[h >> LN2BPRHWP].wp + (h & (BPRHWP-1));
    if (pw->pointer <= 0) break;
    h = pw->pointer - 1;
  }
  pw->weight -= w;
  return h;
}

/* connects the connected component H2[h] to H2[h0],
   where H2[h0] stores a weight */
static void
insert2a (uint64_t h, uint64_t h0)
{
  double *w0;
  tab2_t *pw;

  for (;;) {
    pw = prhwp[h >> LN2BPRHWP].wp + (h & (BPRHWP-1));
    if (pw->pointer <= 0) break;
    h = pw->pointer - 1;
  }
  if (h == h0) /* already connected */
    return;
  w0 = &(prhwp[h0 >> LN2BPRHWP].wp[h0 & (BPRHWP-1)].weight);
  if (*w0 > 0 || pw->weight > 0)
    return; /* avoid an assert which might fail */
  *w0 += pw->weight;
  pw->pointer = h0 + 1;
}

#define SAMPLE 1000       /* number of intervals for component weight */
#define RESOLUTION 100.0  /* 1/width of each interval */

static void
stat2 ()
{
  uint64_t h, nc = 0; /* number of connected components */
  double wmax = 0.0, wsum = 0.0, w;
  uint64_t count[SAMPLE+1], n;

  memset(count, 0, sizeof(*count) * (SAMPLE + 1));
  for (h = 0; h < M; h++)
    if ((w = prhwp[h >> LN2BPRHWP].wp[h & (BPRHWP-1)].weight) < 0)
      {
	w = -w;
        nc++;
        wsum += w;
        if (w > wmax)
          wmax = w;
        n = (uint64_t) (RESOLUTION * w);
        if (n >= SAMPLE)
          n = SAMPLE - 1;
        count[n]++;
      }
  fprintf (stderr,
	   "Number of connected components: %lu\n"
	   "Weight: %.2f (avg), %.2f (max)\n", nc, wsum / (double) nc, wmax);
  n = SAMPLE;
  while (n > 0)
    {
      count[n-1] += count[n];
      if (count[n] > 0)
        fprintf (stderr, "Weight >= %.2f: %lu components\n",
                 (double) (n - 1) / RESOLUTION, count[n-1]);
      if (count[n-1] > target)
        break;
      n--;
    }
  threshold_weight = (double) n / RESOLUTION;
  fprintf (stderr, "Will remove %lu components of weight >= %.2f\n",
           count[n], threshold_weight);
}

/* return a*b mod p, assuming p < 2^40 */
static inline uint64_t
mulmod (uint64_t a, uint64_t b, uint64_t p)
{
  uint64_t r;

  a %= p; /* now 0 <= a < 2^40 */
  b %= p; /* now 0 <= b < 2^40 */

  /* let b = b1 * 2^20 + b0: we compute (b1*a)*2^20 + b0*a */
  r = ((b >> 20) * a) % p; /* b1*a < 2^60 and r < 2^40 */
  return ((r << 20) + (b & 0xfffff) * a) % p; /* r*2^20 < 2^60, b0*a < 2^60 */
}

/* return a/b mod p */
static inline unsigned long
root (long a0, long b0, long p)
{
  long u, w, q, r, t, a, b;

  b = b0 % p;
  assert (b != 0);
  a = b;
  u = 1;
  w = 0;
  b = p;
  /* invariant: a = u*b0 mod p, b = w*b0 mod p */
  while (b) {
      q = a / b;
      r = a - q * b; /* r = a % b; */
      a = b;
      b = r;
      t = u;
      u = w;
      w = (t - q * w) % p;
  }
  /* compute a0*u % p */
  if (u < 0) u += p;
  return (a0 >= 0) ? mulmod (a0, u, p) : p - mulmod (-a0, u, p);
}

/* first pass, return the number of lines read */
void*
pass1_one_thread (void* args)
{
  char s[MAXLINE], *t;
  tab_t *tab = (tab_t*) args;
  FILE *f;
  uint64_t p, p2;
  unsigned long line, b, r;
  long a;
  uint8_t m;
  
  f = (tab[0]->fin) ? tab[0]->fin : gzip_open(tab[0]->g, "r");
  assert (f);
  line = 0;
  while (fgets (s, MAXLINE, f)) {
    line++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = p * 10 + m);
    a *= p;
    assert (a != 0);
    assert (*t == ',');
    t++;
    for (b = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, b = b * 10 + m);
    assert (*t == ':');
    do {
      t++;
      for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = (p << 4) + m);
      if (p >= minpr) {
	insert (INDEX_RAT (p));
      }
    } while (*t == ',');
    if (b) {
      do {
	t++;
	for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = (p << 4) + m);
	if (p >= minpa) {
	  r = root (a, b, p);
	  insert (INDEX_ALG (p, r));
	}
      } while (*t == ',');
    } 
    else
      if (p >= minpa)
	do {
	  t++;
	  for (p2 = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p2 = (p2 << 4) + m);
	  insert (INDEX_ALG (p, p2));
	} while (*t == ',');
  }
  tab[0]->nlines = line;
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  gzip_close (f, tab[0]->g);
  return NULL;
}

/* pass 2: insert in H2 ideals with weight 2 */
void*
pass2_one_thread (void *args)
{
  uint64_t rp[MAX_STOCK_RAT_PER_LINE], ap[MAX_STOCK_ALG_PER_LINE], ar[MAX_STOCK_ALG_PER_LINE];
  char s[MAXLINE], *t;
  tab_t *tab = (tab_t*) args;
  FILE *f;
  uint64_t p, p2, i, h;
  double W;
  unsigned long line, output, b, r;
  long a;
  unsigned int nr, na;
  uint8_t m, w;

  f = (tab[0]->fin) ? tab[0]->fin : gzip_open(tab[0]->g, "r");
  line = output = 0;
 next_line:
  while (fgets (s, MAXLINE, f)) {
    line++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = p * 10 + m);
    a *= p;
    t++;
    for (b = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, b = b * 10 + m);
    nr = na = 0;
    W = 0.0;
    do {
      t++;
      for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
      if (p >= minpr) {
	i = INDEX_RAT(p);
	w = WEIGHT (i);
	if (w < 2) goto next_line;
	if (w == 2)
	  rp[nr++] = p;
	else
	  W += 1.0 / (double) w;
      }
    } while (*t == ',');
    if (b)
      do {
	t++;
	for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
	if (p >= minpa) {
	  r = root (a, b, p);
	  i = INDEX_ALG (p, r);
	  w = WEIGHT (i);
	  if (w < 2) goto next_line;
	  if (w == 2) {
	    ap[na] = p;
	    ar[na++] = r;
	  } else
	    W += 1.0 / (double) w;
	}
      } while (*t == ',');
    else
      if (p >= minpa)
	do {
	  t++;
	  for (p2 = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p2 = (p2 << 4) + m);
	  i = INDEX_ALG (p, p2);
	  w = WEIGHT (i);
	  if (w < 2) goto next_line;
	  if (w == 2) {
	    ap[na] = p;
	    ar[na++] = p2;
	  } else
	    W += 1.0 / (double) w;
	} while (*t == ',');
    output++;
    /* Warning: it can be that an ideal appears twice (or more).
       To avoid reducing exponents mod 2, which would be expensive,
       we insert this ideal twice, but then we should avoid for infinite
       loops. */
    if (nr-- > 0) {
      h = insert2 (INDEX_RAT (rp[nr]), W);
      while (nr-- > 0) insert2a (INDEX_RAT (rp[nr]), h);
      while (na-- > 0) insert2a (INDEX_ALG (ap[na], ar[na]), h);
    }
    else if (na-- > 0) {
      h = insert2 (INDEX_ALG (ap[na], ar[na]), W);
      while (na-- > 0) insert2a (INDEX_ALG (ap[na], ar[na]), h);
    }
  }
  gzip_close (f, tab[0]->g);
#ifdef HAVE_SYNC_FETCH
  __sync_fetch_and_add(&remains, output);
  __sync_fetch_and_add(&nrels, line);
#else
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
#endif
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  return NULL;
}

/* pass 3: writes remaining relations */
void*
pass3_one_thread (void* args)
{
  uint64_t rp[MAX_STOCK_RAT_PER_LINE], ap[MAX_STOCK_ALG_PER_LINE], ar[MAX_STOCK_ALG_PER_LINE];
  char s[MAXLINE], newg[MAXNAME<<1], *t;
  tab_t *tab = (tab_t*) args;
  FILE *f, *newf;
  uint64_t p, p2, w, h;
  unsigned long line = 0, b, r, output = 0;
  long a;
  unsigned int nr, na;
  uint8_t m;

  if ((f = tab[0]->fin))
     newf = tab[0]->fout;
  else {
    f = gzip_open(tab[0]->g, "r");
    strcpy (newg, NEW_DIR);
    strcat (newg, tab[0]->g);
    newf = gzip_open (newg, "w");
  }
  assert (newf);
  line = output = 0;
  if (0) {
  next_line: /* update 2-bit counters of removed ideals of the suppress line */
    while (nr-- > 0)
      delete (INDEX_RAT (rp[nr]));
    while (na-- > 0)
      delete (INDEX_ALG (ap[na], ar[na]));
  }
  while (fgets (s, MAXLINE, f)) {
    line++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = p * 10 + m);
    a *= p;
    t++;
    for (b = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, b = b * 10 + m);
    nr = na = 0;
    do {
      t++;
      for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
      if (p >= minpr) {
	rp[nr++] = p;
	h = INDEX_RAT (p);
	w = WEIGHT (h);
	if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	  goto next_line;
      }
    } while (*t == ',');
    if (b)
      do {
	t++;
	for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
	if (p >= minpa) {
	  r = root (a, b, p);
	  ap[na] = p;
	  ar[na++] = r;
	  h = INDEX_ALG (p, r);
	  w = WEIGHT (h);
	  if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	    goto next_line;
	}
      } while (*t == ',');
    else
      if (p >= minpa)
	do {
	  t++;
	  for (p2 = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p2 = (p2 << 4) + m);
	  ap[na] = p;
	  ar[na++] = p2;
	  h = INDEX_ALG (p, p2);
	  w = WEIGHT (h);
	  if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	    goto next_line;
	} while (*t == ',');
    fputs (s, newf);
    output++;
  }
  gzip_close (f, tab[0]->g);
  gzip_close (newf, tab[0]->g);
  fprintf (stderr, "   new/%s done: remains %lu rels out of %lu\n",
           tab[0]->g, output, line);
#ifdef HAVE_SYNC_FETCH
  __sync_fetch_and_add(&remains, output);
  __sync_fetch_and_add(&nrels, line);
#else
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
#endif
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  return NULL;
}

void*
one_thread3 (void* args)
{
  tab_t *tab = (tab_t*) args;
  uint64_t *hd, *he, h, wt[16];

  memset (wt, 0, sizeof(*wt) * 16);
  hd = &(prhwp[(tab[0]->thread * Hsize) / tab[0]->nthreads].h);
  he = &(prhwp[((tab[0]->thread + 1) * Hsize) / tab[0]->nthreads].h);
  while (hd < he)
    {  /* deal with 64 bits, i.e., 16 entries at a time */
      h = *hd;
      hd = (uint64_t *) ((void *) hd + sizeof(*prhwp));
      __builtin_prefetch(hd);
      while (h)
        {
          wt[h & 0xf] ++;
          h >>= 4;
        }
    }
  tab[0]->nideal = tab[0]->nsingl = wt[1];
  tab[0]->ndupli = wt[2];
  for (hd = &(wt[2]); hd < &(wt[16]); tab[0]->nideal += *hd++);
  return NULL;
}

/* given a hash table of n entries, of which ek entries are already filled,
   return an estimate of the number k of insertions in the table.
   We assume that e(k+1) = e(k) + 1 - e(k)/n with e(0) = 1,
   which gives e(k) = n - n*(1-1/n)^k,
   thus k = log(1-e(k)/n)/log(1-1/n) */
static double
estimate_ideals (double n, double ek)
{
  return ek;
}

/* estimate the number of prime ideals below a */
static double
est_excess (double a)
{
  return a / log (a);
}

/* mono-thread version */
static void*
stat () {
  uint64_t *hd, *he, h, wt[16];
  long st = seconds (), rt = realtime ();
  unsigned long nideal, nsingl, ndupli;
  double est_ideals;
  int dum;

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,&dum);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,&dum);
  /* pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,&dum); */
  memset (wt, 0, sizeof(*wt) * 16);
  for (hd = &(prhwp[0].h), he = &(prhwp[Hsize].h); hd < he; )
    {  /* deal with 64 bits, i.e., 16 entries at a time */
      h = *hd;
      hd = (uint64_t *) ((void *) hd + sizeof(*prhwp));
      __builtin_prefetch(hd);
      while (h)
        {
          wt[h & 0xf] ++;
          h >>= 4;
        }
    }
  nideal = nsingl = wt[1];
  ndupli = wt[2];
  for (hd = &(wt[2]); hd < &(wt[16]); nideal += *hd++);
  fprintf (stderr, "   Stat() took %lds (cpu), %lds (real)\n",
           (unsigned long) (seconds () - st), (unsigned long) (realtime () - rt));
  est_ideals = estimate_ideals ((double) M, (double) nideal);
  fprintf (stderr, "   Read at least %lu rels (%lu/s), %lu ideals (%.2f%%), %lu singletons, %lu of weight 2\n",
           nrels, (unsigned long) (nrels / (rt - wct0)),
           nideal, (100.0 * nideal) / M, nsingl, ndupli);
  fprintf (stderr, "   Estimated excess %1.0f, need %1.0f\n",
           (double) nrels - est_ideals,
           est_excess ((double) minpr) + est_excess ((double) minpa));
  return NULL;
}

/* multi-thread version */
static void
stat_mt (int nthreads)
{
  int i;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  unsigned long nideal = 0, nsingl = 0, ndupli = 0;
  long st = seconds (), rt = realtime ();
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
  fprintf (stderr, "   Stat_mt() took %lds (cpu), %lds (real)\n",
           (unsigned long) (seconds () - st), (unsigned long) (realtime () - rt));
  est_ideals = estimate_ideals ((double) M, (double) nideal);
  fprintf (stderr, "   Read %lu rels (%lu/s), %lu ideals (%.2f%%), %lu singletons, %lu of weight 2\n",
           nrels, (unsigned long) (nrels / (rt - wct0)),
           nideal, (100.0 * nideal) / M, nsingl, ndupli);
  fprintf (stderr, "   Estimated excess %1.0f, need %1.0f\n",
           (double) nrels - est_ideals,
           est_excess ((double) minpr) + est_excess ((double) minpa));
}

/* read all files and fills the hash table */
static void
pass1 (int nthreads, char *filelist)
{
  FILE *f;
  char g[ARG_MAX], *pg;
  int i, j = 0, nbt = 0, notload, nbf, notfirst = 0, avoidwarning = 0;
  tab_t *T;
  pthread_t tid[MAX_THREADS+1], sid;
  double st = seconds(), rt = realtime(), crt, art,
    delay_stat = compute_delay_stat();
  size_t lpg;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 1/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  if (f == NULL) {
    fprintf (stderr, "Error, cannot open %s\n", filelist);
    exit (1);
  }
  T = malloc ((nthreads + 1) * sizeof (tab_t));
  memset(T, 0, (nthreads + 1) * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads);
  crt = art = rt;
  while ((notload = (!feof (f))) || nbt) {
    if (notload) {
      pg = g;
      nbf = 0;
      do {
	if (fgets (pg, MAXNAME, f) && (lpg = strlen(pg))) {
	  nbf++;
	  pg += lpg;
	  pg[-1] = ' ';
	}
	else break;
      } while (nbf < MAX_FILES_PER_THREAD && pg < g + ARG_MAX - MAXNAME);
      if (nbf) {
	pg[-1] = 0;
	for (i = 0; T[i]->busy && i <= nthreads; i++);
	assert(i <= nthreads);
	strcpy (T[i]->g, g);
	T[i]->busy = 2;
	T[i]->thread = i;
	if (nbt++ < nthreads) {
	  T[i]->fin = NULL;
	  pthread_create (&tid[i], NULL, pass1_one_thread, (void *) (T + i));
	} else {
	  avoidwarning = i;
	  T[i]->fin = gzip_open(g, "r");
	}
      }
    }
    sem_wait(&sem_pt);
    if (nbt > nthreads)
      pthread_create (&tid[avoidwarning], NULL, pass1_one_thread, (void *) (T + avoidwarning));
    for (i = 0; i <= nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	nrels += T[i]->nlines;
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j == nthreads) {
      fprintf (stderr, "Pass1; rels: load %lu; krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (seconds() - st),
	       (unsigned long) (art - rt));
      j = 0;
    }
    if (art - crt > delay_stat) {
      if (notfirst) {
	pthread_cancel(sid);
	pthread_join (sid, NULL);
      } else
	notfirst = 1;
      pthread_create (&sid, NULL, stat, NULL);
      crt = art;
    }
  }
  if (notfirst) {
    pthread_cancel(sid);
    pthread_join (sid, NULL);
  }
  fprintf (stderr, "\n*** Pass 1, final stats: *** \n");
  stat_mt (nthreads);
  fprintf (stderr, "Pass 1 took %lds (cpu), %lds (real)\n",
           (unsigned long) (seconds () - st), (unsigned long) (realtime () - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 1/3 DONE *\n"
	   "*****************\n\n");
}

/* second pass: stores ideals of weight 2 */
static void
pass2 (int nthreads, char *filelist)
{
  FILE *f;
  char g[ARG_MAX], *pg;
  tab_t *T;
  pthread_t tid[MAX_THREADS+1];
  double st = seconds(), rt = realtime(), art;
  size_t lpg;
  uint64_t h;
  int i, j = 0, notload, nbt = 0, nbf, avoidwarning = 0;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 2/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  T = malloc ((nthreads + 1) * sizeof (tab_t));
  memset(T, 0, (nthreads + 1) * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads);
  art = rt;
  while ((notload = (!feof (f))) || nbt) {
    if (notload) {
      pg = g;
      nbf = 0;
      do {
	if (fgets (pg, MAXNAME, f) && (lpg = strlen(pg))) {
	  nbf++;
	  pg += lpg;
	  pg[-1] = ' ';
	}
	else
	  break;
      } while (nbf < MAX_FILES_PER_THREAD && pg < g + ARG_MAX - MAXNAME);
      if (nbf) {
	pg[-1] = 0;
	for (i = 0; T[i]->busy && i <= nthreads; i++);
	assert(i <= nthreads);
	strcpy (T[i]->g, g);
	T[i]->busy = 2;
	T[i]->thread = i;
	if (nbt++ < nthreads) {
	  T[i]->fin = NULL;
	  pthread_create (&tid[i], NULL, pass2_one_thread, (void *) (T + i));
	} else {
	  T[i]->fin = gzip_open(g, "r");
	  avoidwarning = i;
	}
      }
    }
    sem_wait(&sem_pt);
    if (nbt > nthreads)
      pthread_create (&tid[avoidwarning], NULL, pass2_one_thread, (void *) (T + avoidwarning));
    for (i = 0; i <= nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j == nthreads) {
      fprintf (stderr, "Pass2; rels: load %lu, used %lu (%.2f%%); "
	       "krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels, remains,
	       (100.0 * remains) / nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (seconds() - st),
	       (unsigned long) (art - rt));
      j = 0;
    }
  }
  art = realtime();
  fprintf (stderr, "\n*** Pass 2, final stats: *** \n"
	   "Relations: load %lu, used %lu (%.2f%%); "
	   "krels/s: %lu; time: %lus cpu, %lus real\n",
	   nrels, remains,
	   (100.0 * remains) / nrels,
	   (unsigned long) (nrels / ((art - rt) * 1000)),
	   (unsigned long) (seconds() - st),
	   (unsigned long) (art - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  stat2 ();
  /* propagate weights to speed up clique_weight */
  for (h = 0; h < M; h++)
    propagate_weight (h);
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 2/3 DONE *\n"
	   "*****************\n\n");
}

/* third pass: outputs remaining relations */
/* Problem: I cannot agglomerate the files because 1 file load = 1 file written.
   So I have to take them one by one and write them synchronously, one by one too.
   It's a pity for disk access...
*/
static void
pass3 (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], newg[MAXNAME<<1];
  int i, j = 0, k = 0, nbt = 0, notload, avoidwarning = 0;
  tab_t *T;
  pthread_t tid[MAX_THREADS+1];
  double st = seconds(), rt = realtime(), art;
  size_t lg;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 3/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  T = malloc ((nthreads + 1) * sizeof (tab_t));
  memset(T, 0, (nthreads + 1) * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads);
  art = rt;
  while ((notload = (!feof (f))) || nbt) {
    if (notload && fgets (g, MAXNAME, f) && (lg = strlen(g))) {
      g[lg - 1] = 0;
      for (i = 0; T[i]->busy && i <= nthreads; i++);
      assert(i <= nthreads);
      strcpy (T[i]->g, g);
      T[i]->busy = 2;
      T[i]->thread = i;
      if (nbt++ < nthreads) {
	T[i]->fin = NULL;
	pthread_create (&tid[i], NULL, pass3_one_thread, (void *) (T + i));
      }
      else {
	T[i]->fin = gzip_open(g, "r");
	strcpy (newg, NEW_DIR);
	strcat (newg, T[i]->g);
	T[i]->fout = gzip_open (newg, "w");
	avoidwarning = i;
      }
    }
    sem_wait(&sem_pt);
    if (nbt > nthreads)
      pthread_create (&tid[avoidwarning], NULL, pass3_one_thread, (void *) (T + avoidwarning));
    for (i = 0; i <= nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j == nthreads && ++k == MAX_FILES_PER_THREAD) {
      fprintf (stderr, "Pass3; rels: load %lu, used %lu (%.2f%%); "
	       "krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels, remains,
	       (100.0 * remains) / nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (seconds() - st),
	       (unsigned long) (art - rt));
      j = k = 0;
    }
  }
  art = realtime();
  fprintf (stderr, "\n*** Pass 3, final stats: *** \n"
	   "Relations: load %lu, used %lu (%.2f%%); "
	   "krels/s: %lu; time: %lus cpu, %lus real\n",
	   nrels, remains,
	   (100.0 * remains) / nrels,
	   (unsigned long) (nrels / ((art - rt) * 1000)),
	   (unsigned long) (seconds() - st),
	   (unsigned long) (art - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 3/3 DONE *\n"
	   "*****************\n\n");
}

int
isprime (uint64_t n)
{
  uint64_t p;

  if (!(n & 1)) return 0;
  for (p = 3; p * p <= n; p += 3)
    if (!(n % p))
      return 0;
  return 1;
}

int
main (int argc, char *argv[])
{
  int nthreads = 1, i;
  char *filelist = NULL;
  size_t size_malloc;

  /* print command-line arguments */
  fprintf (stderr, "%s", argv[0]);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

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
      else if (argc > 2 && strcmp (argv[1], "-remove") == 0)
        {
          target = atoi (argv[2]);
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

  if (target == 0)
    {
      fprintf (stderr, "Missing option -remove\n");
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
  assert (!(M & 1));
  assert (isprime (M>>1));
  /* an *hpr stocks 16 entries. */
  Hsize = (M + BPRHWP-1) >> LN2BPRHWP; /* each entry uses 4 bits AND h is uint64_t */
  size_malloc = Hsize * sizeof (*prhwp);
  fprintf (stderr, "Creating & clearing hashtable + (p,r) array + weight array : %lu mB...", size_malloc >> 20);
  if (posix_memalign ((void **) &prhwp, 1<<14, size_malloc)) {
    fprintf (stderr, "Malloc error.\n");
    exit (1);
  }
  memset (prhwp, 0, size_malloc);
  fprintf (stderr, " Done. \n");
  fprintf (stderr, "Creating possible " NEW_DIR "* directories...");
  create_directories(filelist);
  fprintf (stderr, " Done.\n");

  fprintf (stderr, "Each thread processing %d file(s)\n",
           MAX_FILES_PER_THREAD);

  maxpr = prhwp[M>>LN2BPRHWP].pr + (M&(BPRHWP-1));
  wct0 = realtime ();
  nrels = 0;
  pass1 (nthreads, filelist);

  nrels = remains = 0;
  pass2 (nthreads, filelist);

  nrels = remains = 0;
  pass3 (nthreads, filelist);

  free (prhwp);

  return 0;
}
