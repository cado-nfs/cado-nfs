#ifndef NBUCKETK
#define NBUCKETK 7
#endif
#define NBUCKET (1<<NBUCKETK)
#define CLSIZE ((size_t)64)
#define ENTRIES_PER_LINE (CLSIZE / sizeof(unsigned int))
#define BUCKET_SIZE ((size_t)1<<(27-NBUCKETK))
/* Pseudo-random sequence */
#define RAND(x) (5*(x)+1)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <emmintrin.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include "portability.h"
              
/*
  http://www.mjmwired.net/kernel/Documentation/vm/pagemap.txt
*/


int
check_cache_friendliness(const void * const alloc, const size_t size, 
                         const size_t pagesize, 
                         const size_t cache_size,
                         const unsigned int associativity)
{
  const size_t npages = (size - 1) / pagesize + 1;
  const size_t associative_set_pages = cache_size / associativity / pagesize;
  uint64_t *memmap;
  FILE *f;
  const int verbose = 0;
  unsigned int *tlbsethits;
  size_t i;
  pid_t pid;
  char fn[128];

  pid = getpid();
  snprintf (fn, sizeof(fn), "/proc/%lu/pagemap", (unsigned long) pid);
  f = fopen(fn, "r");
  if (f == NULL) {
    /* Cannot check, return ok */
    fprintf (stderr, "Could not %s for reading\n", fn);
    return 1;
  }
  
  memmap = malloc(8 * npages);
  fseek(f, (uintptr_t) alloc / pagesize * 8, SEEK_SET);
  i = fread(memmap, 8, npages, f);
  fclose(f);
  if (i != npages){
    fprintf (stderr, "Could not read %zu qwords from %s\n", npages, fn);
    return 1;
  }
  
  tlbsethits = malloc(associative_set_pages * sizeof(unsigned int));
  memset (tlbsethits, 0, associative_set_pages * sizeof(unsigned int));
  
  for (i = 0; i < npages; i++) {
    /* The low 54 bits are the page frame number, i.e., the physical
       memory address is PFN * pagesize */
    uint64_t PFN = memmap[i] % ((uint64_t)1 << 55);
    if (verbose) {
      printf ("Page %zu has PFN %lu\n", i, PFN);
      fflush(stdout);
    }
    tlbsethits[PFN % associative_set_pages]++;
  }
  
  free(memmap);
  
  for (i = 0; i < associative_set_pages; i++) {
    if (tlbsethits[i] > associativity) {
      if (verbose) {
        printf ("TLB set %zu: %u hits\n", i, tlbsethits[i]);
        fflush(stdout);
      }
      break;
    }
  }
  free(tlbsethits);
  if (i < associative_set_pages) {
    if (verbose) {
      printf ("Memory not cache friendly\n");
      fflush(stdout);
    }
    return 0;
  } else {
    if (verbose) {
      printf ("Memory cache friendly\n");
      fflush(stdout);
    }
    return 1;
  }
}

void * alloc_friendly(const size_t size, const size_t pagesize)
{
  const int verbose = 1;
  void *alloc;
  const size_t cachesize = 65536;
  const int associativity = 2;
  unsigned int tries = 0;

  do {
    int i;
    tries ++;
    i = posix_memalign (&alloc, pagesize, size);
    if (i != 0)
      exit(EXIT_FAILURE);
    /* Make sure pages are mapped */
    memset (alloc, 0, size);
    if (1 && /* disabled, makes no difference on Core 2 with 8-way associative cache */
        size <= cachesize && 
        !check_cache_friendliness(alloc, size, pagesize, cachesize, associativity)) {
      // free (cacheline_alloc);
      alloc = NULL;
    }
  } while (alloc == NULL);
  if (verbose) {
    printf ("Getting cache-friendly pages took %u tries\n", tries);
    fflush(stdout);
  }
  return alloc;
}


/* Since BUCKET_SIZE is a large power of 2, the start of the buckets
   will be a power-of-two size apart in memory, possibly causing a lot
   of cache line collisions. Thus we leave gaps of size bucket_skip
   between the end of a bucket and the start of the next, to make the 
   start address of the buckets not all map to the same few cache lines.
   This did not seem to make much difference in speed.
   Furthermore, the TLB are also set associative, and a large power-of-2
   distance between write pointers causes very frequent mapping of the 
   pages to which the bucket write pointers point to the same TLB set.
   Leaving a page-sized gap fixes this, and makes a noticable difference 
   in speed. */
static const size_t bucket_skip = 4096;
static const size_t fill_iter = BUCKET_SIZE * NBUCKET * 7 / 8; /* 2^27 * 7/8 */


/* Fill buckets, without cacheline buffer. Simply one write to the bucket 
   memory area per update */
void fill(unsigned int **buckets, unsigned int start) {
  unsigned int loc = start;
  size_t i;
  
  for (i = 0; i < fill_iter; i++) {
    const unsigned int bucket = loc >> (32 - NBUCKETK);
    *buckets[bucket]++ = (unsigned int) loc;
    loc = RAND(loc);
  }
}


void fill_cacheline(unsigned int ** restrict buckets, unsigned int start, 
                    void * restrict const cacheline_alloc) {
  unsigned char * restrict const nr = cacheline_alloc;
  unsigned int * restrict const cacheline = (unsigned int *)((unsigned char *) cacheline_alloc + NBUCKET);
  unsigned int loc = start;
  size_t i;

  // unsigned int cacheline[NBUCKET][ENTRIES_PER_LINE] __attribute__((__aligned__(64)));
  
  /* Putting % 128 here makes the code fail half of the time, which is expected,
     but in the cases where it works, it is about 5% faster, which is not expected.
     Oddly, putting %128 here and __aligned__(128) above does NOT produce the speedup */
  // assert ((uintptr_t)&cacheline % 64 == 0);

  // printf ("cacheline = %p\n", cacheline);
  
  for (i = 0; i < NBUCKET; i++)
    nr[i] = 0;
  
  for (i = 0; i < fill_iter; i++) {
    const unsigned int bucket = (loc >> (32 - NBUCKETK));
    cacheline[bucket * ENTRIES_PER_LINE + (nr[bucket]++)] = (unsigned int) loc;
    
    if (nr[bucket] == ENTRIES_PER_LINE) {
      __m128i *cl = (__m128i *) &cacheline[bucket * ENTRIES_PER_LINE];
      __m128i *b = (__m128i *) buckets[bucket];
      buckets[bucket] += ENTRIES_PER_LINE;
      _mm_stream_si128 (&(b[0]), cl[0]);
      if (CLSIZE > 16)
        _mm_stream_si128 (&(b[1]), cl[1]);
      if (CLSIZE > 32)
        _mm_stream_si128 (&(b[2]), cl[2]);
      if (CLSIZE > 48)
        _mm_stream_si128 (&(b[3]), cl[3]);
      if (CLSIZE > 64)
        _mm_stream_si128 (&(b[4]), cl[4]);
      if (CLSIZE > 80)
        _mm_stream_si128 (&(b[5]), cl[5]);
      if (CLSIZE > 96)
        _mm_stream_si128 (&(b[6]), cl[6]);
      if (CLSIZE > 112)
        _mm_stream_si128 (&(b[7]), cl[7]);
      nr[bucket] = 0;
    }
    loc = RAND(loc);
  }
}

unsigned long 
fill_cache(void * const alloc, const size_t size, const unsigned int iter)
{
  unsigned long i, k;
  volatile unsigned long * const arr = alloc;
  unsigned long sum = 0;

  if (size / 8 == 0 || iter == 0)
    return 0;
  
  for (i = 0; i < 256 * iter; i++) {
    for (k = 0; k < size / 64 * KMULT / 512; k++) {
      sum += arr[k * 8];
    }
  }
  printf ("%zu reads total, sum = %lu\n", 256 * iter * size / 64, sum);
  return sum;
}

double time()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (ru.ru_utime.tv_sec * 1000000 + ru.ru_utime.tv_usec) / 1.e6;
}

unsigned int * 
bucket_start(unsigned int * const bucket_storage, size_t i)
{
  return bucket_storage + i * (BUCKET_SIZE + bucket_skip);
}

void printpointers(unsigned int ** const buckets, unsigned int * const bucket_storage)
{
  int i;

  printf("Bucket pointers:\n");
  for (i = 0; i < NBUCKET; i++) {
    char *start = (char *) bucket_start(bucket_storage, i);
    char *end = (char *) (buckets[i]);
    printf ("%p - %p, len = %zu\n", start, end, end - start);
    assert ((size_t)(end-start) < BUCKET_SIZE * sizeof(unsigned int));
  }
  printf("\n");
}


int main(int argc, char **argv) {
  size_t pagesize;
  int iter = 10;
  unsigned int **buckets;
  unsigned int *bucket_storage;
  void *cacheline_alloc;
  unsigned int *cacheline;
  unsigned char *nr;
  int i, j, which, do_which[2] = {1, 1}, do_cache = 0, friendly = 0;
  double starttime, endtime;
  const char *which_str[2] = {"Without", "With"};
  
  while (argc > 1) {
    if (strcmp (argv[1], "-cl") == 0) {
      do_which[0] = 0;
      argc--;
      argv++;
      continue;
    }
    
    if (strcmp (argv[1], "-ncl") == 0) {
      do_which[1] = 0;
      argc--;
      argv++;
      continue;
    }

    if (strcmp (argv[1], "-cache") == 0) {
      do_which[0] = do_which[1] = 0;
      do_cache = 1;
      argc--;
      argv++;
      continue;
    }

    if (strcmp (argv[1], "-friendly") == 0) {
      friendly = 1;
      argc--;
      argv++;
      continue;
    }

    if (argc > 2 && strcmp (argv[1], "-i") == 0) {
      iter = atoi (argv[2]);
      argc -= 2;
      argv += 2;
      continue;
    }
    
    fprintf (stderr, "Unknown argument: %s\n", argv[1]);
    exit (EXIT_FAILURE);
  }

  pagesize = sysconf(_SC_PAGESIZE);

  if (do_cache) {
    size_t alloc_size = NBUCKET * CLSIZE;
    if (friendly) {
      cacheline_alloc = alloc_friendly(alloc_size, pagesize);
    } else {
      i = posix_memalign (&cacheline_alloc, pagesize, alloc_size);
      if (i != 0)
        exit(EXIT_FAILURE);
    }
    starttime = time();
    fill_cache (cacheline_alloc, alloc_size, iter);
    endtime = time();
                                                                               
    free(cacheline_alloc);
    printf ("Data set size: %zu, %d iter, %f seconds\n", 
            alloc_size, iter, endtime - starttime);
    exit (EXIT_SUCCESS);
  }
  

  if (friendly) {
    cacheline_alloc = alloc_friendly(NBUCKET * (CLSIZE + 1), pagesize);
  } else {
      i = posix_memalign (&cacheline_alloc, pagesize, NBUCKET * (CLSIZE + 1));
      if (i != 0)
        exit(EXIT_FAILURE);
  }
  if (1) {
    nr = cacheline_alloc;
    cacheline = (unsigned int *)((unsigned char *) cacheline_alloc + NBUCKET);
  } else {
    cacheline = cacheline_alloc; 
    nr = (unsigned char *) cacheline_alloc + NBUCKET * CLSIZE;
  }

  buckets = (unsigned int **) malloc(NBUCKET * sizeof (unsigned int *));
  {
    void *alloc;
    i = posix_memalign (&alloc, 4096, NBUCKET * (BUCKET_SIZE + bucket_skip) * sizeof(unsigned int));
    if (i != 0)
      exit(EXIT_FAILURE);
    bucket_storage = (unsigned int *) alloc;
  }
  /* Have system allocate everything before we start timing */
  memset (bucket_storage, 0, NBUCKET * (BUCKET_SIZE + bucket_skip) * sizeof(unsigned int));
  
  for (which = 0; which < 2; which++) {
    if (do_which[which]) {
      starttime = time();
      for (j = 0; j < iter; j++) {
        for (i = 0; i < NBUCKET; i++)
          buckets[i] = bucket_start (bucket_storage, i);

        if (which == 0)
          fill(buckets, j + 100);
        else 
          fill_cacheline(buckets, j + 100, cacheline_alloc);
      }
      endtime = time();
      printf ("%s cache line buffers: iter=%d, fill_iter=%lu, %f seconds\n", 
              which_str[which], iter, fill_iter, endtime - starttime);

      // printpointers(buckets, bucket_storage);
    }
  }

  free (cacheline_alloc);
  exit (EXIT_SUCCESS);
}
