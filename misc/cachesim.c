#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "portability.h"

#define UNUSED 0xFFFFFFFFFFFFFFFFUL

/*
Terminology:
With a 2-way set associative cache of 64 kbytes and 64 bytes cache line size,
there are 1024 lines and 512 cache indices, with 2 lines mapping to one index.

+-----------+-----------+
| line 0    | line 512  | index 0
+-----------+-----------+
| line 1    | line 513  | index 1
+-----------+-----------+
...
+-----------+-----------+
| line 511  | line 1023 | index 511
+-----------+-----------+

Each line stores the 64 bytes of cached data, the 'tag' which is the the 
memory address of the cached data (stored without the bits that are already 
determined by cache line size and index), and least-recently-used info which 
stores the relative 'age' since last access of each of the 2 cache lines at 
a given index.

For a memory access to address x, we first get l = floor(x/64), the line 
aligned part of the address; this l can map to any of the 2 lines at index 
i = l % 512. 

Which of the 2 lines, if any, hold cached data for the memory address l*64 
can be checked by comparing floor(l/512) with the cache tag for each of 
those 2 lines. 

If there is a match, the access hits in cache, the LRU data is updated, 
the cached data is sent to the CPU, and nothing more happens to the state 
of the cache. 

If the access misses in cache, one cacheline worth of data is transferred 
from memory addresses [l*64, l*64+63] to some cache line at index i. Which 
of the 2 possible cachelines at index i gets chosen is determined by a least-
recently-used (LRU) algorithm. Any access to a cache line at index i, hit or 
miss, updates the LRU data for all the cachelines at index i, i.e., the LRU 
data tells the order in which the cache lines at a given index were last 
accessed. When a miss occurs, the least-recently-used (oldest) one of the 
cache lines at index i gets replaced, the LRU info is updated, and l/512 is 
written to the replaced line's tag.
*/

typedef struct {
  int associativity, lines, linesize, indices, verbose, size;
  size_t *tag;
  int *lru;
  unsigned long hits, misses;
} cache_t;

typedef struct {
  size_t pagesize, nr_pages;
  size_t *pages;
} pages_t;

/* For a given index in the cache, return the line at this index that was 
   least recently used */
int find_lru(const cache_t *cache, const int index)
{
  int i;
  
  assert (index < cache->indices);
  /* Check the lines from the associative set to find the least recently used one.
     The age is increasing, so age == associativity - 1 is the oldest */
  for (i = 0; i < cache->associativity; i++) {
    int line = index + i * cache->indices;
    int age = cache->lru[line];
    assert (0 <= age && age < cache->associativity);
    if (age == cache->associativity - 1)
      return line;
  }
  abort();
}


/* Set a line's age to 0; lines in the same index that were younger 
   before, are aged by 1  */
void update_lru(cache_t *cache, const int line)
{
  const int index = line % cache->indices;
  int i, age;
  
  assert (line < cache->lines);
  age = cache->lru[line];
  assert (age < cache->associativity);
  
  for (i = 0; i < cache->associativity; i++) {
    int update_line = index + i * cache->indices;
    if (cache->lru[update_line] < age)
      cache->lru[update_line]++;
  }
  cache->lru[line] = 0;
}


/* Check LRU data for integrity. Among the lines at each index,
   each age 0, .., associativity-1 must occur exactly once */
void check_lru(const cache_t *cache)
{
  int i, j;
  int *counts;
  
  counts = malloc(cache->associativity * sizeof(int));
  
  for (i = 0; i < cache->indices; i++) {
    memset (counts, 0, cache->associativity * sizeof(int));
    for (j = 0; j < cache->associativity; j++) {
      int line = i + j * cache->indices;
      int age = cache->lru[line];
      assert (0 <= age && age < cache->associativity);
      counts[age]++;
    }
    for (j = 0; j < cache->associativity; j++) {
      assert (counts[j] == 1);
    }
  }
  
  free(counts);
}

/* Search the cache tags in the associative set to which address maps */
int is_in_cache(cache_t *cache, size_t line_address)
{
  const int index = line_address % cache->indices;
  const size_t address_tag = line_address / cache->indices;
  int i;
  
  for (i = 0; i < cache->associativity; i++) {
    int line = index + i * cache->indices;
    if (cache->tag[line] == address_tag)
      return line;
  }
  return cache->lines;
}

/* Returns 1 if an address was found in cache. If not, the address is added 
   to the cache, replacing another entry */
void access_cache(cache_t *cache, size_t memory_address)
{
  const size_t line_address = memory_address / cache->linesize;
  const int index = line_address % cache->indices;
  const size_t address_tag = line_address / cache->indices;
  size_t oldaddr;
  int line;
  
  line = is_in_cache(cache, line_address);
  if (line < cache->lines) {
    int age = cache->lru[line];
    cache->hits++;
    update_lru(cache, line);
    if (cache->verbose >= 2) {
      printf ("Access to address %zu (line address %zu) hit in cache line %d, age %d\n", 
              memory_address, line_address, line, age);
    }
    return;
  }

  /* Get the oldest line in this index */
  line = find_lru(cache, index);
  oldaddr = cache->tag[line];
  /* and set it to this address' tag */
  cache->tag[line] = address_tag;
  update_lru(cache, line);
  cache->misses++;
  if (cache->verbose >= 2) {
    if (oldaddr == UNUSED)
      printf ("Access to address %zu (line address %zu) missed in cache, entered in line %d (was unused)\n", 
              memory_address, line_address, line);
    else
      printf ("Access to address %zu (line address %zu) missed in cache, entered in line %d, replacing line address %zu\n", 
              memory_address, line_address, line, oldaddr);
  }
}

void cache_init(cache_t *cache, const int lines, const int associativity, 
                const int linesize, const int verbose) 
{
  int i;
  assert (lines % associativity == 0);
  
  cache->associativity = associativity;
  cache->lines = lines;
  cache->linesize = linesize;
  cache->size = lines * linesize;
  cache->indices = lines / associativity;
  cache->tag = malloc(lines * sizeof (size_t));
  cache->lru = malloc(lines * sizeof (int));
  cache->hits = 0;
  cache->misses = 0;
  cache->verbose = verbose;

  /* Init the LRU buffer to a state that satisfies our invariant:
     among the lines in each index, each age in [0, associativity-1] 
     occurs exactly once */
  for (i = 0; i < lines; i++) {
    cache->lru[i] = i / cache->indices;
  }
  
  for (i = 0; i < lines; i++) {
    cache->tag[i] = UNUSED;
  }
}

void cache_printstats(cache_t *cache)
{
  int i;
  
  if (cache->verbose >= 1) {
    printf ("Addresses in cache:\n");
    for (i = 0; i < cache->lines; i++) {
      printf ("%zu", cache->tag[i]);
      if (i % 20 == 19)
        printf ("\n");
      else
        printf ("\t");
    }
  }
  
  printf ("\n%lu hits, %lu misses\n", cache->hits, cache->misses);
}

void cache_clear(cache_t *cache) 
{
  free (cache->tag);
  cache->tag = NULL;
  free (cache->lru);
  cache->lru = NULL;
}

void pages_init(pages_t *pages, size_t pagesize, size_t nr_pages)
{
  size_t i;
  pages->pagesize = pagesize;
  pages->nr_pages = nr_pages;
  
  pages->pages = malloc(nr_pages * sizeof (size_t));
  for (i = 0; i < nr_pages; i++)
    pages->pages[i] = i; /* Simple linear layout for now */
}

void pages_clear(pages_t *pages)
{
  free (pages->pages);
  pages->pages = NULL;
}


int 
page_translate(pages_t *pages, size_t address)
{
  size_t vpage, PFN;
  
  vpage = address / pages->pagesize;
  assert (vpage < pages->nr_pages);
  PFN = pages->pages[vpage];
  
  return PFN * pages->pagesize + address % pages->pagesize;
}


void 
simulation(cache_t *cache, pages_t *pages, int argc, char **argv) {
  int i;
  int bucket_start = 123 * pages->pagesize;
  int region_start = 234 * pages->pagesize;
  int update, update_size = 4;
  int interleaved = 0;
  size_t bucket_read = bucket_start;
  size_t logical_address, physical_address;
  
  if (argc > 1) {
    interleaved = atoi(argv[1]);
  }
  if (argc > 2) {
    update_size = atoi(argv[2]);
  }
  if (interleaved)
    printf ("Using interleaved buckets\n");
  else
    printf ("Not using interleaved buckets\n");
  printf ("Using updates of size %d\n", update_size);
  
  for (i = 0; i < 1000000; i++) {
    /* Read an update from the bucket */
    
    logical_address = bucket_read;
    bucket_read += update_size;
    if (interleaved && bucket_read % cache->linesize == 0) {
      bucket_read = bucket_read + cache->size - cache->linesize;
    }
    physical_address = page_translate(pages, logical_address);
    access_cache (cache, physical_address);
    
    /* The address where update hits is assumed to be random in the bucket region */
    update = random() % (cache->size);
    
    /* Apply that update to the bucket region */
    logical_address = update + region_start;
    physical_address = page_translate(pages, logical_address);
    access_cache (cache, region_start + update);
  }
}


int main(int argc, char **argv) {
  cache_t cache[1];
  pages_t pages[1];
  int verbose = 0;
  int lines, associativity;
  int linesize = 64;
  
  while (argc > 1 && strcmp(argv[1], "-v") == 0) {
    argc--;
    argv++;
    verbose++;
  }
  
  if (argc < 4) {
    printf ("%s [-v] cachelines associativity linesize [parameters for simulation]\n", argv[0]);
    exit (EXIT_FAILURE);
  }
  
  lines = atoi (argv[1]); 
  argc--; 
  argv++;
  associativity = atoi (argv[1]); 
  argc--; 
  argv++;
  linesize = atoi (argv[1]); 
  argc--; 
  argv++;
  
  cache_init (cache, lines, associativity, linesize, verbose);
  check_lru(cache);
  pages_init(pages, 4096, 1UL<<20);

  simulation(cache, pages, argc, argv);
  
  check_lru(cache);
  cache_printstats(cache);

  cache_clear (cache);
  pages_clear(pages);
  exit (EXIT_FAILURE);
}
