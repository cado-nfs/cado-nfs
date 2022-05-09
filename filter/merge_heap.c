#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "merge_heap.h"
#include "merge_replay_matrix.h"  // typerow_t, index_t, 
#include "sparse.h" // rowCell
#include "omp_proxy.h"
#include "timing.h"  // seconds

/*************************** heap structures *********************************/

/* 
  Threads allocate PAGES of memory of a fixed size to store rows, and each
  thread has a single ACTIVE page in which it writes new rows. When the active
  page is FULL, the thread grabs an EMPTY page that becomes its new active page.
  In each page, there is a pointer [[ptr]] to the begining of the free space. To
  allocate [[b]] bytes for a new row, it suffices to note the current value of
  [[ptr]] and then to increase it by [[b]] --- if this would overflow the current
  page, then it is marked as "full" and a new active page is obtained. Rows are
  stored along with their number, their size and the list of their coefficients.
  To delete a row, we just mark it as deleted by setting its number to -1. Thus,
  row allocation and deallocation are thread-local operations that are very
  fast. The last allocated row can easily be shrunk (by diminishing [[ptr]]).

  After each pass, memory is garbage-collected. All threads do the following
  procedure in parallel, while possible: grab a full page that was not created
  during this pass; copy all non-deleted rows to the current active page; mark the
  old full page as "empty". Note that moving row $i$ to a different address
  in memory requires an update to the ``row pointer'' in [[rows]]; this is
  why rows are stored along with their number. When they need a new page, threads
  first try to grab an existing empty page. If there is none, a new page is
  allocated from the OS. There are global doubly-linked lists of full and empty
  pages, protected by a lock, but these are infrequently accessed.
*/

#define PAGE_DATA_SIZE ((1<<18) - 4) /* seems to be optimal for RSA-512 */

struct page_t {
        struct pagelist_t *list;     /* the pagelist_t structure associated with this page */
        int i;                       /* page number, for debugging purposes */
        int generation;              /* pass in which this page was filled. */
        int ptr;                     /* data[ptr:PAGE_DATA_SIZE] is available*/
        typerow_t data[PAGE_DATA_SIZE];
};

// linked list of pages (doubly-linked for the full pages, simply-linked for the empty pages)
struct pagelist_t {
        struct pagelist_t *next;
        struct pagelist_t *prev;
        struct page_t *page;
};


int heap_config_get_PAGE_DATA_SIZE()
{
    return PAGE_DATA_SIZE;
}

static int n_pages, n_full_pages, n_empty_pages;
static struct pagelist_t headnode;                  // dummy node for the list of full pages
static struct pagelist_t *full_pages, *empty_pages; // head of the page linked lists
static struct page_t **active_page; /* active page of each thread */
static long long *heap_waste;       /* space wasted, per-thread. Can be negative! The sum over all threads is correct. */
static int current_generation;


/* provide an empty page */
static struct page_t *
heap_get_free_page()
{
        struct page_t *page = NULL;
        #pragma omp critical(pagelist)
        {
                // try to grab it from the simply-linked list of empty pages.
                if (empty_pages != NULL) {
                        page = empty_pages->page;
                        empty_pages = empty_pages->next;
                        n_empty_pages--;
                } else {
                        n_pages++;   /* we will malloc() it, update count while still in critical section */
                }
        }
        if (page == NULL) {
                // we must allocate a new page from the OS.
                page = malloc(sizeof(struct page_t));
                struct pagelist_t *item = malloc(sizeof(struct pagelist_t));
                page->list = item;
                item->page = page;
                page->i = n_pages;
        }
        page->ptr = 0;
        page->generation = current_generation;
        return page;
}

/* Provide the oldest full page with generation < current_generation, or NULL if none is available,
   and remove it from the doubly-linked list of full pages */
static struct page_t *
heap_get_full_page()
{
        struct pagelist_t *item = NULL;
        struct page_t *page = NULL;
        #pragma omp critical(pagelist)
        {
                item = full_pages->next;
                if (item->page != NULL && item->page->generation < current_generation) {
                        page = item->page;
                        item->next->prev = item->prev;
                        item->prev->next = item->next;
                        n_full_pages--;
                }
        }
        return page;
}


/* declare that the given page is empty */
static  void
heap_clear_page(struct page_t *page)
{
        struct pagelist_t *item = page->list;
        #pragma omp critical(pagelist)
        {
                item->next = empty_pages;
                empty_pages = item;
                n_empty_pages++;
        }
}

/* declare that the given page is full. Insert to the left of the list of full pages.
   The list is sorted (following next) by increasing generation. */
static void
heap_release_page(struct page_t *page)
{
        struct pagelist_t *list = page->list;
        struct pagelist_t *target;
        #pragma omp critical(pagelist)
        {
                target = full_pages->prev;
                while (target->page != NULL && page->generation < target->page->generation)
                        target = target->prev;
                list->next = target;
                list->prev = target->prev;
                list->next->prev = list;
                list->prev->next = list;
                n_full_pages++;
        }
}

// set up the page linked lists
void
heap_setup()
{
        current_generation = 0;

        // setup the doubly-linked list of full pages.
        full_pages = &headnode;
        full_pages->page = NULL;
        full_pages->next = full_pages;
        full_pages->prev = full_pages;

        empty_pages = NULL;
        int T = omp_get_max_threads();
        active_page = malloc(T * sizeof(*active_page));
        heap_waste = malloc(T * sizeof(*heap_waste));

        #pragma omp parallel for
        for(int t = 0 ; t < T ; t++) {
            active_page[t] = heap_get_free_page();
            heap_waste[t] = 0;
        }
}

void
heap_clear ()
{
  /* clear active pages */
  int T = omp_get_max_threads ();
  for (int t = 0 ; t < T ; t++) {
    free(active_page[t]->list);
    free(active_page[t]);
  }

  /* clear empty pages */
  while (empty_pages != NULL) {
    struct pagelist_t *item = empty_pages;
    empty_pages = item->next;
    free(item->page);
    free(item);
  }

  /* clear full pages. 1. Locate dummy node */
  while (full_pages->page != NULL)
    full_pages = full_pages->next;

  // 2. Skip dummy node
  full_pages = full_pages->next;

  // 3. Walk list until dummy node is met again, free everything.
  while (full_pages->page != NULL) {
    struct pagelist_t *item = full_pages;
    full_pages = full_pages->next;
    free(item->page);
    free(item);
  }
  free (active_page);
  free (heap_waste);
}


/* Returns a pointer to allocated space holding a size-s array of typerow_t.
   This function is thread-safe. */
static inline typerow_t * heap_malloc (size_t s)
{
  ASSERT(s <= PAGE_DATA_SIZE);
  int t = omp_get_thread_num();
  struct page_t *page = active_page[t];
  // ASSERT(page != NULL);
  /* enough room in active page ?*/
  if (page->ptr + s >= PAGE_DATA_SIZE) {
        heap_release_page(page);
        page = heap_get_free_page();
        active_page[t] = page;
  }
  typerow_t *alloc = page->data + page->ptr;
  page->ptr += s;
  return alloc;
}


typerow_t *
heap_alloc_row (index_t i, size_t s)
{
  typerow_t *alloc = heap_malloc(s + 2);
  rowCell(alloc, 0) = i;
  rowCell(alloc, 1) = s;
  return alloc + 1;
}

void
heap_resize_last_row (typerow_t *row, index_t new_size)
{
  int t = omp_get_thread_num();
  struct page_t *page = active_page[t];
  index_t old_size = rowCell(row, 0);
  ASSERT(row + old_size + 1 == page->data + page->ptr);
  int delta = old_size - new_size;
  rowCell(row, 0) = new_size;
  page->ptr -= delta;
}



void
heap_destroy_row(typerow_t *row)
{
  int t = omp_get_thread_num();
  rowCell(row, -1) = (index_signed_t) -1;
  heap_waste[t] += rowCell(row, 0) + 2;
}

/* Copy non-garbage data to the active page of the current thread, then
   return the page to the freelist. Thread-safe.
   Warning: this may fill the current active page and release it. */
static int
collect_page(typerow_t **rows, struct page_t *page)
{
        int garbage = 0;
        int bot = 0;
        int top = page->ptr;
        typerow_t *data = page->data;
        while (bot < top) {
                index_signed_t i = rowCell(data, bot);
                typerow_t *old = data + bot + 1;
                index_t size = rowCell(old, 0);
                if (i == (index_signed_t) -1) {
                        garbage += size + 2;
                } else {
                        ASSERT(rows[i] == old);
                        typerow_t * new = heap_alloc_row(i, size);
                        memcpy(new, old, (size + 1) * sizeof(typerow_t));
                        setCell(new, -1, rowCell(old, -1), 0);
                        rows[i] = new;
                }
                bot += size + 2;
        }
        int t = omp_get_thread_num();
        heap_waste[t] -= garbage;
        heap_clear_page(page);
        return garbage;
}

static double
heap_waste_ratio()
{
        int T = omp_get_max_threads();
        long long total_waste = 0;
        for (int t = 0; t < T; t++)
                total_waste += heap_waste[t];
        double waste = ((double) total_waste) / (n_pages - n_empty_pages) / PAGE_DATA_SIZE;
        return waste;
}

void
heap_garbage_collection(typerow_t **rows)
{
        double waste = heap_waste_ratio();
        printf("Starting collection with %.0f%% of waste...", 100 * waste);
        fflush(stdout);

        // I don't want to collect pages just filled during the collection
        current_generation++;

        int i = 0;
        int initial_full_pages = n_full_pages;
        struct page_t *page;
        long long collected_garbage = 0;
        #pragma omp parallel reduction(+:i, collected_garbage) private(page)
        while ((page = heap_get_full_page()) != NULL) {
                collected_garbage += collect_page(rows, page);
                i++;
        }

        double page_ratio = (double) i / initial_full_pages;
        double recycling = 1 - heap_waste_ratio() / waste;
        if (i == 0)
                i = 1; // avoid division by zero
        printf("Examined %.0f%% of full pages, recycled %.0f%% of waste. %.0f%% of examined data was garbage\n",
        	100 * page_ratio, 100 * recycling, 100.0 * collected_garbage / i / PAGE_DATA_SIZE);
}