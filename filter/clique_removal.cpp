#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

#include "portability.h"

#include "utils_with_io.h"
#include "filter_config.h"
#include "purge_matrix.h"
#include "clique_removal.h"

/********************* comp_t struct (clique) ********************************/

typedef struct {
  uint64_t i; /* smallest row appearing in the connected component */
  float w;   /* Weight of the connected component */
} comp_t;

/* Compare weights of two connected components.
   Return 1 if the weight of c1 is less than the weight of c2.
   Return -1 if the weight of c1 is larger than the weight of c2.
   If the weights are the same
      return 1 if the first row of c1 is less than the first row of c2;
      return -1 if the first row of c1 is not less than the first row of c2.
   We do that in order to be deterministic (even with different number of
   threads).
 */
static inline int
comp_cmp_weight (comp_t c1, comp_t c2)
{
  if (c1.w < c2.w)
    return 1;
  else if (c1.w > c2.w)
    return -1;
  else /* same weight, use first row to be deterministic */
  {
    if (c1.i < c2.i)
      return 1;
    else
      return -1;
  }
}

static inline int
comp_weight_is_smaller (comp_t c1, comp_t c2)
{
  return (comp_cmp_weight (c1, c2) > 0);
}

static inline int
comp_weight_is_larger (comp_t c1, comp_t c2)
{
  return (comp_cmp_weight (c1, c2) < 0);
}

/* Compare weights of two connected components (to use with qsort) */
static int
comp_cmp_weight_for_qsort (const void *p, const void *q)
{
  return comp_cmp_weight (*((comp_t *)p), *((comp_t *)q));
}

/* Contribution of each column of weight w to the weight of the connected
   component.
   Reference [1]: "The filtering step of discrete logarithm and integer
   factorization algorithms", Cyril Bouvier, https://hal.inria.fr/hal-00734654,
   2013, 27 pages.
   Each weight is identified by LAMBDA (0 to 6) and NU (0 to 3),
   the default one is \Omega_{31} (LAMBDA=3, NU=1).
   Cavallar's weight function is \Omega_{23} (LAMBDA=2, NU=3).
*/

#ifndef USE_WEIGHT_LAMBDA
#define USE_WEIGHT_LAMBDA 3
#endif
#ifndef USE_WEIGHT_NU
#define USE_WEIGHT_NU 1
#endif

static inline float
comp_weight_function (weight_t w)
{
  if (w >= 3)
#if USE_WEIGHT_LAMBDA == 0
    return 1.0;
#elif USE_WEIGHT_LAMBDA == 1
    return powf (2.0 / 3.0, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 2
    return powf (0.5, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 3
    return powf (0.8, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 4
    return 1.0 / log2f ((float) w);
#elif USE_WEIGHT_LAMBDA == 5
    return 2.0 / (float) w;
#elif USE_WEIGHT_LAMBDA == 6
    return 4.0 / (float) (w * w);
#else
#error "Invalid value of USE_WEIGHT_LAMBDA"
#endif
  else if (w == 2) /* since each ideal is counted twice, we put here half
                      the values of \nu_i from reference [1] */
#if USE_WEIGHT_NU == 0
    return 0.0;
#elif USE_WEIGHT_NU == 1
    return 0.125;
#elif USE_WEIGHT_NU == 2
    return 0.25;
#elif USE_WEIGHT_NU == 3
    return 0.5;
#else
#error "Invalid value of USE_WEIGHT_NU"
#endif
  else
    return 0.0;
}

/* print info on the weight function that is used */
void
comp_print_info_weight_function ()
{
  fprintf(stdout, "Weight function used during clique removal:\n"
                  "  0     1     2     3     4     5     6     7\n"
                  "0.000 0.000 ");
  for (unsigned int k = 2; k < 8; k++)
    fprintf(stdout, "%0.3f ", comp_weight_function((weight_t) k));
  fprintf(stdout, "\n");
}

/***************** Sorted binary tree of comp_t struct ***********************/

/* This is used for sorting the connected components according to their weigth.
   These are internal functions.
 */

/* Binary tree of comp_t sorted by weight of the connected component.
 * The two sons of a node of index n are the nodes of indexes 2n+1 and 2n+1.
 * The property is: the weight of the comp_t father is less than its 2 sons.
 * so root of the tree is the connected component with the smallest weight.
 * When the tree is full, a new node is inserted in the tree iff its weight is
 * greater than the weight of the root (and the root is removed from the tree).
 * At the end the tree contained the 'size' (= 'alloc') heaviest comp_t.
 */
struct comp_sorted_bin_tree_s {
  size_t alloc;
  size_t size;
  comp_t * tree;
};
typedef struct comp_sorted_bin_tree_s comp_sorted_bin_tree_t[1];
typedef struct comp_sorted_bin_tree_s * comp_sorted_bin_tree_ptr;
typedef const struct comp_sorted_bin_tree_s * comp_sorted_bin_tree_srcptr;

static inline void
comp_sorted_bin_tree_init (comp_sorted_bin_tree_ptr T, size_t max_size)
{
  /* If max_size is even, T->alloc = max_size + 1, in order to avoid a
   * node-father with only one node-son.a We will have one more connected
   * component than needed but we do not care (the computing cost in
   * negligible) */
  if (max_size & 1)
    T->alloc = max_size;
  else
    T->alloc = max_size + 1;

  T->tree = (comp_t *) malloc_check (sizeof(comp_t) * T->alloc);
  T->size = 0;
}

static inline void
comp_sorted_bin_tree_clear (comp_sorted_bin_tree_ptr T)
{
  free(T->tree);
  T->size = T->alloc = 0;
}

static inline void
comp_sorted_bin_tree_qsort (comp_sorted_bin_tree_ptr T,
                            int (*cmp_fct)(const void *, const void *))
{
  qsort (T->tree, T->size, sizeof(comp_t), cmp_fct);
}

static inline void
comp_sorted_bin_tree_insert (comp_sorted_bin_tree_ptr T, comp_t new_node)
{
  /* If the tree is not full, add new_node at the end and go up. */
  if (UNLIKELY (T->size < T->alloc))
  {
    size_t cur = T->size;
    while (cur)
    {
      size_t father = (cur - 1) >> 1;
      if (UNLIKELY (comp_weight_is_larger (new_node, T->tree[father])))
        break;
      T->tree[cur] = T->tree[father];
      cur = father;
    }
    T->size++;
    T->tree[cur] = new_node;
  }
  /* If the tree is full, insert new_node iff the weight of new_node if greater
   * than the weight of the root. In this case, replace the root by new node
   * and go down */
  else if (UNLIKELY(comp_weight_is_larger (new_node, T->tree[0])))
  {
    size_t cur = 0;
    for (;;) /* In this loop, always, T->tree[cur].w < new_node.w */
    {
      size_t lightest_son = cur * 2 + 1; /* son is son1 */
      /* If there is no son: cur is a leaf. Found! */
      if (UNLIKELY (lightest_son >= T->alloc))
        break;
      if (comp_weight_is_smaller(T->tree[lightest_son+1],T->tree[lightest_son]))
      {
        lightest_son++; /* lightest_son is now son 2 */
      }
      /* Now lightest_son and lightest_weight are correct */
      if (UNLIKELY (comp_weight_is_smaller (new_node, T->tree[lightest_son])))
        break; /* Found! */
      /* The lightest_son has a weight lighter than new_node, so we have to go
       * down the tree, but before the son replaces its father. */
      T->tree[cur] = T->tree[lightest_son];
      cur = lightest_son;
    }
    T->tree[cur] = new_node;
  }
  /* else do nothing */
}

/*********** Functions to compute and delete connected component **************/

template <class T>
bool vector_contains(std::vector<T> const & v, const T& val)
{
    return std::find(v.begin(), v.end(), val) != v.end();
}

/* Compute connected component beginning at row clique->i
 * The weight of the connected component is written in clique->w
 *
 * Multithread safe (if each thread has its own row_buffer).
 *
 * Return the number of rows in the connected component or 0 if it exists a
 * row i1 with i1 < clique->i in the connected component (it implies that this
 * component has already been found when the function was called with a smaller
 * clique->i).
 */
uint64_t
compute_one_connected_component (comp_t *clique, purge_matrix_srcptr mat,
                                 std::vector<uint64_t> & row_buffer)
{
  row_buffer.clear();
  row_buffer.push_back(clique->i);

  clique->w = 0.; /* Set initial weight of the connected component to 0. */

  /* Loop on all connected rows */
  for(size_t k = 0 ; k < row_buffer.size() ; k++) {
    uint64_t cur_row = row_buffer[k];

    /* Loop on all columns of the current row */
    for (index_t * h = mat->row_compact[cur_row]; *h != UMAX(*h); h++)
    {
      index_t cur_h = *h;
      weight_t cur_h_weight = mat->cols_weight[cur_h];
      clique->w += comp_weight_function (cur_h_weight);
      if (UNLIKELY(cur_h_weight == 2))
      {
        ASSERT_ALWAYS(cur_row <= mat->sum2_row[cur_h]);
        uint64_t the_other_row = mat->sum2_row[cur_h] - cur_row;
        /* First, if the_other_row < clique.i, the connected component was
         * already found (by this thread or another). return 0 */
        if (the_other_row < clique->i)
          return 0;
        /* If the_other_row is not already in the buffer, add it as a todo. */
        if (!vector_contains(row_buffer, the_other_row))
            row_buffer.push_back(the_other_row);
      }
    }
  }

  return row_buffer.size();
}


/* Delete the connected component containing the row "current_row"
 *
 * WARNING: this code itself is multithread compatible (if each thread has its
 * own row_buffer), but it calls purge_matrix_delete_row, which is NOT
 * compatible!
 */
void
delete_one_connected_component (purge_matrix_ptr mat, uint64_t cur_row,
                                std::vector<uint64_t> & row_buffer)
{
  row_buffer.clear();
  row_buffer.push_back(cur_row);

  /* Loop on all connected rows */
  for(size_t k = 0 ; k < row_buffer.size() ; k++) {
    uint64_t cur_row = row_buffer[k];
    /* Loop on all columns of the current row */
    for (index_t * h = mat->row_compact[cur_row]; *h != UMAX(*h); h++)
    {
      index_t cur_h = *h;
      weight_t cur_h_weight = mat->cols_weight[cur_h];
      /* We might have some H->hashcount[h] = 3, which is decreased to 2, but
       * we don't want to treat that case. Thus we check in addition that
       * mat->sum2_row[h] <> 0, which only occurs when H->hashcount[h] = 2
       * initially.*/
      if (UNLIKELY(cur_h_weight == 2 && mat->sum2_row[cur_h]))
      {
        ASSERT_ALWAYS(cur_row <= mat->sum2_row[cur_h]);
        uint64_t the_other_row = mat->sum2_row[cur_h] - cur_row;
        /* If the_other_row is not already in the buffer, add it as a todo. */
        if (!vector_contains(row_buffer, the_other_row))
            row_buffer.push_back(the_other_row);
      }
    }
  }

  /* Now, we deleted all rows explored */
  for (uint64_t pt : row_buffer)
    purge_matrix_delete_row (mat, pt);
    /* mat->nrows and mat->ncols are updated by purge_matrix_delete_row. */
}


/************* Code for clique (= connected component) removal ***************/

/* Print stats on the length of the cliques (connected components). The stats
 * are expensive to compute, so these functions should not be called by default.
 * Assume mat->sum2_row is already computed.
 */
void
purge_matrix_print_stats_on_cliques (FILE *out, purge_matrix_srcptr mat,
                                     int verbose)
{
  uint64_t *len = NULL;

  len = (uint64_t *) malloc (mat->nrows_init * sizeof (uint64_t));
  ASSERT_ALWAYS (len != NULL);
  memset (len, 0, mat->nrows_init * sizeof (uint64_t));

  /* use a single buffer, avoid constant reallocations */
  std::vector<uint64_t> buf;

  for (uint64_t i = 0; i < mat->nrows_init; i++)
  {
    if (purge_matrix_is_row_active (mat, i))
    {
      comp_t c = {.i = i, .w = 0.0};
      uint64_t nrows = compute_one_connected_component (&c, mat, buf);
      len[i] = nrows;
    }
  }

  print_stats_uint64 (out, len, mat->nrows_init, "cliques", "length", verbose);

  free (len);
}


/* Code for clique_removal --- multithread version */

/* The main structure for the working pthreads pool which compute the
   heaviest cliques */
typedef struct comp_mt_thread_data_s {
  /* Read only part */
  uint64_t begin_first_block;
  uint64_t end_first_block;
  uint64_t jump_to_next_block;
  purge_matrix_srcptr mat;
  /* Read-write part */
  comp_sorted_bin_tree_t comp_tree; /* sorted tree computed by the thread */
} comp_mt_data_t;

/* Size of the block of rows treated by a thread during the computation of the
 * connected component (has to be of the form: BV_BITS << x)*/
#define COMP_MT_ROWS_BLOCK (BV_BITS<<4)

/* This multithread function fills the tree data->comp_tree with the
 * data->comp_tree->alloc heaviest connected components whose smallest row
 * belongs in
 *   [ data->begin_first_block + j * jump_to_next_block,
                             data->end_first_block + j * jump_to_next_block ]
 * with j=0,1,2.. until the interval does not intersect [0..mat->nrows_init-1].
 */
static void *
clique_removal_core_mt_thread (void *pt)
{
  comp_mt_data_t *data = (comp_mt_data_t *) pt;
  uint64_t begin_cur_block, end_cur_block;
  comp_t clique;

  std::vector<uint64_t> buf;

  begin_cur_block = data->begin_first_block;
  end_cur_block = data->end_first_block;
  while (LIKELY(begin_cur_block < data->mat->nrows_init))
  {
    if (UNLIKELY (end_cur_block > data->mat->nrows_init))
      end_cur_block = data->mat->nrows_init;

    for (uint64_t i = begin_cur_block; i < end_cur_block; i++)
    {
      if (purge_matrix_is_row_active (data->mat, i))
      {
        clique.i = i;
        unsigned int nb_rows = compute_one_connected_component (&(clique),
                                              data->mat, buf);
        if (nb_rows == 0) /* this component was already found earlier */
          continue;
        comp_sorted_bin_tree_insert (data->comp_tree, clique);
      }
    }

    begin_cur_block += data->jump_to_next_block;
    end_cur_block += data->jump_to_next_block;
  }

  /* Re-order the connected component by decreasing weight. */
  comp_sorted_bin_tree_qsort (data->comp_tree, comp_cmp_weight_for_qsort);

  return NULL;
}

static uint64_t
clique_removal_core_mt (purge_matrix_ptr mat, int64_t target_excess,
                        size_t max_nb_comp, unsigned int nthreads, int verbose)
{
  pthread_t *threads = NULL;
  comp_mt_data_t *th_data = NULL;
  /* nrows_per_block must be a multiple of BV_BITS */
  uint64_t k, nrows_per_block = COMP_MT_ROWS_BLOCK;
  int err;
  uint64_t nb_cliques_del = 0;

  th_data = (comp_mt_data_t *) malloc (nthreads * sizeof(comp_mt_data_t));
  ASSERT_ALWAYS(th_data != NULL);
  threads = (pthread_t *) malloc (nthreads * sizeof(pthread_t));
  ASSERT_ALWAYS(threads != NULL);

  k = 0;
  for (unsigned int i = 0; i < nthreads; i++)
  {
    th_data[i].mat = mat;
    th_data[i].begin_first_block = k;
    k += nrows_per_block;
    th_data[i].end_first_block = k;
    th_data[i].jump_to_next_block = nthreads * nrows_per_block;
    comp_sorted_bin_tree_init (th_data[i].comp_tree, max_nb_comp);
  }

  for (unsigned int i = 0; i < nthreads; i++)
  {
    if ((err = pthread_create (&(threads[i]), NULL,
                     &clique_removal_core_mt_thread, (void *) &(th_data[i]))))
    {
      fprintf(stderr, "Error, pthread_create failed in clique_removal_core_mt: "
                      "%d. %s\n", err, strerror(errno));
      abort();
    }
  }
  for (unsigned int i = 0; i < nthreads; i++)
    pthread_join (threads[i], NULL);

  fprintf(stdout, "Cliq. rem.: computed heaviest connected components at "
                  "%2.2lf\n", seconds());
  fflush (stdout);

  if (verbose > 0)
    purge_matrix_print_stats_on_cliques (stdout, mat, verbose);

  /* At this point, in each pth[i].comp_tree we have pth[i].comp_tree->size
     connected components ordered by decreasing weight. */
  size_t *next_clique = NULL;
  std::vector<uint64_t> buf;
  next_clique = (size_t *) malloc (nthreads * sizeof (next_clique));
  ASSERT_ALWAYS (next_clique != NULL);
  memset (next_clique, 0, nthreads * sizeof (next_clique));

  while (mat->nrows > target_excess + mat->ncols)
  {
    comp_t * max_comp = NULL;
    size_t max_thread = 0;

    for (unsigned int i = 0; i < nthreads; ++i)
    {
      if (next_clique[i] < th_data[i].comp_tree->size)
      {
        comp_t * cur_comp = &(th_data[i].comp_tree->tree[next_clique[i]]);
        if (max_comp == NULL || comp_weight_is_larger (*cur_comp, *max_comp))
        {
          max_comp = cur_comp;
          max_thread = i;
        }
      }
    }
    if (max_comp == NULL)
    {
      fprintf (stderr, "Cliq. rem.: Warning, all lists of connected components"
                       " are empty\n");
      break;
    }

    delete_one_connected_component (mat, max_comp->i, buf);
    next_clique[max_thread]++;
    nb_cliques_del++;
  }

  free (next_clique);

  /* We can free th_data[i].comp_tree and th_data itself */
  for (unsigned int i = 0; i < nthreads; ++i)
    comp_sorted_bin_tree_clear (th_data[i].comp_tree);
  free (th_data);
  free (threads);

  return nb_cliques_del;
}

/* Code for clique_removal --- monothread version */

static uint64_t
clique_removal_core_mono (purge_matrix_ptr mat, int64_t target_excess,
                          size_t max_nb_comp, int verbose)
{
  comp_sorted_bin_tree_t comp_tree;
  std::vector<uint64_t> buf;
  comp_t clique;

  comp_sorted_bin_tree_init (comp_tree, max_nb_comp);

  for (clique.i = 0; clique.i < mat->nrows_init; clique.i++)
  {
    if (purge_matrix_is_row_active (mat, clique.i))
    {
      unsigned int nb_rows = compute_one_connected_component (&(clique),
                                            mat, buf);
      if (nb_rows == 0) /* this component was already found earlier */
        continue;
      comp_sorted_bin_tree_insert (comp_tree, clique);
    }
  }

  /* Re-order the connected component by decreasing weight. */
  comp_sorted_bin_tree_qsort (comp_tree, comp_cmp_weight_for_qsort);

  fprintf(stdout, "Cliq. rem.: computed heaviest connected components at "
                  "%2.2lf\n", seconds());
  fflush (stdout);

  if (verbose > 0)
    purge_matrix_print_stats_on_cliques (stdout, mat, verbose);

  /* At this point, comp_tree contains max_nb_comp connected components order
   * by decreasing weight.
   */
  size_t next_clique = 0;
  uint64_t nb_clique_deleted = 0;

  while (mat->nrows > target_excess + mat->ncols)
  {
    if (next_clique >= comp_tree->size)
    {
      fprintf (stderr, "Cliq. rem.: Warning, the list of connected components"
                       " is empty\n");
      break;
    }

    delete_one_connected_component (mat, comp_tree->tree[next_clique].i, buf);
    next_clique++;
    nb_clique_deleted++;
  }

  comp_sorted_bin_tree_clear (comp_tree);

  return nb_clique_deleted;
}

/***************************** Clique removal ********************************/

void
cliques_removal (purge_matrix_ptr mat, int64_t target_excess,
                 unsigned int nthreads, int verbose)
{
  int64_t excess = purge_matrix_compute_excess (mat);
  uint64_t nb_cliques_del = 0;
  size_t max_nb_comp;

  /* If the excess is smaller than target_excess, then we have nothing to do. */
  if (excess <= target_excess)
    return;

  /* First, collect, for each column of weight 2, the sum i1+i2, where i1 and i2
   * are the two rows containing the column of weight 2.
   */
  purge_matrix_compute_sum2_row (mat, nthreads);
  fprintf(stdout, "Cliq. rem.: computed mat->sum2_row at %2.2lf\n", seconds());
  fflush (stdout);

  /* max_nb_comp is the maximum number of connected components that each threads
   * is going to store. If we are in monothread, we increase max_nb_comp by 25%
   * to take into account the fact that something a little bit more connected
   * components are needed to achieve the targeted excess.
   */
  max_nb_comp = (size_t) (excess - target_excess);
  if (nthreads == 1)
    max_nb_comp+= max_nb_comp/4;

  /* Then, call the core function (either mono or mt) that compute and delete
   * connected components until excess is equal to target_excess.
   */
  if (nthreads > 1)
    nb_cliques_del = clique_removal_core_mt (mat, target_excess, max_nb_comp,
                                             nthreads, verbose);
  else
    nb_cliques_del = clique_removal_core_mono (mat, target_excess, max_nb_comp,
                                               verbose);


  fprintf(stdout, "Cliq. rem.: deleted %" PRIu64 " heaviest connected "
                  "components at %2.2lf\n", nb_cliques_del, seconds());
  if (verbose > 0)
    fprintf(stdout, "# INFO: max_nb_comp_per_thread=%zu target_excess="
                    "%" PRId64 "\n", max_nb_comp, target_excess);
  fflush (stdout);
}

