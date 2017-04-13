#include "cado.h"

// #include <stdio.h>
// #include <stdlib.h>
// #include <fcntl.h>		/* for _O_BINARY */
// #include <string.h>		/* for strcmp */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>
#include <vector>
#include <iostream>
#include <sstream>

#include "portability.h"

/* We're going to use this transitionally */
#define MKZTYPE_CAVALLAR 0
#define MKZTYPE_PURE 1
#define MKZTYPE_LIGHT 2

#include "filter_config.h"
#include "utils_with_io.h"
// #include "merge_replay_matrix.h" /* for filter_matrix_t */
// #include "report.h"     /* for report_t */
// #include "markowitz.h" /* for MkzInit */
// #include "merge_mono.h" /* for mergeOneByOne */
// #include "sparse.h"

#include "medium_int.hpp"
#include "indexed_priority_queue.hpp"
#include "compressible_heap.hpp"
#include "get_successive_minima.hpp"
#define DEBUG_SMALL_SIZE_POOL
#include "small_size_pool.hpp"

/* NOTE: presently this value has a very significant impact on I/O speed
 * (about a factor of two when reading the purged file), and we have
 * reasons to believe that just about every operation will suffer from
 * the very same problems...
 */
static const int compact_column_index_size = 5;
static const int merge_row_heap_batch_size = 16384;

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "mat", "input purged file");
    param_list_decl_usage(pl, "out", "output history file");
    /*
    param_list_decl_usage(pl, "resume", "resume from history file");
    param_list_decl_usage(pl, "forbidden-cols",
			  "list of columns that cannot be "
			  "used for merges");
    param_list_decl_usage(pl, "path_antebuffer",
			  "path to antebuffer program");
    */
    param_list_decl_usage(pl, "force-posix-threads", "(switch)");
    param_list_decl_usage(pl, "v", "verbose level");
    param_list_decl_usage(pl, "t", "number of threads");
}

static void usage(param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


struct merge_matrix {
  size_t nrows;
  size_t ncols;
  size_t rem_nrows;     /* number of remaining rows */
  size_t rem_ncols;     /* number of remaining columns, including the buried
                             columns */
  size_t nburied;       /* the number of buried columns */
  size_t weight;        /* non-zero coefficients in the active part */
  size_t total_weight;  /* Initial total number of non-zero coefficients */
  size_t keep;          /* target for nrows-ncols */
  int maxlevel;         /* says it */
  unsigned int mkztype;
  unsigned int wmstmax;
  double target_density;

  merge_matrix()
  {
      rem_ncols = 0;
      total_weight = 0;
      maxlevel = DEFAULT_MERGE_MAXLEVEL;
      keep = DEFAULT_FILTER_EXCESS;
      nburied = DEFAULT_MERGE_SKIP;
      target_density = DEFAULT_MERGE_TARGET_DENSITY;
      mkztype = DEFAULT_MERGE_MKZTYPE;
      wmstmax = DEFAULT_MERGE_WMSTMAX;
  }
  /* {{{ merge_row_ideal type details */
#ifdef FOR_DL
  template<int index_size = 8>
  struct merge_row_ideal {
      typedef medium_int<index_size> index_type;
      // typedef size_t index_type;
      typedef int32_t exponent_type;
      /* exponent_type could be compressed, although we acknowledge the
       * fact that merge entails some increase of the exponent values,
       * obviously */
      private:
      index_type h;
      exponent_type e;
      public:
      merge_row_ideal() {}
      merge_row_ideal(prime_t p) : h(p.h), e(p.e) {}
      index_type const & index() const { return h; }
      index_type & index() { return h; }
      exponent_type const & exponent() const { return e; }
      exponent_type & exponent() { return e; }
      bool operator<(merge_row_ideal const& a) const { return h < a.h; }
      template<int otherwidth> friend class merge_row_ideal;
      template<int otherwidth>
      merge_row_ideal(merge_row_ideal<otherwidth> const & p) : h(p.h), e(p.e) {}
  };
#else
  template<int index_size = 8>
  struct merge_row_ideal {
      typedef medium_int<index_size> index_type;
      private:
      index_type h;
      public:
      merge_row_ideal() {}
      merge_row_ideal(prime_t p) : h(p.h) {}
      index_type const & index() const { return h; }
      index_type & index() { return h; }
      bool operator<(merge_row_ideal const& a) const { return h < a.h; }
      template<int otherwidth> friend class merge_row_ideal;
      template<int otherwidth>
      merge_row_ideal(merge_row_ideal<otherwidth> const & p) : h(p.h) {}
  };
#endif
  /* }}} */

  typedef unsigned int rowweight_t;
  typedef uint32_t colweight_t;  /* we used to require that it be signed.
                                    It's no longer the case. In fact,
                                    it's the contrary... */

  /* This structure holds the rows. We strive to avoid fragmentation
   * here. */
  typedef compressible_heap<
      merge_row_ideal<compact_column_index_size>,
      rowweight_t,
      merge_row_heap_batch_size> rows_t;

  rows_t rows;

  /* weight per column index. negative measns inactive */
  std::vector<colweight_t> colweights;

  static void declare_usage(param_list_ptr pl);
  bool interpret_parameters(param_list_ptr pl);

  void push_relation(earlyparsed_relation_ptr rel);
  void read_rows(const char *purgedname);
  size_t count_columns_below_weight(size_t *nbm, size_t wmax);
  void bury_heavy_columns();
  void renumber_columns();
  void remove_row(size_t i, size_t keepcolumn = SIZE_MAX);
  void remove_column(size_t j);
  size_t remove_singletons();
  size_t remove_singletons_iterate();

  void merge();

  typedef int32_t pivot_priority_t;
  sparse_indexed_priority_queue<size_t, pivot_priority_t> markowitz_table;

  pivot_priority_t markowitz_count(size_t j) const;

  /* The integer below decides how coarse the allocation should be.
   * Larger means potentially faster reallocs, but also more memory usage
   *
   * Note that this is inherently tied to the colweights vector !
   */
  typedef small_size_pool<size_t, colweight_t, 1> R_pool_t;
  R_pool_t R_pool;
  std::vector<R_pool_t::size_type> R_table;

  void clear_R_table();
  void prepare_R_table(size_t maxlevel);
};

/*{{{ parameters */
void merge_matrix::declare_usage(param_list_ptr pl) {
    param_list_decl_usage(pl, "keep",
            "excess to keep (default "
            STR(DEFAULT_FILTER_EXCESS) ")");
    param_list_decl_usage(pl, "skip",
            "number of heavy columns to bury (default "
            STR(DEFAULT_MERGE_SKIP) ")");
    param_list_decl_usage(pl, "maxlevel",
            "maximum number of rows in a merge " "(default "
            STR(DEFAULT_MERGE_MAXLEVEL) ")");
    param_list_decl_usage(pl, "target_density",
            "stop when the average row density exceeds this value"
            " (default " STR(DEFAULT_MERGE_TARGET_DENSITY) ")");
    param_list_decl_usage(pl, "mkztype",
			  "controls how the weight of a merge is "
			  "approximated (default " STR(DEFAULT_MERGE_MKZTYPE)
			  ")");
    param_list_decl_usage(pl, "wmstmax",
			  "controls until when a mst is used with "
			  "-mkztype 2 (default " STR(DEFAULT_MERGE_WMSTMAX)
			  ")");
}

bool merge_matrix::interpret_parameters(param_list_ptr pl) {
    param_list_parse_size_t(pl, "keep", &keep);
    param_list_parse_size_t(pl, "skip", &nburied);
    param_list_parse_int(pl, "maxlevel", &maxlevel);
    param_list_parse_double(pl, "target_density", &target_density);
    param_list_parse_uint(pl, "mkztype", &mkztype);
    param_list_parse_uint(pl, "wmstmax", &wmstmax);
    if (maxlevel <= 0 || maxlevel > MERGE_LEVEL_MAX) {
        fprintf(stderr,
                "Error: maxlevel should be positive and less than %d\n",
                MERGE_LEVEL_MAX);
        return false;
    }
    if (mkztype > 2) {
	fprintf(stderr, "Error: -mkztype should be 0, 1, or 2.\n");
        return false;
    }
    return true;
}
/*}}}*/

/* Put in nbm[w] for 0 <= w < wmax, the number of ideals of weight w.
 * Return the number of active columns (w > 0) -- this return value is
 * independent of wmax.
 */
size_t merge_matrix::count_columns_below_weight (size_t * nbm, size_t wmax)/*{{{*/
{
  size_t active = 0;
  for (size_t h = 0; h < wmax; nbm[h++] = 0) ;
  for (size_t h = 0; h < ncols; h++) {
      colweight_t w = colweights[h];
      if ((size_t) w < wmax) nbm[w]++;
      active += w > 0;
  }
  return active;
}
/*}}}*/
void merge_matrix::bury_heavy_columns()/*{{{*/
{
    if (!nburied) {
        printf("# No columns were buried.\n");
        weight = total_weight;
        return;
    }

    static const colweight_t colweight_max=std::numeric_limits<colweight_t>::max();

    double tt = seconds();
    std::vector<size_t> heaviest = get_successive_minima(colweights, nburied, std::greater<colweight_t>());

    colweight_t max_buried_weight = colweights[heaviest.front()];
    colweight_t min_buried_weight = colweights[heaviest.back()];

    /* Compute weight of buried part of the matrix. */
    size_t weight_buried = 0;
    for (size_t i = 0; i < heaviest.size() ; i++)
    {
        colweight_t& w = colweights[heaviest[i]];
        weight_buried += w;
        /* since we saturate the weights at 2^31-1, weight_buried might be less
           than the real weight of buried columns, however this can occur only
           when the number of rows exceeds 2^32-1 */
#if DEBUG >= 1
        fprintf(stderr, "# Burying j=%zu (wt = %zu)\n",
                heaviest[i], (size_t) w);
#endif
        /* reset to 0 the weight of the buried columns */
        w = 0;
    }
    printf("# Number of buried columns is %zu (%zu >= weight >= %zu)\n",
            nburied, (size_t) max_buried_weight, (size_t) min_buried_weight);

    bool weight_buried_is_exact = max_buried_weight < colweight_max;
    printf("# Weight of the buried part of the matrix is %s%zu\n",
            weight_buried_is_exact ? "" : ">= ", weight_buried);

    /* Remove buried columns from rows in mat structure */
    printf("# Start to remove buried columns from rels...\n");
    fflush (stdout);
    // // #ifdef HAVE_OPENMP
    // // #pragma omp parallel for
    // // #endif
    for (auto it = rows.begin() ; it != rows.end() ; ++it) {
        merge_row_ideal<compact_column_index_size> * ptr = it->first;
        rowweight_t nl = 0;
        for (rowweight_t i = 0; i < it->second; i++) {
            size_t h = ptr[i].index();
            colweight_t w = colweights[h];
#if 0
            bool bury = w > min_buried_weight;
            bury = bury || (w == min_buried_weight && h > heaviest.back());
            if (!bury)
                ptr[nl++] = ptr[i];
#else
            if (w)
                ptr[nl++] = ptr[i];
#endif
        }
        rows.shrink_value(it, nl);
    }

    printf("# Done. Total bury time: %.1f s\n", seconds()-tt);
    /* compute the matrix weight */
    weight = 0;
    for (auto const & R : rows)
        weight += R.second;

    if (weight_buried_is_exact)
        ASSERT_ALWAYS (weight + weight_buried == total_weight);
    else /* weight_buried is only a lower bound */
        ASSERT_ALWAYS (weight + weight_buried <= total_weight);
}/*}}}*/
/* {{{ read_rows */
void merge_matrix::push_relation(earlyparsed_relation_ptr rel)
{
    merge_row_ideal<> temp[rel->nb];
    std::copy(rel->primes, rel->primes + rel->nb, temp);
    std::sort(temp, temp + rel->nb);
    /* we should also update the column weights */
    for(weight_t i = 0 ; i < rel->nb ; i++) {
        colweight_t & w(colweights[temp[i].index()]);
        rem_ncols += !w;
        w += w != std::numeric_limits<colweight_t>::max();
    }
    total_weight += rel->nb;
    rows.push_back(temp, temp+rel->nb);
}

/* callback function called by filter_rels */
void * merge_matrix_push_relation (void *context_data, earlyparsed_relation_ptr rel)
{
    merge_matrix & M(*(merge_matrix*) context_data);
    M.push_relation(rel);
    return NULL;
}

/* Renumber the non-zero columns to contiguous values [0, 1, 2, ...]
 * We can use this function only for factorization because in this case we do
 * not need the indexes of the columns (contrary to DL where the indexes of the
 * column are printed in the history file). */
void
merge_matrix::renumber_columns()
{
    double tt = seconds();

    printf("# renumbering columns\n");
    std::vector<size_t> p(ncols);
    /* compute the mapping of column indices, and adjust the colweights
     * table */
    size_t h = 0;
    for (size_t j = 0; j < ncols; j++) {
        if (!colweights[j]) continue;
        p[j]=h;
        colweights[h] = colweights[j];
#ifdef TRACE_COL
        if (TRACE_COL == h || TRACE_COL == j)
            printf ("TRACE_COL: column %zu is renumbered into %zu\n", j, h);
#endif
        h++;
    }
    colweights.erase(colweights.begin() + h, colweights.end());
    /* h should be equal to rem_ncols, which is the number of columns with
     * non-zero weight */
    ASSERT_ALWAYS(h + nburied == rem_ncols);

    ncols = h;

    /* apply mapping to the rows. As p is a non decreasing function, the
     * rows are still sorted after this operation. */
    for (auto & R : rows) {
        merge_row_ideal<compact_column_index_size> * ptr = R.first;
        for (rowweight_t i = 0 ; i < R.second ; i++)
            ptr[i].index()=p[ptr[i].index()];
    }
    printf("# renumbering done in %.1f s\n", seconds()-tt);
}

void merge_matrix::read_rows(const char *purgedname)
{
    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline(purgedname, &nrows, &ncols);
    colweights.assign(ncols, colweight_t());

    char *fic[2] = {(char *) purgedname, NULL};

    /* read all rels. */
    rem_nrows = filter_rels(fic,
            (filter_rels_callback_t) &merge_matrix_push_relation, 
            (void*)this,
            EARLYPARSE_NEED_INDEX, NULL, NULL);
    ASSERT_ALWAYS(rem_nrows == nrows);
    weight = total_weight;


    /* print weight count */
    size_t nbm[256];
    size_t active = count_columns_below_weight(nbm, 256);
    for (int h = 1; h <= maxlevel; h++)
        printf ("# There are %zu column(s) of weight %d\n", nbm[h], h);
    printf ("# Total %zu active columns\n", active);
    ASSERT_ALWAYS(rem_ncols == ncols - nbm[0]);

    printf("# Total weight of the matrix: %zu\n", total_weight);

    bury_heavy_columns();

    printf("# Weight of the active part of the matrix: %zu\n", weight);

    /* XXX for DL, we don't want to do this, do we ? */
    renumber_columns();
    
    for (size_t j = 0; j < ncols; j++) {
    }

#if 0
    /* Allocate mat->R[h] if necessary */
    initMatR (mat);

    /* Re-read all rels (in memory) to fill-in mat->R */
    printf("# Start to fill-in columns of the matrix...\n");
    fflush (stdout);
    fillR (mat);
    printf ("# Done\n");
#endif
}
/* }}} */

void merge_matrix::clear_R_table()/*{{{*/
{
    R_table.clear();
    R_pool.clear();
    R_table.assign(colweights.size(), R_pool_t::size_type());
}
/* }}} */

void merge_matrix::prepare_R_table(size_t cwmax)/*{{{*/
{
    double tt = seconds();
    clear_R_table();
    for(size_t j = 0 ; j < colweights.size() ; j++) {
        if (colweights[j] > (colweight_t) cwmax || !colweights[j]) continue;
        R_table[j] = R_pool.alloc(colweights[j]);
    }
    printf("# Time for allocating R: %.2f s\n", seconds()-tt);
    tt=seconds();
    ASSERT_ALWAYS(cwmax <= std::numeric_limits<uint8_t>::max());
    std::vector<uint8_t> pos(colweights.size(), 0);
    for(auto it = rows.begin() ; it != rows.end() ; ++it) {
        merge_row_ideal<compact_column_index_size> * ptr = it->first;
        for (rowweight_t k = 0; k < it->second; k++) {
            size_t h = ptr[k].index();
            if (!R_table[h]) continue;
            colweight_t w = colweights[h];
            R_pool_t::value_type * Rx = R_pool(w, R_table[h]);
            Rx[pos[h]++]=it.i;
        }
    }
    std::stringstream s;
    s << R_pool.printer("# ");
    fputs(s.str().c_str(), stdout);
    printf("# Time for filling R: %.2f s\n", seconds()-tt);
}/*}}}*/

/* This functions returns the difference in matrix element when we add the
 * lightest row with ideal j to all other rows. If ideal j has weight w,
 * and the lightest row has weight w0:
 * * we remove w elements corresponding to ideal j
 * * we remove w0-1 other elements for the lightest row
 * * we add w0-1 elements to the w-1 other rows
 * Thus the result is (w0-1)*(w-2) - w = (w0-2)*(w-2) - 2.
 * (actually we negate this, because we want a high score given that we
 * put all that in a max heap)
 *  
 * Note: we could also take into account the "cancelled ideals", i.e.,
 * the ideals apart the pivot j that got cancelled "by chance". However
 * some experiments show that this does not improve much (if at all) the
 * result of merge.  A possible explanation is that the contribution of
 * cancelled ideals follows a normal distribution, and that on the long
 * run we mainly see the average.
 */
merge_matrix::pivot_priority_t merge_matrix::markowitz_count(size_t j) const
{
    colweight_t w = colweights[j];

    /* empty cols and singletons have max priority */
    if (w <= 1) return std::numeric_limits<pivot_priority_t>::max();
    if (w == 2) return 2;

    /* inactive columns don't count of course, but we shouldn't reach
     * here anyway */
    ASSERT_ALWAYS(R_table[j]);
    if (!R_table[j]) return std::numeric_limits<pivot_priority_t>::min();

    const R_pool_t::value_type * Rx = R_pool(w, R_table[j]);
    rowweight_t w0 = std::numeric_limits<rowweight_t>::max();
    for(colweight_t k = 0 ; k < w ; k++)
        w0 = std::min(w0, rows[Rx[k]].second);

    return 2 - (w0 - 2) * (w - 2);
}

void merge_matrix::remove_row(size_t i, size_t keepcolumn)
{
    rows_t::reference R = rows[i];
    merge_row_ideal<compact_column_index_size> * ptr = R.first;
    for( ; ptr != R.first + R.second ; ptr++) {
        size_t j = ptr->index();
        if (j == keepcolumn) continue;
        if (!R_table[j]) continue;
        /* remove mention of this row in the R entry associated to j */
        colweight_t w = colweights[j];
        R_pool_t::value_type * Rx = R_pool(w, R_table[j]);

        R_pool_t::size_type smaller_key = R_pool.alloc(w-1);
        R_pool_t::value_type * nRx = R_pool(w-1, smaller_key);

        R_pool_t::value_type * me = std::find(Rx, Rx + w, i);
        ASSERT_ALWAYS(me != Rx + w);

        std::copy(Rx, me, nRx);
        std::copy(me + 1, Rx + w, nRx + (me - Rx));

        R_pool.free(w, R_table[j]);
        R_table[j] = smaller_key;
        colweights[j] = w-1;
    }
    rows.kill(i);
}

void merge_matrix::remove_column(size_t j)
{
    colweight_t & w = colweights[j];
    R_pool_t::value_type * Rx = R_pool(w, R_table[j]);
    for(colweight_t k = 0 ; k < w ; k++) {
        remove_row(Rx[k], j);   // do not touch column j now.
    }
    R_pool.free(w, R_table[j]);
    w = 0;
}

size_t merge_matrix::remove_singletons()/*{{{*/
{
    size_t nr = 0;
    for(size_t j = 0 ; j < colweights.size() ; j++) {
        if (colweights[j] == 1) {
            remove_column(j);
            nr++;
        }
    }
    return nr;
}/*}}}*/

size_t merge_matrix::remove_singletons_iterate()/*{{{*/
{
    double tt = seconds();
    size_t tnr = 0;
    for(size_t nr ; (nr = remove_singletons()) ; ) {
        printf("# removed %zu singletons in %.2f s\n", nr, seconds() - tt);
        tnr += nr;
        tt = seconds();
    }
    return tnr;
}
/*}}}*/


void merge_matrix::merge()
{
#if 0
    index_t j;
    index_signed_t mkz;
    uint64_t WN_prev, WN_cur, WN_min;
    double WoverN;
    unsigned int ncost = 0;
    int m;
    merge_stats_t *merge_data;

    printf ("# Using %s to compute the merges\n", __func__);
#endif

    struct compare_by_row_weight {
        merge_matrix const & M;
        compare_by_row_weight(merge_matrix const & M):M(M){}
        bool operator()(size_t a, size_t b) const {
            return M.rows[a].second < M.rows[b].second;
        }
    };
    std::priority_queue<size_t, std::vector<size_t>, compare_by_row_weight> heavy_rows(compare_by_row_weight(*this));

    for(auto const& R : rows)
        heavy_rows.push(R.second);

    for(int cwmax = 2 ; cwmax <= maxlevel ; cwmax++) {
        size_t excess = rem_nrows - rem_ncols;
        printf("# Entering loop for level %d ; %zu rows, %zu cols, %zu excess\n", cwmax, rem_nrows, rem_ncols, excess);

        clear_R_table();

        prepare_R_table(cwmax);

        remove_singletons_iterate();

    }
#if 0
    // clean things
    removeSingletons (rep, mat);

    merge_data = merge_stats_malloc (maxlevel);

    WN_min = WN_cur = compute_WN(mat);
    WoverN = compute_WoverN (mat);

    /* report lines are printed for each multiple of report_incr */
    double report_incr = REPORT_INCR;
    double report_next = ceil (WoverN / report_incr) * report_incr;

    int64_t excess;
    /* stop when excess = keep and (WoverN >= target or queue is empty) */
    while ((excess = mat->rem_nrows - mat->rem_ncols) > mat->keep ||
            (WoverN < target_density || MkzQueueCardinality (mat) > 0))
    {
        if (excess > mat->keep &&
                (WoverN >= target_density || MkzQueueCardinality (mat) == 0))
        { /* we hope removing the remaining excess will decrease the average
             row density and thus enable a few more merges */
            printf ("Removing final excess, nrows=%" PRIu64 "\n",
                    mat->rem_nrows);
            deleteSuperfluousRows (rep, mat, excess - mat->keep, INT_MAX);
            printf ("Removing singletons, nrows=%" PRIu64 "\n", mat->rem_nrows);
            removeSingletons (rep, mat);
            recomputeR (mat);
            WoverN = compute_WoverN (mat);
            continue; /* this cannot loop forever since the excess decreases */
        }

        /* now excess = keep or (WoverN < target_density and #Q > 0) */

        if (excess == mat->keep &&
                (WoverN >= target_density || MkzQueueCardinality (mat) == 0))
            break; /* we cannot do any more merge */

        /* now WoverN < target_density and #Q > 0 */

        /* Do one merge */
        int ret = MkzPopQueue (&j, &mkz, mat);
        ASSERT_ALWAYS(ret != 0);

        m = mat->wt[j];
        ASSERT_ALWAYS(m >= 0);

        /* m=0 can happen if two ideals have weight 2, and are in the same 2
relations: when we merge one of them, the other one will have weight
0 after the merge */
        if (m == 0)
            goto next_merge;

        uint64_t weight0 = mat->weight;
        if (m == 1) /* singleton ideal */
            removeColDefinitely(rep, mat, j);
        else /* m >= 2 */
        {
            fprintf(rep->outfile, "#\n");
            mergeForColumn (rep, mat, m, j);
        }
        index_signed_t real_mkz = mat->weight - weight0;
        if (m > 1 && merge_stats_is_first_merge (merge_data, m))
        {
            fprintf (rep->outfile, "## First %d-merge\n", m);
            print_report (mat);
            printf ("First %d-merge, estimated cost %ld, real cost %ld\n",
                    m, (long) mkz, (long) real_mkz);
            fflush (stdout);
        }

        /* Update values and report if necessary */
        merge_stats_update (merge_data, m, real_mkz);
        WoverN = compute_WoverN (mat);
        WN_prev = WN_cur;
        WN_cur = compute_WN(mat);
        if (WN_cur > WN_prev)
            ncost++;
        else
            ncost = 0;
        if (WN_cur < WN_min)
            WN_min = WN_cur;

        if (WoverN >= report_next)
        {
            int ni2rem;

            print_report (mat);
            report_next += report_incr;
            ni2rem = number_of_superfluous_rows (mat);
            deleteSuperfluousRows (rep, mat, ni2rem, m);
            removeSingletons (rep, mat);
            if (mat->cwmax == mat->mergelevelmax)
                recomputeR (mat);
        }

next_merge:
        /* if we have not yet reached maxlevel:
           (a) if we have reached the target excess, we increase cwmax to maxlevel
           (b) if we have not yet reached keep, and the queue is empty, we increase
           cwmax by 1 only */
        if (mat->cwmax < mat->mergelevelmax &&
                (MkzQueueCardinality (mat) == 0 || excess == mat->keep))
        {
            if (excess == mat->keep)
                mat->cwmax = mat->mergelevelmax;
            else /* MkzQueueCardinality (mat) == 0 */
                mat->cwmax ++;
            recomputeR (mat);
        }
    }

    /* now excess = keep and (WoverN >= target or queue is empty) */

    if (WoverN >= target_density)
        printf ("Reached target density W/N=%.2f.\n", WoverN);
    else
        printf ("Heap is empty, stopping. Rerun with larger maxlevel if "
                "more merges are needed\n");

#if DEBUG >= 1
    checkWeights (mat);
#endif

    merge_stats_print (merge_data, maxlevel);
    free (merge_data);

    printf ("Total number of removed columns: %" PRId64 "\n",
            mat->ncols - mat->rem_ncols);

#ifdef TIMINGS
    for (int m = 2; m < MERGE_LEVEL_MAX; m++)
    {
        if (tfill[m] != 0 || tmst[m] != 0)
            printf ("m=%d: ncalls=%.0f tfill=%.2f tmst=%.2f\n",
                    m, nfill[m], tfill[m], tmst[m]);
        tfill[0] += tfill[m];
        tmst[0] += tmst[m];
        nfill[0] += nfill[m];
    }
    printf ("Total: ncalls=%.0f tfill=%.2f tmst=%.2f\n", nfill[0], tfill[0],
            tmst[0]);
    printf ("Time for recomputeR: %.2f\n", trecomputeR);
#endif
#endif
}

int main(int argc, char *argv[])
{
    char *argv0 = argv[0];

#if 0
    filter_matrix_t mat[1];
    report_t rep[1];
#endif

    int nthreads = 1;
    /* use real MST minimum for wt[j] <= wmstmax */

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    double tt;
    double wct0 = wct_seconds();
    param_list pl;
    int verbose = 0;
    param_list_init(pl);
    declare_usage(pl);
    merge_matrix::declare_usage(pl);
    argv++, argc--;

    param_list_configure_switch(pl, "force-posix-threads",
				&filter_rels_force_posix_threads);
    param_list_configure_switch(pl, "v", &verbose);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    if (argc == 0)
	usage(pl, argv0);

    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv))
	    continue;
	fprintf(stderr, "Unknown option: %s\n", argv[0]);
	usage(pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    const char *purgedname = param_list_lookup_string(pl, "mat");
    const char *outname = param_list_lookup_string(pl, "out");

#if 0
    /* -resume can be useful to continue a merge stopped due  */
    /* to a too small value of -maxlevel                      */
    const char *resumename = param_list_lookup_string(pl, "resume");
    const char *path_antebuffer =
	param_list_lookup_string(pl, "path_antebuffer");
    const char *forbidden_cols =
	param_list_lookup_string(pl, "forbidden-cols");
#endif

    param_list_parse_int(pl, "t", &nthreads);
#ifdef HAVE_OPENMP
    omp_set_num_threads(nthreads);
#endif

    merge_matrix M;

    if (!M.interpret_parameters(pl)) usage(pl, argv0);

    /* Some checks on command line arguments */
    if (param_list_warn_unused(pl)) {
	fprintf(stderr, "Error, unused parameters are given\n");
	usage(pl, argv0);
    }

    if (purgedname == NULL) {
	fprintf(stderr, "Error, missing -mat command line argument\n");
	usage(pl, argv0);
    }
    if (outname == NULL) {
	fprintf(stderr, "Error, missing -out command line argument\n");
	usage(pl, argv0);
    }

    // set_antebuffer_path(argv0, path_antebuffer);

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    /*
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed(outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);
    */

    /* Init structure containing the matrix and the heap of potential merges */
    // initMat(mat, maxlevel, keep, skip);

    /* Read all rels and fill-in the mat structure */
    tt = seconds();

    M.read_rows(purgedname);

    printf("# Time for filter_matrix_read: %2.2lfs\n", seconds() - tt);


    M.merge();


#if 0
    /* resume from given history file if needed */
    if (resumename != NULL)
	resume(rep, mat, resumename);
#endif

#if 0
    /* TODO: this is an untested feature, so not easy to put back in
     * operation */
    /* Some columns can be disable so merge won't use them as pivot */
    if (forbidden_cols != NULL) {
	printf("Disabling columns from %s\n", forbidden_cols);
	matR_disable_cols(mat, forbidden_cols);
    }
#endif

#if 0
    tt = seconds();
    MkzInit(mat, 1);
    printf("Time for MkzInit: %2.2lfs\n", seconds() - tt);

    mergeOneByOne(rep, mat, maxlevel, target_density);

    fclose_maybe_compressed(rep->outfile, outname);
    printf("Final matrix has N=%" PRIu64 " nc=%" PRIu64 " (%" PRId64 ") "
	   "W=%" PRIu64 " W*N=%.2e W/N=%.2f\n",
	   mat->rem_nrows, mat->rem_ncols,
	   ((int64_t) mat->rem_nrows) - ((int64_t) mat->rem_ncols),
	   mat->weight, compute_WN(mat), compute_WoverN(mat));
    fflush(stdout);
    MkzClear(mat, 1);
    clearMat(mat);
#endif

    param_list_clear(pl);

    printf("Total merge time: %.2f seconds\n", seconds());

    print_timing_and_memory(stdout, wct0);

    return 0;
}
