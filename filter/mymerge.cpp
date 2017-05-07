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
#include <queue>
#include <set>

#include "portability.h"

/* We're going to use this transitionally */
#define MKZTYPE_CAVALLAR 0
#define MKZTYPE_PURE 1
#define MKZTYPE_LIGHT 2

#define REPORT_INCR 3.0
#define DEFAULT_EXCESS_INJECT_RATIO 0.04

#include "filter_config.h"
#include "utils_with_io.h"
// #include "merge_replay_matrix.h" /* for filter_matrix_t */
// #include "report.h"     /* for report_t */
// #include "markowitz.h" /* for MkzInit */
// #include "merge_mono.h" /* for mergeOneByOne */
// #include "sparse.h"

#include "select_mpi.h"
#include "medium_int.hpp"
#include "indexed_priority_queue.hpp"
#include "compressible_heap.hpp"
#include "get_successive_minima.hpp"
#define DEBUG_SMALL_SIZE_POOL
#include "small_size_pool.hpp"
#include "minimum_spanning_tree.hpp"

/* NOTE: presently this value has a very significant impact on I/O speed
 * (about a factor of two when reading the purged file), and we have
 * reasons to believe that just about every operation will suffer from
 * the very same problems...
 */
static const int compact_column_index_size = 8;
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
    /* initial values, and aggregated over all nodes */
    size_t initial_nrows=0;
    size_t initial_ncols=0;
    size_t initial_weight=0;  // all non-zero coefficients (not only active)

    /* these are _local_ values */
    size_t nrows=0;     // number of remaining rows
    size_t ncols=0;     // number of remaining columns, including the buried
    size_t weight=0;    // non-zero coefficients in the active part
    /* We also have aggregated versions of these (summed over all nodes in
     * MPI context) */
    size_t global_weight=0;
    size_t global_ncols=0;


    int cwmax=0;
    int maxlevel = DEFAULT_MERGE_MAXLEVEL;
    size_t keep = DEFAULT_FILTER_EXCESS; /* target for initial_nrows-initial_ncols */
    size_t nburied = DEFAULT_MERGE_SKIP;
    unsigned int mkztype = DEFAULT_MERGE_MKZTYPE;
    unsigned int wmstmax = DEFAULT_MERGE_WMSTMAX;
    double target_density = DEFAULT_MERGE_TARGET_DENSITY;
    double excess_inject_ratio = DEFAULT_EXCESS_INJECT_RATIO;

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
            merge_row_ideal operator*(exponent_type const& v) const {
                merge_row_ideal res;
                res.h = h;
                res.e = e * v;
                return res;
            }
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

    /*{{{ parameters */
    static void declare_usage(param_list_ptr pl);
    bool interpret_parameters(param_list_ptr pl);
    /*}}}*/

    /* {{{ row-oriented structures */
    typedef unsigned int row_weight_t;
    /* This structure holds the rows. We strive to avoid fragmentation
     * here. */

    typedef merge_row_ideal<compact_column_index_size> row_value_t;

    typedef compressible_heap<
        row_value_t,
        row_weight_t,
        merge_row_heap_batch_size> rows_t;
    rows_t rows;
    struct row_weight_iterator : public decltype(rows)::iterator {
        typedef decltype(rows)::iterator super;
        typedef std::pair<size_t, row_weight_t> value_type;
        row_weight_iterator(super const & x) : super(x) {}
        value_type operator*() {
            super& s(*this);
            ASSERT_ALWAYS(s->first);
            return std::make_pair(index(), s->second);
        }
    };

    /* Simple and straightforward. Short-lived anyway. Notice that rows_t
     * is not made of shortlived_row_t elements. */
    typedef std::vector<row_value_t> shortlived_row_t;


    void expensive_check();

    void read_rows(const char *purgedname);
    int push_relation(earlyparsed_relation_ptr rel);
    /* }}} */

    /* {{{ column-oriented structures. We have several types of columns:
     *  - "buried" columns. They're essentially deleted from the matrix
     *    forever.
     *  - "inactive" columns, for which R_table[j]==0
     *  - "active" columns, for which col_weights[j] <= cwmax (cwmax loops
     *    from 2 to maxlevel).
     */
    typedef uint32_t col_weight_t;  /* we used to require that it be signed.
                                       It's no longer the case. In fact,
                                       it's even the contrary... */
    /* weight per local column index */
    std::vector<col_weight_t> col_weights;

    /* The integer below decides how coarse the allocation should be.
     * Larger means potentially faster reallocs, but also more memory usage
     *
     * Note that this is inherently tied to the col_weights vector !
     */
    typedef small_size_pool<size_t, col_weight_t, 1> R_pool_t;
    R_pool_t R_pool;
    std::vector<R_pool_t::size_type> R_table;
    void clear_R_table();
    void prepare_R_table();
    size_t count_columns_below_weight(size_t *nbm, size_t wmax);
    void bury_heavy_columns();
    void renumber_columns();
    /* }}} */

    /* {{{ priority queue for active columns to be considered for merging
    */
    typedef int32_t pivot_priority_t;
    struct markowitz_comp_t {
        typedef std::pair<size_t, pivot_priority_t> value_t;
        bool operator()(value_t const& a, value_t const& b) {
            if (a.second < b.second) return true;
            if (b.second < a.second) return false;
            if (a.first < b.first) return true;
            if (b.first < a.first) return false;
            return false;
        }
    };
    typedef sparse_indexed_priority_queue<size_t, pivot_priority_t> markowitz_table_t;
    markowitz_table_t markowitz_table;
    pivot_priority_t markowitz_count(size_t j);
    void fill_markowitz_table();
    /* }}} */

    /*{{{ general high-level operations (local) */
    void addRj(size_t i);
    void subRj(size_t i, size_t keepcolumn = SIZE_MAX);
    void remove_column(size_t j);
    void remove_row_attached(size_t i, size_t keepcolumn = SIZE_MAX);
    void remove_row_detached(size_t i);
    size_t remove_singletons();
    size_t remove_singletons_iterate();
    std::pair<int, std::vector<std::pair<int, int>>> mst_for_column(size_t);
    void pivot_with_column(size_t j);
#if 0
#ifdef FOR_DL
    void multipliers_for_row_combination(
            row_value_t::exponent_type & e0, row_value_t::exponent_type & e1,
            size_t i0, size_t i1, size_t j);
#endif
    /* The set is optional. The routine adds to it the set of column
     * indices (not counting j) which appear more than once in rows i0 and
     * i1).
     */
    row_weight_t weight_of_sum_of_rows(size_t i0, size_t i1, size_t j, std::set<size_t> * S = NULL) const;
#endif
#if 0
    struct compare_by_row_weight {
        merge_matrix const & M;
        compare_by_row_weight(merge_matrix const & M):M(M){}
        bool operator()(size_t a, size_t b) const {
            if (!M.rows[a].first) return true;
            return M.rows[a].second < M.rows[b].second;
        }
    };
#endif
    /* scores in the heavy_rows table are global, and so is the
     * heavy_weight ! */
    high_score_table<size_t, row_weight_t> heavy_rows;
    size_t heavy_weight = 0;
    size_t remove_excess(size_t count = SIZE_MAX);
    void merge();
    /*}}}*/

    /* {{{ reports and statistics */
    uint64_t WN_cur;
    uint64_t WN_min;
    double WoverN;
    /* report lines are printed for each multiple of report_incr */
    double report_incr = REPORT_INCR;
    double report_next;
    void report_init() {
        // probably useless, I think I got the tracking of these right
        // already
        aggregate_weights();
        uint64_t wx = global_weight - heavy_weight;
        WN_cur = (uint64_t) global_ncols * wx;
        WN_min = WN_cur;
        WoverN = (double) wx / (double) nrows;
        report_next = ceil (WoverN / report_incr) * report_incr;
        printf("rank %d report_next = %.1f ; global_ncols=%zu ; wx=%zu\n",
                comm_rank, report_next,
                global_ncols, (size_t) wx);
    }
    double report(bool force = false) {
        aggregate_weights();
        uint64_t wx = global_weight - heavy_weight;
        WN_cur = (uint64_t) global_ncols * wx;
        WN_min = std::min(WN_min, WN_cur);
        WoverN = (double) wx / (double) nrows;

        bool disp = false;
        if ((disp = WoverN >= report_next))
            report_next += report_incr;
        if ((disp || force) && !comm_rank) {
            long vmsize = Memusage();
            long vmrss = Memusage2();
            size_t explained =0;
            explained += rows.allocated_bytes();
            explained += R_pool.allocated_bytes();
            explained += col_weights.capacity()*sizeof(col_weight_t);
            explained += R_table.capacity()*sizeof(R_pool_t::size_type);
            explained += markowitz_table.allocated_bytes();
            explained += heavy_rows.allocated_bytes();

            printf ("N=%zu (%zd) m=%d W=%zu W*N=%.2e W/N=%.2f #Q=%zu [%.1f/%.1f/%.1f]\n",
                    nrows, nrows-global_ncols, cwmax, global_weight,
                    (double) WN_cur, WoverN, 
                    /* Beware: those are local only */
                    markowitz_table.size(),
                    explained/1048576., vmrss/1024.0, vmsize/1024.0);
            printf("done %zu %d-merges, discarded %zu (%.1f%%)\n",
                    done_merges[cwmax], cwmax,
                    discarded_merges[cwmax],
                    100.0 * discarded_merges[cwmax] / (done_merges[cwmax] + discarded_merges[cwmax]));
            /*
               printf("# rows %.1f\n", rows.allocated_bytes() / 1048576.);
               printf("# R %.1f + %.1f + %.1f\n",
               R_pool.allocated_bytes() / 1048576.,
               col_weights.capacity()*sizeof(col_weight_t) / 1048576.,
               R_table.capacity()*sizeof(R_pool_t::size_type) / 1048576.
               );
               printf("# mkz %.1f\n", markowitz_table.allocated_bytes() / 1048576.);
               printf("# heavy %.1f\n", heavy_rows.allocated_bytes() / 1048576.);
               */
            fflush (stdout);
        }
        return WoverN;
    }
    /* }}} */

    /*{{{ MPI-related structures */
    MPI_Comm comm;
    int comm_rank;
    int comm_size;
    void mpi_init(MPI_Comm c) {
        comm = c;
        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);
    }
    bool is_my_col(size_t j) const { return (int) (j % comm_size) == comm_rank; }
    void aggregate_weights() {
        MPI_Allreduce(&weight, &global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    }

    /*}}}*/

    /*{{{ general high-level operations (MPI) */
    void parallel_pivot(std::vector<size_t> const & my_col_indices, size_t n);
    void batch_Rj_update(std::vector<size_t> const & out, std::vector<size_t> const & in);
    void parallel_merge(size_t batch_size);
    size_t collectively_remove_rows(std::vector<size_t> const & killed);
    std::map<int, size_t> discarded_merges;
    std::map<int, size_t> done_merges;
    /*}}}*/
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
    param_list_decl_usage(pl, "excess_inject_ratio",
            "fraction of excess to prune when stepping mergelevel"
            " (default " STR(DEFAULT_EXCESS_INJECT_RATIO) ")");
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
    param_list_parse_double(pl, "excess_inject_ratio", &excess_inject_ratio);
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
    if (excess_inject_ratio < 0 || excess_inject_ratio > 1) {
	fprintf(stderr, "Error: -excess_inject_ratio must be in [0,1]\n");
        return false;
    }
    return true;
}
/*}}}*/

void merge_matrix::expensive_check()
{
    /* First check the column weights */
    std::vector<col_weight_t> check(col_weights.size());
    size_t active = 0;
    size_t nr=0;
    for(size_t i = 0 ; i < rows.size() ; i++) {
        auto const & row(rows[i]);
        int empty = !row.first;
        int all_empty;
        int one_empty;
        MPI_Allreduce(&empty, &one_empty, 1, MPI_INT, MPI_LOR, comm);
        MPI_Allreduce(&empty, &all_empty, 1, MPI_INT, MPI_LAND, comm);
        ASSERT_ALWAYS(one_empty == all_empty);
        if (!row.first) continue;
        nr++;
        active += row.second;
        for(const row_value_t * ptr = row.first ; ptr != row.first + row.second ; ptr++) {
            size_t lj = ptr->index() / comm_size;
            if (!R_table.empty() && R_table[lj]) {
                col_weight_t w = col_weights[lj];
                R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
                ASSERT_ALWAYS(check[lj] < w);
                ASSERT_ALWAYS(Rx[check[lj]] == i);
            }
            check[lj]++;
        }
    }
    ASSERT_ALWAYS(nr == nrows);
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        ASSERT_ALWAYS(check[lj] == col_weights[lj]);
    }
    size_t test_global_weight;
    size_t test_global_ncols;
    MPI_Allreduce(MPI_IN_PLACE, &active, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&weight, &test_global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&ncols, &test_global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    ASSERT_ALWAYS(global_weight == test_global_weight);
    ASSERT_ALWAYS(active == global_weight);
    ASSERT_ALWAYS(global_ncols == test_global_ncols);
    size_t test_nrows = nrows;
    MPI_Bcast(&test_nrows, 1, MPI_MY_SIZE_T, 0, comm);
    ASSERT_ALWAYS(nrows == test_nrows);
}

/* {{{ merge_matrix::push_relation and ::read_rows */
int merge_matrix::push_relation(earlyparsed_relation_ptr rel)
{
    size_t nb = comm_rank ? 0 : rel->nb;
    MPI_Bcast(&nb, 1, MPI_MY_SIZE_T, 0, comm);
    if (nb == SIZE_MAX) return 0;

    merge_row_ideal<> temp[nb];

    if (!comm_rank) {
        std::copy(rel->primes, rel->primes + rel->nb, temp);
        std::sort(temp, temp + rel->nb);
    }
    MPI_Bcast(temp, nb * sizeof(merge_row_ideal<>), MPI_BYTE, 0, comm);

    /* we should also update the column weights */
    weight_t row_weight = 0;
    for(weight_t i = 0 ; i < nb ; i++) {
        if (!is_my_col(temp[i].index())) continue;
        col_weight_t & w(col_weights[temp[i].index() / comm_size]);
        ncols += !w;
        w += w != std::numeric_limits<col_weight_t>::max();
        temp[row_weight++]=temp[i];
    }
    weight += row_weight;
    rows.push_back(temp, row_weight);
    return 1;
}

/* callback function called by filter_rels */
void * merge_matrix_push_relation (void *context_data, earlyparsed_relation_ptr rel)
{
    merge_matrix & M(*(merge_matrix*) context_data);
    M.push_relation(rel);
    return NULL;
}

void merge_matrix::read_rows(const char *purgedname)
{
    /* Read number of rows and cols on first line of purged file */
    if (!comm_rank)
        purgedfile_read_firstline(purgedname, &initial_nrows, &initial_ncols);
    MPI_Bcast(&initial_nrows, 1, MPI_MY_SIZE_T, 0, comm);
    MPI_Bcast(&initial_ncols, 1, MPI_MY_SIZE_T, 0, comm);

    col_weights.assign(iceildiv(initial_ncols, comm_size), col_weight_t());

    /* The way we read the relations is somewhat quirky. Only one node is
     * doing the filter_rels trick, while for all other nodes we'll spin
     * with the internal functions. We append a magic marker at the end,
     * so that both are able to synchronize when the data file has been
     * read completely.
     */
    if (!comm_rank) {
        char *fic[2] = {(char *) purgedname, NULL};
        /* read all rels. */
        nrows = filter_rels(fic,
                (filter_rels_callback_t) &merge_matrix_push_relation, 
                (void*)this,
                EARLYPARSE_NEED_INDEX, NULL, NULL);
        size_t done = SIZE_MAX;
        MPI_Bcast(&done, 1, MPI_MY_SIZE_T, 0, comm);
    } else {
        for( ; push_relation(NULL) ; ) ;
    }
    MPI_Bcast(&nrows, 1, MPI_MY_SIZE_T, 0, comm);
    ASSERT_ALWAYS(rows.size() == nrows);
    printf("# rank %d nrows %zu ncols %zu weight %zu\n",
            comm_rank, nrows, ncols, weight);

    ASSERT_ALWAYS(nrows == initial_nrows);
    MPI_Allreduce(&weight, &initial_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);

    heavy_rows.set_depth(initial_nrows - global_ncols);
    if (comm_size > 1) {
        // create a table of row weights, and allreduce() it.
        std::vector<size_t> all_aggr_row_weights(nrows);
        for(size_t i = 0 ; i < nrows ; i++)
            all_aggr_row_weights[i] = rows[i].second;
        MPI_Allreduce(MPI_IN_PLACE, &all_aggr_row_weights[0], nrows, MPI_MY_SIZE_T, MPI_SUM, comm);
        for(size_t i = 0 ; i < nrows ; i++) {
            heavy_weight += heavy_rows.push(i, all_aggr_row_weights[i]);
            if (!(i % 1000)) {
                printf("rank %d gw %zu hw %zu\n", comm_rank, global_weight, heavy_weight);
            }
        }
        // now everyone has the same heavy_rows pqueue.
    } else {
        for(size_t i = 0 ; i < nrows ; i++) {
            heavy_weight += heavy_rows.push(i, rows[i].second);
        }
    }

    /* print weight count */
    size_t nbm[256];
    size_t active = count_columns_below_weight(nbm, 256);
    size_t temp = ncols + nbm[0];
    printf("# rank %d ncols %zu ncols0 %zu initial_ncols %zu\n",
            comm_rank, ncols, nbm[0], initial_ncols);
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    ASSERT_ALWAYS(initial_ncols == temp);
    for(int rank = 0 ; rank < comm_size ; rank++) {
        MPI_Barrier(comm);
        if (comm_rank != rank) continue;
        for (int h = 1; h <= maxlevel; h++)
            printf ("# Rank %d has %zu column(s) of weight %d\n", comm_rank, nbm[h], h);
        printf ("# Rank %d: total %zu active columns\n", comm_rank, active);
    }
    MPI_Barrier(comm);
    temp = active;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    if (!comm_rank) {
        printf ("# all ranks: total %zu active columns\n", temp);
        printf("# Total weight of the matrix: %zu\n", initial_weight);
    }

    bury_heavy_columns();

    temp = weight;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    if (!comm_rank)
        printf("# Weight of the active part of the matrix: %zu\n", temp);

    /* XXX for DL, we don't want to do this, do we ? */
#ifndef FOR_DL
    renumber_columns();
#endif

    aggregate_weights();
    printf("rank %d gw %zu hw %zu\n", comm_rank, global_weight, heavy_weight);
}
/* }}} */

/*{{{ size_t merge_matrix::count_columns_below_weight */
/* Put in nbm[w] for 0 <= w < wmax, the number of ideals of weight w.
 * Return the number of active columns (w > 0) -- this return value is
 * independent of wmax.
 */
size_t merge_matrix::count_columns_below_weight (size_t * nbm, size_t wmax)
{
  size_t active = 0;
  for (size_t h = 0; h < wmax; nbm[h++] = 0) ;
  for (size_t h = 0; h < initial_ncols; h++) {
      if (!is_my_col(h)) continue;
      col_weight_t w = col_weights[h / comm_size];
      if ((size_t) w < wmax) nbm[w]++;
      active += w > 0;
  }
  return active;
}
/*}}}*/
void merge_matrix::bury_heavy_columns()/*{{{*/
{
    /* make nburied a local value */
    nburied = (nburied / comm_size) + (comm_rank < (int) (nburied % comm_size));

    if (!nburied) {
        if (!comm_rank) printf("# No columns were buried.\n");
        return;
    }

    static const col_weight_t colweight_max=std::numeric_limits<col_weight_t>::max();

    double tt = seconds();
    std::vector<size_t> heaviest = get_successive_minima(col_weights, nburied, std::greater<col_weight_t>());

    col_weight_t max_buried_weight = col_weights[heaviest.front()];
    col_weight_t min_buried_weight = col_weights[heaviest.back()];

    /* Compute weight of buried part of the matrix. */
    size_t weight_buried = 0;
    for (size_t i = 0; i < heaviest.size() ; i++)
    {
        col_weight_t& w = col_weights[heaviest[i]];
        weight_buried += w;
        /* since we saturate the weights at 2^31-1, weight_buried might be less
           than the real weight of buried columns, however this can occur only
           when the number of rows exceeds 2^32-1 */
#if DEBUG >= 1
        fprintf(stderr, "# Burying j=%zu (wt = %zu)\n",
                heaviest[i] * comm_size + comm_rank, (size_t) w);
#endif
        /* reset to 0 the weight of the buried columns */
        w = 0;
    }
    bool weight_buried_is_exact = max_buried_weight < colweight_max;

    for(int rank = 0 ; rank < comm_size ; rank++) {
        MPI_Barrier(comm);
        if (comm_rank != rank) continue;
        printf("# rank %d Number of buried columns is %zu"
                " (%zu >= weight >= %zu)\n",
                rank,
                nburied,
                (size_t) max_buried_weight, (size_t) min_buried_weight);
        printf("# rank %d Weight of the buried part of the matrix is %s%zu\n",
                rank,
                weight_buried_is_exact ? "" : ">= ", weight_buried);
    }
    MPI_Barrier(comm);

    if (!comm_rank)
        printf("# Start to remove buried columns from rels...\n");
    fflush (stdout);

    /* Remove buried columns from rows in mat structure */
// #ifdef HAVE_OPENMP
// #pragma omp parallel for
// #endif
    for (size_t i = 0 ; i < rows.size() ; i++) {
        row_value_t * ptr = rows[i].first;
        if (!ptr) continue;
        row_weight_t nl = 0;
        for (row_weight_t k = 0; k < rows[i].second; k++) {
            size_t h = ptr[k].index();
            col_weight_t w = col_weights[h / comm_size];
#if 0
            bool bury = w > min_buried_weight;
            bury = bury || (w == min_buried_weight && h > heaviest.back());
            if (!bury)
                ptr[nl++] = ptr[k];
#else
            if (w)
                ptr[nl++] = ptr[k];
#endif
        }
        rows.shrink_value_unlocked(i, nl);
    }
    rows.compress();

    if (!comm_rank)
        printf("# Done. Total bury time: %.1f s\n", seconds()-tt);

    size_t weight0 = weight;

    /* compute the matrix weight */
    weight = 0;
    for (auto const & R : rows)
        weight += R.second;

    if (weight_buried_is_exact)
        ASSERT_ALWAYS (weight + weight_buried == weight0);
    else /* weight_buried is only a lower bound */
        ASSERT_ALWAYS (weight + weight_buried <= weight0);
}/*}}}*/
/*{{{ void merge_matrix::renumber_columns() */
/* Renumber the non-zero columns to contiguous values [0, 1, 2, ...]
 * We can use this function only for factorization because in this case we do
 * not need the indexes of the columns (contrary to DL where the indexes of the
 * column are printed in the history file). */
void merge_matrix::renumber_columns()
{
    double tt = seconds();

    printf("# renumbering columns\n");
    std::vector<size_t> p(ncols);
    /* compute the mapping of column indices, and adjust the col_weights
     * table */
    size_t lh = 0;
    for (size_t lj = 0; lj < ncols; lj++) {
        if (!col_weights[lj]) continue;
        p[lj]=lh;
        col_weights[lh] = col_weights[lj];
#ifdef TRACE_COL
        size_t jz = lj * comm_size + comm_rank;
        size_t hz = lh * comm_size + comm_rank;
        if (TRACE_COL == hz || TRACE_COL == jz)
            printf ("TRACE_COL: column %zu is renumbered into %zu\n", jz, hz);
#endif
        lh++;
    }
    col_weights.erase(col_weights.begin() + lh, col_weights.end());
    /* lh should be equal to ncols, which is the number of columns with
     * non-zero weight */
    ASSERT_ALWAYS(lh + nburied == ncols);

    /* XXX We used to set initial_ncols here. Now I no longer have a
     * routine use case for renumber_columns, and I'm skeptical about
     * having initial_ncols deviate from what it's meant to be, namely
     * the header data. I'm commenting this out, waiting to see what
     * breaks downstream. At any rate, just initial_ncols = h cannot
     * work, since initial_ncols is global, and h is global.
    initial_ncols = lh;
     */

    /* On the other hand, it's true that I need to update ncols ! */
    ncols = lh;

    /* apply mapping to the rows. As p is a non decreasing function, the
     * rows are still sorted after this operation. */
    for (auto & R : rows) {
        row_value_t * ptr = R.first;
        for (row_weight_t i = 0 ; i < R.second ; i++) {
            /* make sure we remain in the same congruence class */
            ptr[i].index()=comm_size * p[ptr[i].index() / comm_size] + comm_rank;
        }
    }
    printf("# renumbering done in %.1f s\n", seconds()-tt);
}
/*}}}*/

void merge_matrix::clear_R_table()/*{{{*/
{
    R_table.clear();
    R_pool.clear();
    R_table.assign(col_weights.size(), R_pool_t::size_type());
}
/* }}} */
void merge_matrix::prepare_R_table()/*{{{*/
{
    // double tt = seconds();
    clear_R_table();
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        if (!col_weights[lj]) continue;
        if (col_weights[lj] > (col_weight_t) cwmax) continue;
        R_table[lj] = R_pool.alloc(col_weights[lj]);
    }
    // printf("# Time for allocating R: %.2f s\n", seconds()-tt);
    // tt=seconds();
    ASSERT_ALWAYS(cwmax <= std::numeric_limits<uint8_t>::max());
    std::vector<uint8_t> pos(col_weights.size(), 0);
    for(auto it = rows.begin() ; it != rows.end() ; ++it) {
        row_value_t * ptr = it->first;
        for (row_weight_t k = 0; k < it->second; k++) {
            size_t lj = ptr[k].index() / comm_size;
            if (!R_table[lj]) continue;
            col_weight_t w = col_weights[lj];
            R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
            ASSERT_ALWAYS(pos[lj] < w);
            Rx[pos[lj]++]=it.i;
        }
    }
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        if (!col_weights[lj]) continue;
        if (col_weights[lj] > (col_weight_t) cwmax) continue;
        ASSERT_ALWAYS(pos[lj] == col_weights[lj]);
    }
    // printf("# Time for filling R: %.2f s\n", seconds()-tt);
}/*}}}*/


/*{{{ pivot_priority_t merge_matrix::markowitz_count */
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
merge_matrix::pivot_priority_t merge_matrix::markowitz_count(size_t j)
{
    size_t lj = j / comm_size;
    col_weight_t w = col_weights[lj];

    /* empty cols and singletons have max priority */
    if (w <= 1) return std::numeric_limits<pivot_priority_t>::max();
    if (w == 2) return 2;

    /* inactive columns don't count of course, but we shouldn't reach
     * here anyway */
    ASSERT_ALWAYS(R_table[lj]);
    if (!R_table[lj])
        return std::numeric_limits<pivot_priority_t>::min();

    const R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
    row_weight_t w0 = std::numeric_limits<row_weight_t>::max();
    row_weight_t ws = 0;
    for(col_weight_t k = 0 ; k < w ; k++) {
        w0 = std::min(w0, rows[Rx[k]].second);
        ws += rows[Rx[k]].second;
    }

    merge_matrix::pivot_priority_t prio1 = 2 - (w0 - 2) * (w - 2);
    return prio1;

    /* On a c155, this second approach saves about 0.5% on the number of
     * rows, for a given target density. It is also considerably more
     * expensive computationally speaking (but there is some room for
     * improvement, see mst_for_column().
     *
    merge_matrix::pivot_priority_t prio2 = -mst_for_column(j).first;
    if (prio2 + (merge_matrix::pivot_priority_t) ws >= prio2) prio2 += ws;
    ASSERT_ALWAYS(prio2 >= prio1);
    return prio2;
    */
}
/*}}}*/


void merge_matrix::fill_markowitz_table()/*{{{*/
{
    // double tt = seconds();
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        if (!R_table[lj]) continue;
        size_t j = lj * comm_size + comm_rank;
        markowitz_table.push(std::make_pair(j, markowitz_count(j)));
    }
    // printf("# Time for filling markowitz table: %.2f s\n", seconds()-tt);
}/*}}}*/

void merge_matrix::addRj(size_t i) /*{{{*/
{
    /* add to the R table the entries for row i */
    rows_t::const_reference R = rows[i];
    for(const row_value_t * p = R.first ; p != R.first + R.second ; p++) {
        size_t j = p->index();
        size_t lj = j / comm_size;
        if (!R_table[lj]) { col_weights[lj]++; continue; }
        // printf("# addRj for col %zu\n", j);
        /* add mention of this row in the R entry associated to j */
        col_weight_t w = col_weights[lj];
        R_pool.realloc(col_weights[lj], R_table[lj], w+1);
        R_pool(col_weights[lj], R_table[lj])[w] = i;
        markowitz_table.update(std::make_pair(j, markowitz_count(j)));
    }
}/*}}}*/

void merge_matrix::subRj(size_t i, size_t keepcolumn)/*{{{*/
{
    /* subtract from the R table the entries for row i */
    rows_t::const_reference R = rows[i];
    for(const row_value_t * p = R.first ; p != R.first + R.second ; p++) {
        size_t j = p->index();
        size_t lj = j / comm_size;
        if (j == keepcolumn) continue;
        if (!R_table[lj]) { col_weights[lj]--; continue; }
        // printf("# subRj for col %zu\n", j);
        /* remove mention of this row in the R entry associated to j */
        col_weight_t w = col_weights[lj];
        R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
        R_pool_t::value_type * me = std::find(Rx, Rx + w, i);
        ASSERT_ALWAYS(me < Rx + w);
        std::copy(me + 1, Rx + w, me);
        R_pool.realloc(col_weights[lj], R_table[lj], w-1);
        if (!col_weights[lj]) {
            ncols--;
            markowitz_table.remove(j);
        } else {
            markowitz_table.update(std::make_pair(j, markowitz_count(j)));
        }
    }
}
/*}}}*/

void merge_matrix::remove_row_detached(size_t i)/*{{{*/
{
    ASSERT_ALWAYS(rows[i].first);
    weight -= rows.kill(i);
    std::pair<bool, size_t> rem = heavy_rows.remove(i);
    if (rem.first)
        heavy_weight -= rem.second;
    nrows--;
}/*}}}*/
void merge_matrix::remove_row_attached(size_t i, size_t keepcolumn)/*{{{*/
{
    ASSERT_ALWAYS(rows[i].first);
    subRj(i, keepcolumn);   // do not touch column j now.
    remove_row_detached(i);
}/*}}}*/

void merge_matrix::remove_column(size_t j)/*{{{*/
{
    ASSERT_ALWAYS(is_my_col(j));
    size_t lj = j / comm_size;
    col_weight_t & w = col_weights[lj];
    R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
    R_pool_t::value_type Rx_copy[w];
    std::copy(Rx, Rx+w, Rx_copy);
    for(col_weight_t k = 0 ; k < w ; k++)
        remove_row_attached(Rx_copy[k],j);
    ASSERT_ALWAYS(col_weights[lj] == w);
    R_pool.free(col_weights[lj], R_table[lj]);
    ncols--;
}
/*}}}*/

size_t merge_matrix::remove_singletons()/*{{{*/
{
    // size_t nr = 0;
    std::vector<size_t> killed;
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        /* We might have created more singletons without re-creating the
         * R table entries for these columns.
         */
        if (col_weights[lj] == 1 && R_table[lj]) {
            R_pool_t::value_type * Rx = R_pool(1, R_table[lj]);
            killed.push_back(Rx[0]);
            // remove_column(lj * comm_size + comm_rank);
            // nr++;
        }
    }
    // return nr;
    return collectively_remove_rows(killed);
}/*}}}*/

size_t merge_matrix::remove_singletons_iterate()/*{{{*/
{
    // double tt = seconds();
    size_t tnr = 0;
    for(size_t nr ; (nr = remove_singletons()) ; ) {
        // printf("# removed %zu singletons in %.2f s\n", nr, seconds() - tt);
        tnr += nr;
        // tt = seconds();
    }
    aggregate_weights();
    return tnr;
}
/*}}}*/

struct row_combiner /*{{{*/ {
    typedef merge_matrix::rows_t rows_t;
    typedef merge_matrix::row_value_t row_value_t;
    typedef merge_matrix::row_weight_t row_weight_t;
    typedef merge_matrix::shortlived_row_t row_t;
    row_value_t const * p0;
    row_weight_t n0;
    row_value_t const * p1;
    row_weight_t n1;
    size_t j;
    row_weight_t sumweight;
#ifdef FOR_DL
    row_value_t::exponent_type e0;
    row_value_t::exponent_type e1;
    row_value_t::exponent_type emax0;
    row_value_t::exponent_type emax1;
#endif
    row_combiner(
                row_value_t const * p0, row_weight_t n0,
                row_value_t const * p1, row_weight_t n1, size_t j)
        : p0(p0), n0(n0), p1(p1), n1(n1), j(j)
    {
        ASSERT_ALWAYS(p0);
        ASSERT_ALWAYS(p1);
        sumweight = 0;
#ifdef FOR_DL
        e0 = e1 = emax0 = emax1 = 0;
#endif
    }
    row_combiner(merge_matrix & M, size_t i0, size_t i1, size_t j)
        : row_combiner(
                M.rows[i0].first, M.rows[i0].second,
                M.rows[i1].first, M.rows[i1].second, j) {}
    row_combiner(row_t const& r0, row_t const& r1, size_t j)
        : row_combiner(&r0[0], r0.size(), &r1[0], r1.size(), j) {}
#ifdef FOR_DL
    void multipliers(row_value_t::exponent_type& x0, row_value_t::exponent_type& x1)/*{{{*/
    {
        if (!e0) {
            row_weight_t k0, k1;
            /* We scan backwards because we have better
             * hope of finding j early that way */
            for(k0 = n0, k1 = n1 ; k0 && k1 ; k0--,k1--) {
                if (p0[k0-1].index() < p1[k1-1].index()) {
                    k0++; // we will replay this one.
                } else if (p0[k0-1].index() > p1[k1-1].index()) {
                    k1++; // we will replay this one.
                } else if (p0[k0-1].index() == j && p1[k1-1].index() == j) {
                    e0 = p0[k0-1].exponent();
                    e1 = p1[k1-1].exponent();
                    break;
                }
            }
            ASSERT_ALWAYS(e0 && e1);
            row_value_t::exponent_type d = gcd_int64 ((int64_t) e0, (int64_t) e1);
            e0 /= d;
            e1 /= d;
        }
        x0=e0;
        x1=e1;
    }/*}}}*/
#endif
    row_weight_t weight_of_sum(std::set<size_t> * S = NULL)/*{{{*/
    {
        if (sumweight) return sumweight;
        /* Compute the sum of the (non-buried fragment of) rows i0 and i1,
         * based on the fact that these are expected to be used for
         * cancelling column j.
         */
#ifdef FOR_DL
        /* for discrete log, we need to pay attention to valuations as well */
        multipliers(e0,e1);
#endif
        row_weight_t k0, k1;
        for(k0 = 0, k1 = 0 ; k0 < n0 && k1 < n1 ; k0++,k1++) {
#ifdef FOR_DL
            if (std::abs(p0[k0].exponent()) > emax0)
                emax0=std::abs(p0[k0].exponent());
            if (std::abs(p1[k1].exponent()) > emax1)
                emax1=std::abs(p1[k1].exponent());
#endif
            if (p0[k0].index() < p1[k1].index()) {
                k1--; // we will replay this one.
                sumweight++;
            } else if (p0[k0].index() > p1[k1].index()) {
                k0--; // we will replay this one.
                sumweight++;
            } else if (p0[k0].index() != j) {
                /* Otherwise this column index appears several times */
#ifdef FOR_DL
                /* except for the exceptional case of the ideal we're
                 * cancelling, we have expectations that we'll get a new
                 * coefficient */
                row_value_t::exponent_type e = p1[k1].exponent()*e0-p0[k0].exponent()*e1;
                sumweight += (e != 0);
#endif
                if (S) S->insert(p0[k0].index());
            }
        }
        sumweight += n0 - k0 + n1 - k1;
        return sumweight;
    } /*}}}*/
    int mst_score(std::set<size_t> * S = NULL)
    {
        weight_of_sum(S);
#ifdef FOR_DL
        if (std::abs((int64_t) emax0 * (int64_t) e1) > 32)
            return INT_MAX;
        if (std::abs((int64_t) emax1 * (int64_t) e0) > 32)
            return INT_MAX;
#endif
        return sumweight;
    }


    row_t subtract_naively()/*{{{*/
    {
        /* This creates a new row with R[i0]-R[i1], with mulitplying
         * coefficients adjusted so that column j is cancelled. Rows i0 and
         * i1 are not removed from the pile of rows just yet, this comes at a
         * later step */
        
        /* This also triggers the computation of the exponents if needed.  */
        weight_of_sum();

#ifdef FOR_DL
        multipliers(e0,e1);
        static row_value_t::exponent_type emax = 4;
#endif

        row_weight_t n2 = sumweight;
        row_t res(n2);
        row_value_t * p2 = &(res[0]);

        /* It's of course pretty much the same loop as in other cases. */
        row_weight_t k0 = 0, k1 = 0, k2 = 0;
        for( ; k0 < n0 && k1 < n1 ; k0++,k1++) {
            if (p0[k0].index() < p1[k1].index()) {
                k1--; // we will replay this one.
#ifdef FOR_DL
                p2[k2] = p0[k0] * (-e1);
                    if (std::abs(p2[k2].exponent()) > emax) {
                        emax = std::abs(p2[k2].exponent());
                        printf("# New record coefficient: %ld\n", (long) emax);
                    }
#else
                p2[k2] = p0[k0];
#endif
                k2++;
            } else if (p0[k0].index() > p1[k1].index()) {
                k0--; // we will replay this one.
#ifdef FOR_DL
                p2[k2] = p1[k1] * e0;
                    if (std::abs(p2[k2].exponent()) > emax) {
                        emax = std::abs(p2[k2].exponent());
                        printf("# New record coefficient: %ld\n", (long) emax);
                    }
#else
                p2[k2] = p1[k1];
#endif
                k2++;
            } else if (p0[k0].index() != j) {
                /* Otherwise this column index appears several times */
#ifdef FOR_DL
                /* except for the exceptional case of the ideal we're
                 * cancelling, we have expectations that we'll get a new
                 * coefficient */
                row_value_t::exponent_type e = p1[k1].exponent()*e0-p0[k0].exponent()*e1;
                if (e) {
                    p2[k2].index() = p0[k0].index();
                    p2[k2].exponent() = e;
                    if (std::abs(p2[k2].exponent()) > emax) {
                        emax = std::abs(p2[k2].exponent());
                        printf("# New record coefficient: %ld\n", (long) emax);
                    }
                    k2++;
                }
#endif
            }
        }
        for( ; k0 < n0 ; k0++) {
#ifdef FOR_DL
                p2[k2] = p0[k0] * e0;
                    if (std::abs(p2[k2].exponent()) > emax) {
                        emax = std::abs(p2[k2].exponent());
                        printf("# New record coefficient: %ld\n", (long) emax);
                    }
#else
                p2[k2] = p0[k0];
#endif
                k2++;
        }
        for( ; k1 < n1 ; k1++) {
#ifdef FOR_DL
                p2[k2] = p1[k1] * (-e1);
                    if (std::abs(p2[k2].exponent()) > emax) {
                        emax = std::abs(p2[k2].exponent());
                        printf("# New record coefficient: %ld\n", (long) emax);
                    }
#else
                p2[k2] = p1[k1];
#endif
                k2++;
        }
        ASSERT_ALWAYS(k2 == n2);
        return res;
    }/*}}}*/
}; /*}}}*/

std::pair<int, std::vector<std::pair<int, int>>> merge_matrix::mst_for_column(size_t j)/*{{{*/
{
    std::pair<int, std::vector<std::pair<int, int>>> res;
    size_t lj = j / comm_size;
    weight_t w = col_weights[lj];
    if (w == 0) {
        res.first = 0;
        return res;
    }

    R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
    
    if (w == 1) {
        res.first = -rows[Rx[0]].second;
        return res;
    }

    if (w == 2) {
        res.first = -2;
        res.second.push_back(std::make_pair(0,1));
        return res;
    }

    matrix<int> weights(w,w,INT_MAX);

    /* while we're doing so, we may also count which are the column
     * indices that happen multiple times in these rows */
    // std::set<size_t> S;
    for(col_weight_t k = 0 ; k < w ; k++) {
        for(col_weight_t l = k + 1 ; l < w ; l++) {
            int s = row_combiner(*this, Rx[k], Rx[l], j).mst_score();
            weights(k,l) = s;
            weights(l,k) = s;
        }
    }
    /* Note: our algorithm for computing the weight matrix can be
     * improved a lot. Currently, for a level-w merge (w is the weight of
     * column j), and given rows with average weight L, our cost is
     * O(w^2*L) for creating the matrix.
     *
     * To improve it, first notice that if V_k is the w-dimensional
     * row vector giving coordinates for column k, then the matrix giving at
     * (i0,i1) the valuation of column k when row i1 is cancelled from
     * row i1 is given by trsp(V_k)*V_j-trsp(V_j)*V_k (obviously
     * antisymmetric).
     *
     * We need to do two things.
     *  - walk all w (sorted!) rows simultaneously, maintaining a
     *    priority queue for the "lowest next" column coordinate. By
     *    doing so, we'll be able to react on all column indices which
     *    happen to appear in several rows (because two consecutive
     *    pop()s from the queue will give the same index). So at select
     *    times, this gives a subset of w' rows, with that index
     *    appearing.
     *    This step will cost O(w*L*log(w)).
     *  - Next, whenever we've identified a column index which appears
     *    multiple times (say w' times), build the w'*w' matrix
     *    trsp(V'_k)*V'_j-trsp(V'_j)*V'_k, where all vectors have length
     *    w'. Inspecting the 0 which appear in this matrix gives the
     *    places where we have exceptional cancellations. This can count
     *    as -1 in the total weight matrix.
     * Overall the cost will become O(w*L*log(w)+w'*L'^2+w^2), which
     * should be somewhat faster.
     */

    auto mst = minimum_spanning_tree(weights);

    if (mst.first == INT_MAX) {
#ifdef FOR_DL
        printf("# disconnected graph for %d-merge due to valuations, skipping merge\n", (int) w);
#else
        printf("# Fatal: graph is disconnected\n");
        for(weight_t k = 0 ; k < w ; k++) {
            printf("[%d]", (int) Rx[k]);
            rows_t::const_reference R = rows[Rx[k]];
            for(const row_value_t * p = R.first ; p != R.first + R.second ; p++) {
                if ((uint64_t) p->index() == (size_t) j)
                    printf("\033[01;31m");
                printf(" %" PRId64, (uint64_t) p->index());
#ifdef FOR_DL
                printf(":%d", (int) p->exponent());
#endif
                if ((uint64_t) p->index() == (size_t) j)
                    printf("\033[00m");
            }
            printf("\n");
        }
#endif
    }
    return mst;
}/*}}}*/

void merge_matrix::pivot_with_column(size_t j)/*{{{*/
{
    ASSERT_ALWAYS(comm_size == 1);
    size_t lj = j / comm_size;
    weight_t w = col_weights[lj];
    done_merges[w]++;
    auto mst = mst_for_column(j);

    // now do the merge for real. We need to be cautious so that we don't
    // do an excessive number of reallocs for the R tables. However
    // that's for another implementation. We'll make do with a
    // simple-minded one for the moment.
    // std::map<size_t, std::vector<R_pool_t::size_type> > dup_R;
    

    /* note that a modifying operation may change the pointers in the R
     * table, so we need to grab the list of indices beforehand */
    R_pool_t::value_type * Rx = R_pool(w, R_table[lj]);
    R_pool_t::value_type Rx_copy[w];
    std::copy(Rx, Rx+w, Rx_copy);
    for(auto const & edge : mst.second) {
        row_combiner C(*this, Rx_copy[edge.first], Rx_copy[edge.second], j);
        auto newrow = C.subtract_naively();
        rows.push_back(newrow);
        heavy_weight += heavy_rows.push(rows.size()-1, newrow.size());
        nrows++;
        weight += newrow.size();
        addRj(rows.size()-1);
    }
    remove_column(j);
    /* pivot_with_column is to be used in one-job context only, but since
     * the outer control loop uses the global_ values as well, we must
     * pay attention and keep it in sync.
     */
    aggregate_weights();
}
/*}}}*/

size_t merge_matrix::collectively_remove_rows(std::vector<size_t> const & killed)/*{{{*/
{
    /* each node orders a set of removal of rows, and we remove these
     * rows everywhere. A priori we expect these set to be disjoint, but
     * overlap is tolerated.  return the total number of rows removed (at
     * all nodess, and without duplicates)
     */
    std::vector<int> nkilled(comm_size);
    std::vector<int> displs(comm_size);
    nkilled[comm_rank] = killed.size();
    MPI_Allgather(
            MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            (void*) &nkilled[0], 1, MPI_INT, comm);
    size_t allkilled_size = 0;
    for(int rank = 0 ; rank < comm_size ; rank++) {
        displs[rank] = rank ? displs[rank-1] + nkilled[rank-1] : 0;
        allkilled_size += nkilled[rank];
    }
    std::vector<size_t> allkilled(allkilled_size);
    MPI_Allgatherv(&killed[0], killed.size(), MPI_MY_SIZE_T,
            &allkilled[0], &nkilled[0], &displs[0], MPI_MY_SIZE_T,
            comm);
    size_t nk = 0;
    for(auto i : allkilled) {
        if (rows[i].first) {
            remove_row_attached(i);
            nk++;
        }
    }
    MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&weight, &global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    return nk;
}/*}}}*/

size_t merge_matrix::remove_excess(size_t count)/*{{{*/
{
    size_t excess = nrows - global_ncols;

    count = std::min(count, excess - keep);

    if (!count || heavy_rows.size() <= excess - count)
        return 0;
    
    /*
    printf("# Removing %zd excess rows (from excess = %zu-%zu=%zd)\n",
            count, nrows, ncols, excess);
            */

    size_t nb_heaviest = heavy_rows.size() - (excess - count);
    high_score_table<size_t, row_weight_t> heaviest(nb_heaviest);
    heavy_rows.filter_to(heaviest);
    std::vector<size_t> killed;
    for( ; !heaviest.empty() ; heaviest.pop()) {
        size_t i = heaviest.top().first;
        ASSERT_ALWAYS(rows[i].first);
        killed.push_back(i);
    }

    collectively_remove_rows(killed);


    printf("rank %d removes %zu excess rows\n", comm_rank, killed.size());
    heavy_rows.set_depth(nrows - global_ncols);
    heavy_weight = heavy_rows.sum();
    return nb_heaviest;
}/*}}}*/

/*{{{ collective merge heavylifting */
struct collective_merge_operation {/*{{{*/
    typedef merge_matrix::row_value_t row_value_t;
    typedef merge_matrix::col_weight_t col_weight_t;

    /* Simple and straightforward. Short-lived anyway. */
    typedef merge_matrix::shortlived_row_t row_t;

    /* a row_batch_t has (weight-of-column) rows. */
    typedef std::vector<row_t> row_batch_t;

    /* depending on cases, a row_batch_set_t has:
     * - n keys for the column merges _I_ am requesting,
     * - n*comm_size keys for the "everybody's" variant. */
    typedef std::vector<row_batch_t> row_batch_set_t;

    collective_merge_operation(collective_merge_operation&) = delete;

    MPI_Comm comm;
    MPI_Datatype mpi_row_value_t;
    int comm_rank;
    int comm_size;
    size_t n;   /* batch size */
    size_t w;   /* weight promise. Possibly an upper bound. */

    /* how many of these merge got discarded based on row conflicts ? */
    size_t discarded_merges = 0;
    size_t done_merges = 0;

    /* Because we deal with multi-level structures, we wish to enforce
     * some consistency in naming indices.
     *
     * indices to the list of mpi ranks are named either "p", or "peer".
     * indices to the list of columns in a set of merges considered
     * simultaneously are named "q", and those run up to the bound "n".
     * The values of column indices are named "j".
     * indices within the list of rows that comprise a column are named
     * "r", and run either to the bound "wj", which is specific to j, or
     * "w", which is an upper bound on all wj's.
     * row indices are named "i"
     * indices to row values are named "s", and column indices in rows is
     * "j".
     *
     * Indices to flattened lists (e.g. a list of size n*w) are named by
     * the concatenation of the typical indices for these ranges (e.g.
     * "qr").
     */

    typedef std::vector<size_t> index_vec_t;
    typedef std::vector<int> displ_vec_t;

    index_vec_t my_col_indices;
    // std::vector<size_t> all_col_indices;
    // std::vector<std::vector<size_t>> my_row_indices;
    std::vector<index_vec_t> all_row_indices;

    index_vec_t get_my_row_indices_flat() const {
        index_vec_t res;
        res.reserve(n*w);
        for(size_t q = 0 ; q < n ; q++) {
            index_vec_t const & col = all_row_indices[comm_rank * n + q];
            res.insert(res.end(), col.begin(), col.end());
        }
        return res;
    }

    /* Those are temporary, but nevertheless useful */
    displ_vec_t my_contrib_sizes;
    displ_vec_t remote_contrib_sizes;
    collective_merge_operation(MPI_Comm comm, size_t n, size_t w)
        : comm(comm),
        n(n), w(w)
    {
        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);
        MPI_Type_contiguous(sizeof(row_value_t), MPI_BYTE, &mpi_row_value_t);
        MPI_Type_commit(&mpi_row_value_t);
    }
    ~collective_merge_operation()
    {
        MPI_Type_free(&mpi_row_value_t);
    }

    void announce_columns(std::vector<size_t> const &cols);
    void deduce_and_share_rows(merge_matrix const&);
    void discard_repeated_rows();
    row_batch_set_t merge_scatter_rows(merge_matrix const&);
    /* this function flattens the new set of rows */
    row_batch_t allgather_rows(row_batch_set_t const&);
    void share_contribs();
    struct displ {
        displ_vec_t counts;
        displ_vec_t displs;
        int total;
        int * counts_ptr() { return &counts[0]; }
        int * displs_ptr() { return &displs[0]; }
        const int * counts_ptr() const { return &counts[0]; }
        const int * displs_ptr() const { return &displs[0]; }
    };
    displ send, recv;
    void prepare_sendrecv_displs();
};
/*}}}*/
void collective_merge_operation::announce_columns(std::vector<size_t> const &cols)/*{{{*/
{
    ASSERT_ALWAYS(cols.size() <= n);
    /* w is to be understood as a promise on the column weight, that's all */
    my_col_indices.assign(cols.begin(), cols.end());
    // all_col_indices.reserve(n * comm_size);
    // It seems that we don't need to announce the columns, do we ?

}/*}}}*/
void collective_merge_operation::share_contribs()/*{{{*/
{
    remote_contrib_sizes.assign(comm_size * n * w, 0);
    /* now this needs to be transposed, so that I can now how many
     * coefficients I'll receive from the other folks. */
    MPI_Alltoall(
            (void*) &my_contrib_sizes[0], n * w, MPI_INT,
            (void*) &remote_contrib_sizes[0], n * w, MPI_INT,
            comm);
}/*}}}*/
void collective_merge_operation::prepare_sendrecv_displs()/*{{{*/
{
    send.counts.assign(comm_size, 0);
    send.displs.assign(comm_size, 0);
    send.total = 0;
    for(int peer = 0 ; peer < comm_size ; peer++) {
        int x = 0;
        for(size_t qr = 0 ; qr < n * w ; qr++)
            x += my_contrib_sizes[peer * n * w + qr];
        send.total += x;
        send.counts[peer] = x;
        send.displs[peer] = peer ? send.displs[peer-1] + send.counts[peer-1] : 0;
    }

    recv.counts.assign(comm_size, 0);
    recv.displs.assign(comm_size, 0);
    recv.total = 0;
    for(int peer = 0 ; peer < comm_size ; peer++) {
        int x = 0;
        for(size_t qr = 0 ; qr < n * w ; qr++)
            x += remote_contrib_sizes[peer * n * w + qr]; 
        recv.total += x;
        recv.counts[peer] = x;
        recv.displs[peer] = peer ? recv.displs[peer-1] + recv.counts[peer-1] : 0;
    }

}/*}}}*/
void collective_merge_operation::deduce_and_share_rows(merge_matrix const& M)/*{{{*/
{
    // std::vector<std::vector<size_t>> my_row_indices(n);
    // my_row_indices.assign(n, std::vector<size_t>());
    index_vec_t flat_row_table(comm_size * n * w, SIZE_MAX);
    for(size_t q = 0 ; q < my_col_indices.size() ; q++) {
        size_t j = my_col_indices[q];
        size_t lj = j / comm_size;
        col_weight_t wj = M.col_weights[lj];
        /* notice that wj == 0 is allowed here ! */
        ASSERT_ALWAYS(wj <= w);
        if (!wj) continue;
        ASSERT_ALWAYS(M.R_table[lj]);
        const merge_matrix::R_pool_t::value_type * Rx = M.R_pool(wj, M.R_table[lj]);
        // my_row_indices[q].assign(Rx, Rx + wj);
        std::copy(Rx, Rx + wj, flat_row_table.begin() + (comm_rank * n + q) * w);
    }

    MPI_Allgather(
            MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            (void*) &flat_row_table[0], n * w, MPI_MY_SIZE_T, comm);

    all_row_indices.assign(n * comm_size, index_vec_t());

    for(size_t pq = 0 ; pq < n * comm_size ; pq++) {
        auto ptr0 = flat_row_table.begin() + pq * w;
        auto ptr1 = ptr0 + w;
        for( ; ptr1 > ptr0 && ptr1[-1] == SIZE_MAX ; ptr1--);
        all_row_indices[pq].assign(ptr0, ptr1);
    }
}/*}}}*/
void collective_merge_operation::discard_repeated_rows()/*{{{*/
{
    /* get rid of column merges which involve a row which has already
     * been encountered
     */
    std::set<size_t> met;
    for(std::vector<size_t> & col : all_row_indices) {
        bool ismet = false;
        for(size_t i : col) {
            if (met.find(i) != met.end()) {
                ismet = true;
                break;
            }
        }
        if (ismet) {
            discarded_merges++;
            col.clear();
        } else {
            done_merges += !col.empty();
            for(size_t i : col) met.insert(i);
        }
    }
    /*
    if (!comm_rank && discarded_merges)
        printf("# discarded %zu <=%zu-merges out of %zu\n",
                discarded_merges, w,
                discarded_merges + done_merges);
                */
}/*}}}*/
/* {{{ collective_merge_operation::merge_scatter_rows */
/* This is a collective all-to-all operation. Each node collects from
 * other nodes the remote coefficients from the rows it has an interest
 * in. */
collective_merge_operation::row_batch_set_t
collective_merge_operation::merge_scatter_rows(merge_matrix const& M)
{
    // columns involve an irregular number of rows. However, for good
    // organization of the computation, we want to do as is everyone was
    // dealing with columns of the same weight (and disregard the
    // potentially skipped columns).
    
    ASSERT_ALWAYS(all_row_indices.size() == comm_size * n);

    /* {{{ Get the send-to's and deduce the recv-from's, for each column
     * in the batch, and each involved row */
    my_contrib_sizes.assign(comm_size * n * w, 0);
    for(int peer = 0 ; peer < comm_size ; peer++) {
        for(size_t q = 0 ; q < n ; q++) {
            std::vector<size_t> const& col(all_row_indices[peer * n + q]);
            for(size_t r = 0 ; r < col.size() ; r++)
                my_contrib_sizes[(peer * n + q) * w + r] = M.rows[col[r]].second;
        }
    }
    /* }}} */
    share_contribs();

    prepare_sendrecv_displs();

    /* Put all this stuff in a temp buffer. */
    row_t my_coeffs;
    my_coeffs.reserve(send.total);
    for(auto const & col : all_row_indices) {
        for(size_t i : col) {
            int nr = M.rows[i].second;
            const row_value_t * ptr = M.rows[i].first;
            my_coeffs.insert(my_coeffs.end(), ptr, ptr + nr);
        }
    }
    ASSERT_ALWAYS(my_coeffs.size() == (size_t) send.total);

    row_t remote_coeffs(recv.total);

    MPI_Alltoallv(
            (const void *)&my_coeffs[0],
            send.counts_ptr(), send.displs_ptr(), mpi_row_value_t,
            (void*)&remote_coeffs[0],
            recv.counts_ptr(), recv.displs_ptr(), mpi_row_value_t,
            comm);

    row_batch_set_t res(n);
    const row_value_t * ptrs[comm_size];
    for(int peer = 0 ; peer < comm_size ; peer++)
        ptrs[peer] = &remote_coeffs[recv.displs[peer]];
    for(size_t q = 0 ; q < n ; q++) {
        // all_row_indices has size n * comm_size
        index_vec_t const& col(all_row_indices[comm_rank * n + q]);
        res[q].reserve(col.size());
        for(size_t r = 0 ; r < col.size() ; r++) {
            int x = 0;
            for(int peer = 0 ; peer < comm_size ; peer++) {
                int dx = remote_contrib_sizes[(peer * n + q) * w + r];
                x += dx;
            }
            row_t newrow;
            newrow.reserve(x);

            for(int peer = 0 ; peer < comm_size ; peer++) {
                int dx = remote_contrib_sizes[(peer * n + q) * w + r];
                newrow.insert(newrow.end(), ptrs[peer], ptrs[peer] + dx);
                ptrs[peer] += dx;
            }
            std::sort(newrow.begin(), newrow.end());
            res[q].push_back(std::move(newrow));
        }
    }
    return res;
}
/* }}} */
collective_merge_operation::row_batch_t collective_merge_operation::allgather_rows(row_batch_set_t const& R)/*{{{*/
{
    /* input is _my_ batches */
    ASSERT_ALWAYS(R.size() <= n);

    /* {{{ Get the send-to's and deduce the recv-from's, for each column
     * in the batch, and each involved row */
    /* We'll follow basically the same algorithm as for the
     * merge-scatter. So first we'll count the exchange data with all
     * nodes, and we'll do that via a flat array.
     */
    my_contrib_sizes.assign(comm_size * n * w, 0);
    for(size_t q = 0 ; q < R.size() ; q++) {
        row_batch_t const & batch(R[q]);
        size_t wj = batch.size();
        ASSERT_ALWAYS(wj <= w);
        for(size_t r = 0 ; r < wj ; r++) {
            row_t const & row(batch[r]);
            for(auto jx : row) {
                int peer = jx.index() % comm_size;
                my_contrib_sizes[(peer * n + q) * w + r]++;
            }
        }
    }
    /* }}} */

    // my_contrib_sizes[(peer * n + q) * w + r] =
    // my_contrib_sizes[(peer * n * w) + q * w + r] =
    // the number of coefficients I have in the r-th new row _I_ created,
    // from the q-th merge ordered by me, that will go to this peer.

    share_contribs();

    // remote_contrib_sizes[(p * n + q) * w + r] = the number of new
    // coefficients _I_ will have to store as elements of the r-th new
    // row among the ones in the q-th merge ordered by peer p

    prepare_sendrecv_displs();

    /* Put all this stuff in a temp buffer. */
    row_t my_coeffs(send.total);
    row_value_t * ptrs[comm_size];
    for(int peer = 0 ; peer < comm_size ; peer++)
        ptrs[peer] = &my_coeffs[send.displs[peer]];

    for(size_t q = 0 ; q < R.size() ; q++) {
        row_batch_t const & batch(R[q]);
        size_t wj = batch.size();
        /* In reality, it's wj < w for the use case of a merge, because
         * we make w-1 rows from w rows. */
        ASSERT_ALWAYS(wj <= w);
        for(size_t r = 0 ; r < wj ; r++) {
            row_t const & row(batch[r]);
            for(auto jx : row) {
                int peer = jx.index() % comm_size;
                *(ptrs[peer]++) = jx;
            }
        }
    }

    row_t remote_coeffs(recv.total);

    MPI_Alltoallv(
            (const void *)&my_coeffs[0],
            send.counts_ptr(), send.displs_ptr(), mpi_row_value_t,
            (void*)&remote_coeffs[0],
            recv.counts_ptr(), recv.displs_ptr(), mpi_row_value_t,
            comm);

    row_batch_t res;
    res.reserve(n * comm_size * w);
    for(int peer = 0 ; peer < comm_size ; peer++)
        ptrs[peer] = &remote_coeffs[recv.displs[peer]];

    /* The way we write this loop is where the "virtual ordering"
     * referred to above is actually created */
    for(size_t q = 0 ; q < n ; q++) {
        for(int peer = 0 ; peer < comm_size ; peer++) {
            // my chunks of the q-th merge ordered by peer p. This will
            // entail reading values from peer p only.
            std::vector<size_t> const& col(all_row_indices[peer * n + q]);
            if (col.empty()) continue;
            for(size_t r = 0 ; r < col.size() - 1 ; r++) {
                // remote_contrib_sizes[(peer * n + q) * w + r] = the
                // number of new coefficients _I_ will have to store as
                // elements of the r-th new row among the ones in the
                // q-th merge ordered by peer p
                int dx = remote_contrib_sizes[(peer * n + q) * w + r];
                row_t newrow(ptrs[peer], ptrs[peer] + dx);
                ptrs[peer] += dx;
                std::sort(newrow.begin(), newrow.end());
                res.push_back(std::move(newrow));
            }
        }
    }
    return res;
}/*}}}*/
/*}}}*/

/* {{{ parallel_pivot */
/* There's a virtual ordering here, namely that we do the 0-th merge on
 * node 0, then the 0-th merge on node 1, etc, until the (n-1)-th merge
 * on node (comm_size-1).
 *
 * We arrange so that if either of the two following events occur, a
 * merge is deferred (resulting in no action):
 *  - if a merge involves a row which already appears in another
 *    merge of the same batch.
 *  - if a merge is for a column which appeared in another merge.
 *
 * Note that provided the former event does not occur, the latter does
 * not occur either.
 */
void merge_matrix::parallel_pivot(std::vector<size_t> const & my_col_indices, size_t n)
{
    collective_merge_operation CM(comm, n, cwmax);

    CM.announce_columns(my_col_indices);

    CM.deduce_and_share_rows(*this);

    CM.discard_repeated_rows();
    discarded_merges[cwmax] += CM.discarded_merges;
    done_merges[cwmax] += CM.done_merges;


    /* We'll get a vector of n vectors of up to w rows */
    collective_merge_operation::row_batch_set_t row_batches = CM.merge_scatter_rows(*this);

    /* This has the same type, but with two fine subtleties:
     *  - all_merged_new_rows has n * comm_size fragments (all stored
     *    contiguously), because we receive fragments from everyone !
     *    (my_merged_new_rows has n entries).
     *  - each has wj-1 new rows, not wj.
     */
    collective_merge_operation::row_batch_set_t my_merged_new_rows;

    /* I have rows precisely for the column indices I care about. But
     * these are complete rows, which is good. Let's do the merge, but
     * *without* pushing to the Rj table just yet. We'll do that
     * deferred.
     *
     * Note that this loop is free of MPI calls.
     */
    for(size_t q = 0 ; q < my_col_indices.size() ; q++) {
        size_t j = my_col_indices[q];
        size_t lj = j / comm_size;
        collective_merge_operation::row_batch_t const & batch(row_batches[q]);

        col_weight_t wj = col_weights[lj];
        ASSERT_ALWAYS(wj <= (col_weight_t) cwmax);

        if (batch.empty()) {
            /* this merge was discarded, just skip it.  */
            my_merged_new_rows.push_back(collective_merge_operation::row_batch_t());
            continue;
        }

        ASSERT_ALWAYS(batch.size() == (size_t) wj);

        std::pair<int, std::vector<std::pair<int, int>>> mst;

        if (wj == 1) {
            mst.first = -batch[0].size();
        } else if (wj == 2) {
            mst.first = -2;
            mst.second.push_back(std::make_pair(0,1));
        } else {
            matrix<int> weights(wj,wj,INT_MAX);

            for(col_weight_t k = 0 ; k < wj ; k++) {
                for(col_weight_t l = k + 1 ; l < wj ; l++) {
                    int s = row_combiner(batch[k], batch[l], j).mst_score();
                    weights(k,l) = s;
                    weights(l,k) = s;
                }
            }
            mst = minimum_spanning_tree(weights);
        }

        /* Do the merge, but locally first */
        // R_pool_t::value_type * Rx = R_pool(wj, R_table[j]);
        // R_pool_t::value_type Rx_copy[wj];
        // std::copy(Rx, Rx+wj, Rx_copy);
        collective_merge_operation::row_batch_t new_row_batch;
        for(auto const & edge : mst.second) {
            row_combiner C(batch[edge.first], batch[edge.second], j);
            new_row_batch.push_back(std::move(C.subtract_naively()));
        }
        my_merged_new_rows.push_back(std::move(new_row_batch));
    }

    collective_merge_operation::row_batch_t all_merged_new_rows = CM.allgather_rows(my_merged_new_rows);

    std::vector<size_t> removed_row_indices;
    removed_row_indices.reserve(n * comm_size * cwmax);
    for(auto const & x : CM.all_row_indices)
        removed_row_indices.insert(removed_row_indices.end(), x.begin(), x.end());
    sort(removed_row_indices.begin(), removed_row_indices.end());

    std::vector<size_t> inserted_row_indices;
    inserted_row_indices.reserve(all_merged_new_rows.size());
    /* This adds the rows, but skips the Rj update entirely */
    for(auto const & row : all_merged_new_rows) {
        size_t i = rows.size();
        rows.push_back(row);
        size_t newrow_size = row.size();
        MPI_Allreduce(MPI_IN_PLACE, &newrow_size, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        inserted_row_indices.push_back(i);
        nrows++;
        weight += row.size();
        global_weight += newrow_size;
        heavy_weight += heavy_rows.push(i, newrow_size);
    }
    sort(inserted_row_indices.begin(), inserted_row_indices.end());

    batch_Rj_update(removed_row_indices, inserted_row_indices);

    for(size_t q = 0 ; q < my_col_indices.size() ; q++) {
        size_t j = my_col_indices[q];
        size_t lj = j / comm_size;
        collective_merge_operation::row_batch_t const & batch(row_batches[q]);
        if (!batch.empty()) {
            /* batch_Rj_update had the side-effect of purging the column
             * entirely ! */
            ASSERT_ALWAYS(col_weights[lj] == 0);
            remove_column(j);
        }
    }
    for(auto i : removed_row_indices)
        remove_row_detached(i);
    MPI_Allreduce(&weight, &global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
}
/* }}} */

void merge_matrix::batch_Rj_update(std::vector<size_t> const & out, std::vector<size_t> const & in)/*{{{*/
{
    /* Note: the two arrays must be sorted */
    std::vector<size_t> ii(in.size() + out.size());
    /* We'll maintain a priority queue of size out.size() + in.size(),
     * with the "next index to come" for both the in and out rows. We
     * need to know which of the in/out rows are concerned by a given
     * value. We use the following hack: IN rows are stored as
     * out.size()+x in the priority queue, with 0<=x<in.size().
     */
    typedef std::pair<size_t, size_t> pq_t;
    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> next;
    for(size_t r = 0 ; r < out.size() ; r++) {
        size_t i = out[r];
        if (rows[i].second) {
            next.push(std::make_pair(rows[i].first[0].index(), r));
            ii[r]++;
        }
    }
    for(size_t r = 0 ; r < in.size() ; r++) {
        size_t i = in[r];
        if (rows[i].second) {
            next.push(std::make_pair(rows[i].first[0].index(), r + out.size()));
            ii[r+out.size()]++;
        }
    }
    std::vector<size_t> local;
    for( ; !next.empty() ; ) {
        local.clear();
        pq_t me = next.top();
        /* By construction this loop will run at least once */
        for( ; !next.empty() && next.top().first == me.first ; ) {
            me = next.top();
            size_t r = me.second;
            local.push_back(r);
            next.pop();
            size_t i = r < out.size() ? out[r] : in[r-out.size()];
            if (ii[r] < rows[i].second) {
                next.push(std::make_pair(rows[i].first[ii[r]].index(), r));
                ii[r]++;
            }
        }

        /* We have a j in common, right ? It must be ours */
        size_t j = me.first;
        size_t lj = j / comm_size;
        ASSERT_ALWAYS(is_my_col(j));

        /* Now we have several updates to do to Rj */
        ssize_t sizediff = 0;
        for(auto const & x : local) {
            if (x < out.size()) sizediff--; else sizediff++;
        }

        if (!R_table[lj]) {
            col_weights[lj] += sizediff;
            continue;
        }
        // printf("# addRj for col %zu\n", j);
        /* add mention of this row in the R entry associated to j */
        col_weight_t w0 = col_weights[lj];
        col_weight_t w1 = w0 + sizediff;

        /* We used to *not* do this: this had the effect that
         * higher-level merges were still in the queue. But now that we
         * can do the Rj update after having creating the rows, it seems
         * reasonable to do it.
         */
        if (w0 + sizediff > cwmax) {
            R_pool.free(col_weights[lj], R_table[lj]);
            col_weights[lj] = w1;
            /* don't forget to take it off the priority queue, or we'll
             * have a very bizarre situation */
            markowitz_table.remove(j);
            continue;
        }

        std::vector<size_t> Lout;
        std::vector<size_t> Lin;
        for(auto r : local) {
            if (r < out.size())
                Lout.push_back(out[r]);
            else
                Lin.push_back(in[r-out.size()]);
        }

        /* Ok. It's actually fairly annoying. We can't get away with just
         * a realloc, and in-place modification. To start with, a
         * realloc() is just plain wrong, since the indices we're going
         * to take off need not be at the end. But also, given that we
         * may have both adds and removes simultaneously, we can't say
         * for sure in whih direction it's Right to do the copy.
         *
         * So it's quite unfortunate, but we have no other option than to
         * allocate a new area, and free afterwards.
         */

        R_pool_t::size_type T0 = R_table[lj];
        R_pool_t::size_type T1 = R_pool.alloc(w1);

        const R_pool_t::value_type * Rx0 = R_pool(w0, T0);
        R_pool_t::value_type * Rx1 = R_pool(w1, T1);

        /* In fact, we might be able to do the markowitz count at the
         * same time. Whether it's useful or not depends on the relative
         * size of w0 and sizediff.
         */
        // we rely on the fact that the Rj lists are sorted,
        // and so should be the lists of rows in and out.
        std::vector<size_t>::const_iterator xin = Lin.begin();
        std::vector<size_t>::const_iterator xout = Lout.begin();

        const R_pool_t::value_type * p = Rx0;
        R_pool_t::value_type * q = Rx1;
        for( ; p != Rx0 + w0 && xout != Lout.end();) {
            for( ; p != Rx0 + w0 && *p < *xout ; *q++ = *p++) ;
            if (p == Rx0 + w0) break;
            /* here *p >= *xout */
            if (*p == *xout)
                p++;
            /* here *p > *xout since the row is sorted (or we're past its
             * end) */
            /* No need to do an inner loop, the outer one will do just
             * fine. It can conceivably spin with p unchanged until xout
             * moves to the next one with *p < *xout
             */
            xout++;
        }
        ASSERT_ALWAYS(xout == Lout.end());
        for( ; p != Rx0 + w0 ; *q++ = *p++) ;
        for( ; q != Rx1 + w1 ; *q++ = *xin++) ;
        ASSERT_ALWAYS(xin == Lin.end());

        R_pool.free(w0, T0);
        R_table[lj] = T1;
        col_weights[lj] = w1;

        /* last but not least ! */
        markowitz_table.update(std::make_pair(j, markowitz_count(j)));
    }
}/*}}}*/

void merge_matrix::merge()/*{{{*/
{
    ASSERT_ALWAYS(comm_size == 1);

    report_init();
    int spin = 0;
    report(true);
    for(cwmax = 2 ; cwmax <= maxlevel && spin < 10 ; cwmax++) {

        /* Note: we need the R table for remove_singletons. Singletons
         * which are created within remove_singletons_iterate, but for
         * columns which were originally too heavy to have the R_table
         * entry set, will be discarded in a further pass. */
        prepare_R_table();
        remove_singletons_iterate();

        markowitz_table.clear();
        fill_markowitz_table();

        if (spin && markowitz_table.empty())
            break;

        for( ; report() < target_density && !markowitz_table.empty() ; ) {
            // ASSERT_ALWAYS(markowitz_table.is_heap());

            markowitz_table_t::value_type q = markowitz_table.top();

            /* We're going to remove this column, so it's important that
             * we remove it from the queue right now, otherwise our
             * actions may cause a shuffling of the queue, and we could
             * very well end up popping the wrong one if pop() comes too
             * late.
             */
            markowitz_table.pop();

            /* if the best score is a (cwmax+1)-merge, then we'd better
             * schedule a re-scan of all potential (cwmax+1)-merges
             * before we process this one */
            if (col_weights[q.first / comm_size] > (col_weight_t) cwmax) {
                if (cwmax < maxlevel)
                    break;
                else
                    continue;
            }
            pivot_with_column(q.first);
        }
        // empirical.
        remove_excess((nrows-global_ncols) * excess_inject_ratio);
        remove_singletons_iterate();

        report(true);

        if (cwmax == maxlevel) {
            if (!remove_excess())
                spin++;
            cwmax--;
        }
    }
    report(true);

    printf("# merge stats:\n");
    size_t nmerges=0;
    for(auto x: done_merges) {
        printf("# %d-merges: %zu\n", x.first, x.second);
        nmerges += x.second;
    }
    printf("# Total: %zu merges\n", nmerges);

    /* print weight count */
    size_t nbm[256];
    count_columns_below_weight(nbm, 256);
    for (int h = 1; h <= maxlevel; h++)
        if (nbm[h])
            printf ("# There are %zu column(s) of weight %d\n", nbm[h], h);
    printf("# Total weight of the matrix: %zu\n", initial_weight);
}
/*}}}*/

void merge_matrix::parallel_merge(size_t batch_size)/*{{{*/
{
    report_init();
    int spin = 0;
    report(true);
    for(cwmax = 2 ; cwmax <= maxlevel && spin < 10 ; cwmax++) {

        /* Note: we need the R table for remove_singletons */
        prepare_R_table();
        // expensive_check();

        remove_singletons_iterate();

        markowitz_table.clear();
        fill_markowitz_table();

        int wannabreak = spin && markowitz_table.empty();
        MPI_Allreduce(MPI_IN_PLACE, &wannabreak, 1, MPI_INT, MPI_LAND, comm);
        if (wannabreak) break;

        /* the control structure is global, as report() works on global
         * data */

        for( ; report() < target_density ; ) {
            int empty = markowitz_table.empty();
            MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_LAND, comm);
            if (empty) break;

            std::vector<size_t> my_cols;

            for( ; my_cols.size() < batch_size && !markowitz_table.empty() ; ) {
                markowitz_table_t::value_type q = markowitz_table.top();
                /* We're going to remove this column, so it's important that
                 * we remove it from the queue right now, otherwise our
                 * actions may cause a shuffling of the queue, and we could
                 * very well end up popping the wrong one if pop() comes too
                 * late.
                 */
                markowitz_table.pop();

                /* if the best score is a (cwmax+1)-merge, then we'd better
                 * schedule a re-scan of all potential (cwmax+1)-merges
                 * before we process this one */
                if (col_weights[q.first / comm_size] > (col_weight_t) cwmax) {
                    if (cwmax < maxlevel)
                        markowitz_table.clear();
                } else {
                    my_cols.push_back(q.first);
                }
            }

            // printf("rank %d proposes %zu %d-merges\n", comm_rank, my_cols.size(), cwmax);
            // fflush(stdout);

            parallel_pivot(my_cols, batch_size);
            // expensive_check();
        }
        // empirical.
        remove_excess((nrows-ncols) * excess_inject_ratio);
        remove_singletons_iterate();

        report(true);

        if (cwmax == maxlevel) {
            if (!remove_excess())
                spin++;
            cwmax--;
        }
    }
    report(true);

    printf("# merge stats:\n");
    size_t nmerges=0;
    for(auto x: done_merges) {
        printf("# %d-merges: %zu (discarded %zu)\n",
                x.first, x.second, discarded_merges[x.first]);
        nmerges += x.second;
    }
    printf("# Total: %zu merges\n", nmerges);

    /* print weight count */
    size_t nbm[256];
    count_columns_below_weight(nbm, 256);
    for (int h = 1; h <= maxlevel; h++)
        if (nbm[h])
            printf ("# There are %zu column(s) of weight %d\n", nbm[h], h);
    printf("# Total weight of the matrix: %zu\n", initial_weight);
}
/*}}}*/

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    all_mpi_runtime_checks();
    // mpi_runtime_check<size_t>(MPI_MY_SIZE_T);
    // mpi_runtime_check<uint64_t>(MPI_MY_UINT64_T);
    // mpi_runtime_check<int64_t>(MPI_MY_INT64_T);

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

    M.mpi_init(MPI_COMM_WORLD);

    /* Read all rels and fill-in the mat structure */
    tt = seconds();

    M.read_rows(purgedname);

    printf("# Time for filter_matrix_read: %2.2lfs\n", seconds() - tt);

    if (M.comm_size == 1) {
        M.merge();
    } else {
        M.parallel_merge(128);
    }


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
	   mat->nrows, mat->ncols,
	   ((int64_t) mat->nrows) - ((int64_t) mat->ncols),
	   mat->weight, compute_WN(mat), compute_WoverN(mat));
    fflush(stdout);
    MkzClear(mat, 1);
    clearMat(mat);
#endif

    param_list_clear(pl);

    printf("Total merge time: %.2f seconds\n", seconds());

    print_timing_and_memory(stdout, wct0);

    MPI_Finalize();
    return 0;
}
