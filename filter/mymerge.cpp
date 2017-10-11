#include "cado.h"

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

#define REPORT_INCR 1.0
#define DEFAULT_EXCESS_INJECT_RATIO 0.04

#include "filter_config.h"
#include "utils_with_io.h"

#include "select_mpi.h"
#include "medium_int.hpp"
#include "indexed_priority_queue.hpp"
#include "get_successive_minima.hpp"

/* These are two fancy data structures that, to be honest, do not bring
 * much improvement.
 */
#define xxxR_TABLE_USES_SMALL_SIZE_POOL
#define xxxROW_TABLE_USES_COMPRESSIBLE_HEAP
#define xxxROW_TABLE_USES_SMALL_SIZE_POOL

#if defined(R_TABLE_USES_SMALL_SIZE_POOL) || defined(ROW_TABLE_USES_SMALL_SIZE_POOL)
#define DEBUG_SMALL_SIZE_POOL
#include "small_size_pool.hpp"
#endif

#ifdef ROW_TABLE_USES_COMPRESSIBLE_HEAP
#include "compressible_heap.hpp"
#endif

#include "minimum_spanning_tree.hpp"

/* NOTE: presently this value has a very significant impact on I/O speed
 * (about a factor of two when reading the purged file), and we have
 * reasons to believe that just about every operation will suffer from
 * the very same problems...
 */
static const int compact_column_index_size = 8;


#ifdef ROW_TABLE_USES_COMPRESSIBLE_HEAP
static const int merge_row_heap_batch_size = 16384;
#endif

static void declare_usage(param_list pl)/*{{{*/
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
    param_list_decl_usage(pl, "merge-batch-size", "merge batch size (default is 1 for single-job, 32 above)");
}/*}}}*/

static void usage(param_list pl, char *argv0)/*{{{*/
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}/*}}}*/

struct collective_merge_operation;

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
            template<int otherwidth> friend struct merge_row_ideal;
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
            template<int otherwidth> friend struct merge_row_ideal;
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

#ifdef ROW_TABLE_USES_COMPRESSIBLE_HEAP
    typedef compressible_heap<
        row_value_t,
        row_weight_t,
        merge_row_heap_batch_size> rows_t;
#elif defined(ROW_TABLE_USES_SMALL_SIZE_POOL)
    class rows_t {
        typedef small_size_pool<row_value_t, row_weight_t, 1> row_pool_t;
        row_pool_t row_pool;
        std::vector<row_pool_t::size_type> rows_table;
        std::vector<row_weight_t> weight_table;

        typedef std::pair<row_value_t *, row_weight_t> value_type;

        private:
        rows_t(rows_t const&) = delete;
        rows_t operator=(rows_t const&) = delete;
        size_t allocated_weight = 0;
        public:
        rows_t() =  default;
        ~rows_t() = default;
        size_t allocated_bytes() const {
            return 0; // return container_type::capacity() * sizeof(value_type) + allocated_weight * sizeof(row_value_t) + sizeof(*this);
        }
        size_t overhead_bytes() const {
            return 0; // return allocated_bytes() - container_type::size() * sizeof(value_type) - allocated_weight * sizeof(row_value_t);
        }
        void push_back(const row_value_t * p, row_weight_t n) {
            weight_table.push_back(n);
            size_t idx = row_pool.alloc(n);
            rows_table.push_back(idx);
            std::copy(p, p + n, row_pool(n, idx));
            allocated_weight += n;
        }
        void push_back(std::vector<row_value_t> const & w) {
            push_back(&w[0], w.size());
        }
        void shrink_value_unlocked(size_t i, row_weight_t nl) {
            allocated_weight -= weight_table[i] - nl;
            row_pool.realloc(weight_table[i], rows_table[i], nl);
        }
        size_t kill(size_t i) {
            if (!rows_table[i]) return 0;
            allocated_weight -= weight_table[i];
            row_weight_t oldw = weight_table[i];
            row_pool.free(oldw, rows_table[i]);
            return weight_table[i];
        }
        std::pair<row_value_t *, row_weight_t> operator[](size_t i) {
            return std::make_pair(row_pool(weight_table[i], rows_table[i]), weight_table[i]);
        }
        std::pair<const row_value_t *, row_weight_t> operator[](size_t i) const {
            return std::make_pair(row_pool(weight_table[i], rows_table[i]), weight_table[i]);
        }
        size_t size() const { return rows_table.size(); }
        void compress() {}
    };
#else
    class rows_t : public std::vector<std::pair<row_value_t *, row_weight_t>> {
        typedef std::pair<row_value_t *, row_weight_t> value_type;
        typedef std::vector<value_type> container_type;
        typedef container_type::iterator iterator;
        typedef container_type::const_iterator const_iterator;
        private:
        rows_t(rows_t const&) = delete;
        rows_t operator=(rows_t const&) = delete;
        size_t allocated_weight = 0;
        public:
        rows_t() =  default;
        ~rows_t() { for(auto x : *this) if (x.first) delete[] x.first; }
        size_t allocated_bytes() const {
            return container_type::capacity() * sizeof(value_type) + allocated_weight * sizeof(row_value_t) + sizeof(*this);
        }
        size_t overhead_bytes() const {
            return allocated_bytes() - container_type::size() * sizeof(value_type) - allocated_weight * sizeof(row_value_t);
        }
        void push_back(const row_value_t * p, row_weight_t n) {
            row_value_t * np = new row_value_t[n];
            std::copy(p, p + n, np);
            container_type::push_back(std::make_pair(np, n));
            allocated_weight += n;
        }
        void push_back(std::vector<row_value_t> const & w) {
            push_back(&w[0], w.size());
        }
        void shrink_value_unlocked(size_t i, row_weight_t nl) {
            container_type & v(*this);
            ASSERT_ALWAYS(nl <= v[i].second);
            allocated_weight -= v[i].second - nl;
            row_value_t * p = v[i].first;
            row_value_t * np = new row_value_t[nl];
            std::copy(p, p + nl, np);
            delete[] p;
            v[i].first = np;
            v[i].second = nl;
        }
        size_t kill(size_t i) {
            container_type & v(*this);
            if (!v[i].first) return 0;
            delete[] v[i].first;
            v[i].first = NULL;
            allocated_weight -= v[i].second;
            return v[i].second;
        }
        void compress() {}
    };
#endif
    rows_t rows;

    /* Simple and straightforward. Short-lived anyway. Notice that rows_t
     * is not made of shortlived_row_t elements. */
    typedef std::vector<row_value_t> shortlived_row_t;


    void expensive_check();
    size_t print_active_weight_count();
    void print_final_merge_stats();

    void read_rows(const char *purgedname);
    int push_relation(earlyparsed_relation_ptr rel);
    /* }}} */

    /* {{{ column-oriented structures. We have several types of columns:
     *  - "buried" columns. They're essentially deleted from the matrix
     *    forever.
     *  - "inactive" columns, for which has_indirection(lj)==NULL
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
    typedef std::pair<size_t, row_weight_t> indirection_t;
#ifdef R_TABLE_USES_SMALL_SIZE_POOL
    private:
    typedef small_size_pool<indirection_t, col_weight_t, 1> R_pool_t;
    R_pool_t R_pool;
    std::vector<R_pool_t::size_type> R_table;
    public:
    inline bool has_indirection(size_t lj) const {
        return !R_table.empty() && R_table[lj];
    }
    inline indirection_t * get_indirection(size_t lj) {
        return R_pool(col_weights[lj], R_table[lj]);
    }
    inline const indirection_t * get_indirection(size_t lj) const {
        return R_pool(col_weights[lj], R_table[lj]);
    }
    void allocate_indirection(size_t lj) {
        ASSERT(!R_table[lj]);
        R_table[lj] = R_pool.alloc(col_weights[lj]);
    }
    /* w1 can be used to give the new weight if the reason for
     * deactivation is that the weight has gone over cwmax */
    void deactivate_indirection(size_t lj, size_t w1 = 0) {
        R_pool.free(col_weights[lj], R_table[lj]);
        col_weights[lj] = w1;
    }
    struct temporary_indirection {
        size_t x;
        indirection_t * ptr;
        operator indirection_t*() { return ptr; }
        indirection_t * operator->() { return ptr; }
        indirection_t & operator[](int i) { return ptr[i]; }
    };
    temporary_indirection get_temporary_indirection(size_t w) {
        temporary_indirection t;
        t.x = R_pool.alloc(w);
        t.ptr = R_pool(w, t.x);
        return t;
    }
    void replace_indirection(size_t lj, temporary_indirection T, size_t w) {
        ASSERT_ALWAYS(R_table[lj]);
        R_pool.free(col_weights[lj], R_table[lj]);
        col_weights[lj] = w;
        R_table[lj] = T.x;
    }
#else   /* otherwise we use the libc */
    private:
    std::vector<indirection_t *> R_table;
    public:
    inline bool has_indirection(size_t lj) const {
        return !R_table.empty() && R_table[lj];
    }
    inline indirection_t * get_indirection(size_t lj) {
        return R_table[lj];
    }
    inline const indirection_t * get_indirection(size_t lj) const {
        return R_table[lj];
    }
    void allocate_indirection(size_t lj) {
        ASSERT_ALWAYS(!R_table[lj]);
        R_table[lj] = new indirection_t[col_weights[lj]];
    }
    void deactivate_indirection(size_t lj, size_t w1 = 0) {
        delete[] R_table[lj];
        R_table[lj] = NULL;
        col_weights[lj] = w1;
    }
    typedef indirection_t * temporary_indirection;
    temporary_indirection get_temporary_indirection(size_t w) {
        return new indirection_t[w];
    }
    void replace_indirection(size_t lj, temporary_indirection T, size_t w) {
        ASSERT_ALWAYS(R_table[lj]);
        delete[] R_table[lj];
        R_table[lj] = T;
        col_weights[lj] = w;
    }
#endif

    void prepare_R_table();
    void clear_R_table();

    merge_matrix() = default;
    merge_matrix(merge_matrix const&) = delete;
    ~merge_matrix() { clear_R_table(); }

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
    inline merge_matrix::pivot_priority_t markowitz_count(size_t j, size_t min_rw, size_t cw);
    /* }}} */

    /*{{{ general high-level operations (local) */
    void remove_row_attached(size_t i);
    void remove_row_detached(size_t i);
    size_t remove_singletons();
    size_t remove_singletons_iterate();
    /* scores in the heavy_rows table are global, and so is the
     * heavy_weight ! */
    high_score_table<size_t, row_weight_t> heavy_rows;
    size_t heavy_weight = 0;
    size_t remove_excess(size_t count = SIZE_MAX);
    /*}}}*/

    /* {{{ reports and statistics */
    uint64_t WN_cur;
    uint64_t WN_min;
    double WoverN;
    /* report lines are printed for each multiple of report_incr */
    double report_incr = REPORT_INCR;
    double report_next;
    double t0;
    void report_init() {
        // probably useless, I think I got the tracking of these right
        // already
        aggregate_weights();
        uint64_t wx = global_weight - heavy_weight;
        WN_cur = (uint64_t) global_ncols * wx;
        WN_min = WN_cur;
        WoverN = (double) wx / (double) nrows;
        report_next = ceil (WoverN / report_incr) * report_incr;
        t0 = wct_seconds();
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
        if (disp || force) {
            /* gather some stats counts for everyone before reporting */
            size_t dm = done_merges[cwmax];
            size_t xm = discarded_merges[cwmax];
            size_t mq = markowitz_table.size();
            MPI_Allreduce(MPI_IN_PLACE, &dm, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &xm, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &mq, 1, MPI_MY_SIZE_T, MPI_SUM, comm);

#if 0
            std::map<std::string, size_t> contrib;
            contrib["rows"] = rows.allocated_bytes();
            contrib["rows overhead"] = rows.overhead_bytes();
            contrib["col_weights"] = col_weights.capacity()*sizeof(col_weight_t);
            // contrib["R_pool.alloc"] = R_pool.allocated_bytes();
            // contrib["R_table"] = R_table.capacity()*sizeof(R_pool_t::size_type);
            contrib["markowitz_table"] = markowitz_table.allocated_bytes();
            contrib["heavy_rows"] = heavy_rows.allocated_bytes();

            size_t explained = 0;
            for(auto & x: contrib) {
                MPI_Allreduce(MPI_IN_PLACE, &x.second, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
                explained += x.second;
            }
            size_t vmsize = Memusage();
            size_t vmrss = Memusage2();
            MPI_Allreduce(MPI_IN_PLACE, &vmsize, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &vmrss, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
#endif
            size_t vmrss = Memusage2();
            size_t vmsum = vmrss << 10;
            double vmsum2 = pow(vmrss << 10, 2);
            MPI_Allreduce(MPI_IN_PLACE, &vmsum, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &vmsum2, 1, MPI_DOUBLE, MPI_SUM, comm);
            double vmavg = vmsum / comm_size;
            double vmsdev = sqrt(vmsum2 / comm_size - pow(vmavg, 2));
            vmsdev = abs(100.0 * vmsdev / vmavg);

            if (!comm_rank) {
                char sbuf1[16];
                printf ("%.1f N=%zu (%zd) m=%d W*N=%.2e W/N=%.2f #Q=%zu [%s ~ %.1f%%]\n",
                        wct_seconds() - t0,
                        nrows, nrows-global_ncols, cwmax,
                        // size_disp(global_weight, sbuf0),
                        (double) WN_cur, WoverN, 
                        mq,
                        //explained/1048576., vmrss/1024.0, vmsize/1024.0
                        size_disp(vmavg, sbuf1),
                        vmsdev
                        );
                // printf("# done %zu %d-merges, discarded %zu (%.1f%%)\n", dm, cwmax, xm, 100.0 * xm / dm);
#if 0
                for(auto const& x: contrib) {
                    if (x.first == "rows overhead" || x.second < explained / 10) continue;
                    if (x.first == "rows") {
                        size_t overhead = contrib["rows overhead"];
                        size_t corrected = contrib["rows"] - overhead;
                        printf("# %s: %s+%s [%.1f bytes per active coeff, %.1f without overhead]\n",
                                x.first.c_str(),
                                size_disp(corrected, sbuf0),
                                size_disp(overhead, sbuf1),
                                x.second * 1.0 / global_weight,
                                corrected * 1.0 / global_weight);
                    } else {
                        printf("# %s: %s\n", x.first.c_str(),
                                size_disp(x.second, sbuf0));
                    }
                }
#endif
                /*
                 * printf("# rows %.1f\n", rows.allocated_bytes() / 1048576.);
                 * printf("# R %.1f + %.1f + %.1f\n",
                 * R_pool.allocated_bytes() / 1048576.,
                 * col_weights.capacity()*sizeof(col_weight_t) / 1048576.,
                 * R_table.capacity()*sizeof(R_pool_t::size_type) / 1048576.
                 * );
                 * printf("# mkz %.1f\n", markowitz_table.allocated_bytes() / 1048576.);
                 * printf("# heavy %.1f\n", heavy_rows.allocated_bytes() / 1048576.);
                 * */
            }
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
    void parallel_pivot(collective_merge_operation&);
    void batch_Rj_update(std::vector<size_t> const & out, std::vector<std::pair<size_t, row_weight_t>> const & in = std::vector<std::pair<size_t, row_weight_t>>());
    void parallel_merge(size_t batch_size);
    size_t collectively_remove_rows(std::vector<size_t> const & killed);
    std::map<int, size_t> discarded_merges;
    std::map<int, size_t> done_merges;

    typedef std::vector<std::pair<size_t, pivot_priority_t>> proposal_t;
    void filter_merge_proposal(proposal_t & proposal, proposal_t & trash);
    /*}}}*/
};


void merge_matrix::declare_usage(param_list_ptr pl) { /*{{{ */
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
}

bool merge_matrix::interpret_parameters(param_list_ptr pl) {
    param_list_parse_size_t(pl, "keep", &keep);
    param_list_parse_size_t(pl, "skip", &nburied);
    param_list_parse_int(pl, "maxlevel", &maxlevel);
    param_list_parse_double(pl, "target_density", &target_density);
    param_list_parse_double(pl, "excess_inject_ratio", &excess_inject_ratio);
    if (maxlevel <= 0 || maxlevel > MERGE_LEVEL_MAX) {
        fprintf(stderr,
                "Error: maxlevel should be positive and less than %d\n",
                MERGE_LEVEL_MAX);
        return false;
    }
    if (excess_inject_ratio < 0 || excess_inject_ratio > 1) {
	fprintf(stderr, "Error: -excess_inject_ratio must be in [0,1]\n");
        return false;
    }
    return true;
}
/*}}}*/

void merge_matrix::expensive_check()/*{{{*/
{
    /* The expensive check is *really* expensive. We want the user to
     * notice that it's on, and we want to avoid that he/she gets bored
     * by a program which seemingly produces no output */
    printf("nrows=%zu rank=%d expensive_check\n", nrows, comm_rank);
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
        size_t s = row.second;
        MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        for(const row_value_t * ptr = row.first ; ptr != row.first + row.second ; ptr++) {
            size_t lj = ptr->index() / comm_size;
            if (has_indirection(lj)) {
                col_weight_t w = col_weights[lj];
                indirection_t * Rx = get_indirection(lj);
                ASSERT_ALWAYS(check[lj] < w);
                ASSERT_ALWAYS(Rx[check[lj]].first == i);
                ASSERT_ALWAYS(Rx[check[lj]].second == s);
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
}/*}}}*/

size_t merge_matrix::print_active_weight_count()/*{{{*/
{
    /* print weight count */
    size_t nbm[256];
    size_t active = count_columns_below_weight(nbm, 256);
    size_t temp = ncols + nbm[0];
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    /* This is true only at the very beginning, of course. */
    // ASSERT_ALWAYS(initial_ncols == temp);
    for (int h = 1; h <= maxlevel; h++) {
        size_t n = nbm[h];
        size_t n0 = nbm[h];
        size_t n1 = nbm[h];
        MPI_Allreduce(&n, &n0, 1, MPI_MY_SIZE_T, MPI_MIN, comm);
        MPI_Allreduce(&n, &n1, 1, MPI_MY_SIZE_T, MPI_MAX, comm);
        MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        if (n && !comm_rank)
            printf ("# %zu column(s) [%zu-%zu per node] of weight %d\n", n, n0, n1, h);
    }
    MPI_Allreduce(MPI_IN_PLACE, &active, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    return active;
}/*}}}*/

void merge_matrix::print_final_merge_stats()/*{{{*/
{
    if (!comm_rank) printf("# merge stats:\n");
    size_t nmerges=0;
    for(int w = 0 ; w <= maxlevel ; w++) {
        size_t dm = done_merges[w];
        size_t dm0, dm1, dms;
        size_t xm = discarded_merges[w];
        size_t xm0, xm1, xms;
        MPI_Allreduce(&dm, &dms, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        MPI_Allreduce(&xm, &xms, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        MPI_Allreduce(&dm, &dm0, 1, MPI_MY_SIZE_T, MPI_MIN, comm);
        MPI_Allreduce(&xm, &xm0, 1, MPI_MY_SIZE_T, MPI_MIN, comm);
        MPI_Allreduce(&dm, &dm1, 1, MPI_MY_SIZE_T, MPI_MAX, comm);
        MPI_Allreduce(&xm, &xm1, 1, MPI_MY_SIZE_T, MPI_MAX, comm);
        if (!comm_rank && (dms || xms))
            printf("# %d-merges: %zu [%zu-%zu per node],"
                    " discarded %zu [%zu-%zu per node], avg %.1f%%\n",
                    w, dms, dm0, dm1, xms, xm0, xm1, 100.0 * xms / dms);
        nmerges += dms;
    }
    if (!comm_rank) printf("# Total: %zu merges\n", nmerges);
}/*}}}*/

/* {{{ merge_matrix::push_relation and ::read_rows */
int merge_matrix::push_relation(earlyparsed_relation_ptr rel)
{
    size_t nb = comm_rank ? 0 : rel->nb;
    MPI_Bcast(&nb, 1, MPI_MY_SIZE_T, 0, comm);
    if (nb == SIZE_MAX) return 0;

    std::vector<merge_row_ideal<>> temp(nb);

    if (!comm_rank) {
        std::copy(rel->primes, rel->primes + rel->nb, temp.begin());
        std::sort(temp.begin(), temp.end());
    }
    MPI_Bcast(&temp.front(), nb * sizeof(merge_row_ideal<>), MPI_BYTE, 0, comm);

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
    temp.erase(temp.begin() + row_weight, temp.end());
    rows.push_back(temp);
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
    if (!comm_rank) {
        /* for some reason the prototype purgedfile_read_firstline wants
         * 64-bits. I think it's pointless to even remotely consider
         * dealing with more than 4G entries on 32-bit dinosaurs. */
        uint64_t nr64, nc64;
        purgedfile_read_firstline(purgedname, &nr64, &nc64);
        initial_nrows = nr64;
        initial_ncols = nc64;
        ASSERT_ALWAYS(initial_nrows == nr64);// do we fit ?
        ASSERT_ALWAYS(initial_ncols == nc64);// do we fit ?
    }

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
        for(size_t i = 0 ; i < nrows ; i++)
            heavy_weight += heavy_rows.push(i, all_aggr_row_weights[i]);
        // now everyone has the same heavy_rows pqueue.
    } else {
        for(size_t i = 0 ; i < nrows ; i++) {
            heavy_weight += heavy_rows.push(i, rows[i].second);
        }
    }

    size_t active = print_active_weight_count();

    if (!comm_rank)  {
        printf("# Total %zu active columns\n", active);
        printf("# Total weight of the matrix: %zu\n", initial_weight);
    }

    bury_heavy_columns();

    size_t temp = weight;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    if (!comm_rank)
        printf("# Weight of the active part of the matrix: %zu\n", temp);

    /* XXX for DL, we don't want to do this, do we ? */
#ifndef FOR_DL
    renumber_columns();
#endif

    aggregate_weights();
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
            if (w)
                ptr[nl++] = ptr[k];
        }
        rows.shrink_value_unlocked(i, nl);
    }
    rows.compress();

    if (!comm_rank)
        printf("# Done. Total bury time: %.1f s\n", seconds()-tt);

    size_t weight0 = weight;

    /* compute the matrix weight */
    weight = 0;
    for (size_t i = 0 ; i < rows.size() ; i++)
        if (rows[i].first) weight += rows[i].second;

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
    for (size_t i = 0 ; i < rows.size() ; i++) {
        auto R(rows[i]);
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
    for(size_t lj = 0 ; lj < R_table.size() ; lj++)
        deactivate_indirection(lj, col_weights[lj]);
}
/*}}}*/
void merge_matrix::prepare_R_table()/*{{{*/
{
#ifdef R_TABLE_USES_SMALL_SIZE_POOL
    R_table.assign(col_weights.size(), 0);
#else
    if (R_table.empty())
        R_table.assign(col_weights.size(), NULL);
#endif
    const uint8_t max8 = std::numeric_limits<uint8_t>::max();
    ASSERT_ALWAYS(cwmax <= max8);
    std::vector<uint8_t> pos(col_weights.size(), 0);
    std::vector<uint8_t> minpos(col_weights.size(), 0);
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        if (!col_weights[lj]) continue;
        if (col_weights[lj] > (col_weight_t) cwmax)
            continue;
        if (has_indirection(lj))
            pos[lj] = max8;
        else
            allocate_indirection(lj);
    }

    /* We need a window of known row weights. Because it's somewhat
     * inefficient to to the Allreduce for each new row, we batch it
     * a little.
     */
    std::vector<size_t> rw;
    rw.reserve(1024);
    size_t rwok = 0;
    for(size_t i = 0 ; i < rows.size() ; i++) {
        row_value_t * ptr = rows[i].first;
        if (rwok >= rw.size()) {
            rw.assign(1024, 0);
            rwok = 0;
            for(size_t l = 0 ; rw.size() < rw.capacity() && (i + l < rows.size()) ; l++) {
                if (!rows[i+l].first) continue;
                rw[l] = rows[i+l].second;
            }
            MPI_Allreduce(MPI_IN_PLACE, &rw[0], rw.size(), MPI_MY_SIZE_T, MPI_SUM, comm);
        }
        ASSERT_ALWAYS(rwok < rw.size());
        row_weight_t this_rw = rw[rwok++];
        if (!ptr) continue;
        for (row_weight_t k = 0; k < rows[i].second; k++) {
            size_t lj = ptr[k].index() / comm_size;
            if (!R_table[lj]) continue;
            if (pos[lj] == max8) continue;
            col_weight_t w = col_weights[lj];
            ASSERT_ALWAYS(pos[lj] < w);
            indirection_t * Rx = get_indirection(lj);
            Rx[pos[lj]]=std::make_pair(i, this_rw);
            if (this_rw < Rx[minpos[lj]].second)
                minpos[lj]=pos[lj];
            pos[lj]++;
        }
    }
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        col_weight_t w = col_weights[lj];
        if (!w) continue;
        if (w > (col_weight_t) cwmax) continue;
        if (pos[lj] == max8) continue;
        ASSERT_ALWAYS(pos[lj] == w);
        indirection_t * Rx = get_indirection(lj);
        row_weight_t w0 = Rx[minpos[lj]].second;
        size_t j = comm_size * lj + comm_rank;
        markowitz_table.push(std::make_pair(j, markowitz_count(j, w0, w)));
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
 *
 *
 * The approach of computing the weight matrix explicitly, and
 * prioritizing the merges based on that, does not seem to pay (on a
 * c155, we barely win 0.5%).
 */
inline merge_matrix::pivot_priority_t merge_matrix::markowitz_count(size_t j MAYBE_UNUSED, size_t min_rw, size_t cw)
{
    return 2 - (min_rw - 2) * (cw - 2);
}
/*}}}*/


void merge_matrix::remove_row_detached(size_t i)/*{{{*/
{
    /* note that it is conceivable that this row is still alive, but with
     * an empty local weight */
    ASSERT_ALWAYS(rows[i].first || !rows[i].second);
    weight -= rows.kill(i);
    std::pair<bool, size_t> rem = heavy_rows.remove(i);
    if (rem.first)
        heavy_weight -= rem.second;
    nrows--;
}/*}}}*/

void merge_matrix::remove_row_attached(size_t i)/*{{{*/
{
    ASSERT_ALWAYS(rows[i].first || !rows[i].second);
    batch_Rj_update(std::vector<size_t>(1,i));
    remove_row_detached(i);
}/*}}}*/

size_t merge_matrix::remove_singletons()/*{{{*/
{
    // size_t nr = 0;
    std::vector<size_t> killed;
    for(size_t lj = 0 ; lj < col_weights.size() ; lj++) {
        /* We might have created more singletons without re-creating the
         * R table entries for these columns.
         */
        if (col_weights[lj] == 1 && has_indirection(lj)) {
            indirection_t * Rx = get_indirection(lj);
            killed.push_back(Rx[0].first);
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
    // recall formula: (coeff from row 1)*e0-(coeff from row 0)*e1;
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
        // recall formula: (coeff from row 1)*e0-(coeff from row 0)*e1;
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
        if (std::abs((int64_t) emax0 * (int64_t) e1) + std::abs((int64_t) emax1 * (int64_t) e0) > 32)
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
        struct store_record {
            row_value_t::exponent_type emax = 4;
            bool operator()(row_value_t::exponent_type e) {
                bool b = std::abs(e) > emax;
                if (b) {
                    emax = std::abs(e);
                    printf("# New record coefficient: %ld\n", (long) e);
                }
                return b;
            }
        };
        static store_record check_coeff;
#endif

        row_weight_t n2 = sumweight;
        row_t res(n2);
        row_value_t * p2 = &(res[0]);

        /* It's of course pretty much the same loop as in other cases. */
        // recall formula: (coeff from row 1)*e0-(coeff from row 0)*e1;
        row_weight_t k0 = 0, k1 = 0, k2 = 0;
        for( ; k0 < n0 && k1 < n1 ; k0++,k1++) {
            if (p0[k0].index() < p1[k1].index()) {
                k1--; // we will replay this one.
#ifdef FOR_DL
                p2[k2] = p0[k0] * (-e1);
                check_coeff(p2[k2].exponent());
#else
                p2[k2] = p0[k0];
#endif
                k2++;
            } else if (p0[k0].index() > p1[k1].index()) {
                k0--; // we will replay this one.
#ifdef FOR_DL
                p2[k2] = p1[k1] * e0;
                check_coeff(p2[k2].exponent());
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
                    check_coeff(p2[k2].exponent());
                    k2++;
                }
#endif
            }
        }
        for( ; k0 < n0 ; k0++) {
#ifdef FOR_DL
                p2[k2] = p0[k0] * (-e1);
                check_coeff(p2[k2].exponent());
#else
                p2[k2] = p0[k0];
#endif
                k2++;
        }
        for( ; k1 < n1 ; k1++) {
#ifdef FOR_DL
                p2[k2] = p1[k1] * e0;
                check_coeff(p2[k2].exponent());
#else
                p2[k2] = p1[k1];
#endif
                k2++;
        }
        ASSERT_ALWAYS(k2 == n2);
        return res;
    }/*}}}*/
}; /*}}}*/

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
    std::set<size_t> allkilled_set(allkilled.begin(), allkilled.end());
    for(auto i : allkilled_set) {
        remove_row_attached(i);
    }
    MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&weight, &global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    return allkilled.size();
}/*}}}*/

size_t merge_matrix::remove_excess(size_t count)/*{{{*/
{
    size_t excess = nrows - global_ncols;

    ASSERT_ALWAYS(heavy_rows.size() <= excess);

    /* First count how many rows we would remove if we had
     * heavy_rows.size() == excess. That is an upper bound */
    count = std::min(count, excess - keep);

    /* Then take into account the fact that some may be missing (because
     * they might have been removed, being identified as singletons).
     */
    if (count <= (excess - heavy_rows.size()))
        return 0;
    
    /*
    printf("# Removing %zd excess rows (from excess = %zu-%zu=%zd)\n",
            count, nrows, ncols, excess);
            */

    size_t nb_heaviest = count - (excess - heavy_rows.size());
    high_score_table<size_t, row_weight_t> heaviest(nb_heaviest);
    heavy_rows.filter_to(heaviest);
    std::vector<size_t> killed;
    for( ; !heaviest.empty() ; heaviest.pop()) {
        size_t i = heaviest.top().first;
        ASSERT_ALWAYS(rows[i].first);
        killed.push_back(i);
    }

    collectively_remove_rows(killed);

    // printf("rank %d removes %zu excess rows\n", comm_rank, killed.size());
    heavy_rows.set_depth(nrows - global_ncols);
    heavy_weight = heavy_rows.sum();
    return nb_heaviest;
}/*}}}*/

/*{{{ collective merge heavylifting */

struct displ { /*{{{ helper for AlltoAllv */
    typedef std::vector<int> displ_vec_t;
    displ_vec_t counts;
    displ_vec_t displs;
    int total;
    int * counts_ptr() { return &counts[0]; }
    int * displs_ptr() { return &displs[0]; }
    const int * counts_ptr() const { return &counts[0]; }
    const int * displs_ptr() const { return &displs[0]; }
    displ(displ_vec_t & local, size_t spread) {
        ASSERT_ALWAYS(local.size() % spread == 0);
        size_t n = local.size() / spread;
        counts.assign(n, 0);
        displs.assign(n, 0);
        total = 0;
        for(size_t peer = 0 ; peer < n ; peer++) {
            int x = 0;
            for(size_t qr = 0 ; qr < spread ; qr++)
                x += local[peer * spread + qr];
            displs[peer] = total;
            total += counts[peer] = x;
        }
    }
};/*}}}*/


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
    size_t my_discarded_merges = 0;
    size_t my_done_merges = 0;

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

    index_vec_t my_col_indices;

    typedef std::vector<merge_matrix::markowitz_table_t::value_type> proposal_t;
    proposal_t proposal;
    
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

    void announce_columns(proposal_t const &prop);
    void deduce_and_share_rows(merge_matrix const&);
    proposal_t discard_repeated_rows();
    row_batch_set_t merge_scatter_rows(merge_matrix const&);
    /* this function flattens the new set of rows */
    row_batch_t allgather_rows(row_batch_set_t const&);
};
/*}}}*/
void collective_merge_operation::announce_columns(proposal_t const &prop)
{
    ASSERT_ALWAYS(prop.size() <= n);
    /* w is to be understood as a promise on the column weight, that's all */
    proposal = prop;
    my_col_indices.reserve(prop.size());
    for(auto const & x : proposal)
        my_col_indices.push_back(x.first);
    // all_col_indices.reserve(n * comm_size);
    // It seems that we don't need to announce the columns, do we ?

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
        const merge_matrix::indirection_t * Rx = M.get_indirection(lj);
        // my_row_indices[q].assign(Rx, Rx + wj);
        for(size_t s = 0 ; s < (size_t) wj ; s++)
            flat_row_table[(comm_rank * n + q) * w + s] = Rx[s].first;
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
collective_merge_operation::proposal_t collective_merge_operation::discard_repeated_rows()/*{{{*/
{
    /* get rid of column merges which involve a row which has already
     * been encountered
     */
    std::set<size_t> met;
    proposal_t requeue;
    for(int peer = 0 ; peer < comm_size ; peer++) {
        for(size_t q = 0 ; q < n ; q++) {
            std::vector<size_t> & col(all_row_indices[peer * n + q]);
            bool ismet = false;
            for(size_t i : col) {
                if (met.find(i) != met.end()) {
                    ismet = true;
                    break;
                }
            }
            if (ismet) {
                if (peer == comm_rank) {
                    requeue.push_back(proposal[q]);
                    my_discarded_merges++;
                }
                col.clear();
            } else {
                if (peer == comm_rank)
                    my_done_merges += !col.empty();
                for(size_t i : col) met.insert(i);
            }
        }
    }
    ASSERT_ALWAYS(requeue.size() == my_discarded_merges);
    return requeue;
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

    /* {{{ Prepare the send buffer */
    displ::displ_vec_t my_contrib_sizes(comm_size * n * w, 0);
    for(size_t pq = 0 ; pq < n * comm_size ; pq++) {
        std::vector<size_t> const& col(all_row_indices[pq]);
        for(size_t r = 0 ; r < col.size() ; r++) {
            size_t i = col[r];
            my_contrib_sizes[pq * w + r] = M.rows[i].second;
        }
    }
    /* Put all this stuff in a temp buffer. */
    displ send(my_contrib_sizes, n * w);
    row_t my_coeffs;
    my_coeffs.reserve(send.total);
    for(size_t pq = 0 ; pq < n * comm_size ; pq++) {
        std::vector<size_t> const& col(all_row_indices[pq]);
        for(size_t r = 0 ; r < col.size() ; r++) {
            size_t i = col[r];
            const row_value_t * ptr = M.rows[i].first;
            size_t nr = M.rows[i].second;
            my_coeffs.insert(my_coeffs.end(), ptr, ptr + nr);
        }
    }
    ASSERT_ALWAYS(my_coeffs.size() == (size_t) send.total);
    /* }}} */


    displ::displ_vec_t remote_contrib_sizes(comm_size * n * w, 0);
    MPI_Alltoall(
            (void*) &my_contrib_sizes[0], n * w, MPI_INT,
            (void*) &remote_contrib_sizes[0], n * w, MPI_INT,
            comm);

    displ recv(remote_contrib_sizes, n * w);
    row_t remote_coeffs(recv.total);

    MPI_Alltoallv(
            (const void *)&my_coeffs[0],
            send.counts_ptr(), send.displs_ptr(), mpi_row_value_t,
            (void*)&remote_coeffs[0],
            recv.counts_ptr(), recv.displs_ptr(), mpi_row_value_t,
            comm);

    /* we receive locally. So only for n batches of w rows. Note though
     * that for each of these batches, we may receive from all peers */
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

    /* {{{ Prepare the send buffer */
    /* We'll follow basically the same algorithm as for the
     * merge-scatter. So first we'll count the exchange data with all
     * nodes, and we'll do that via a flat array.
     */
    displ::displ_vec_t my_contrib_sizes(comm_size * n * w, 0);
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
                my_contrib_sizes[(peer * n + q) * w + r]++;
            }
        }
    }
    displ send(my_contrib_sizes, n * w);
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
    /* }}} */

    displ::displ_vec_t remote_contrib_sizes(comm_size * n * w, 0);
    MPI_Alltoall(
            (void*) &my_contrib_sizes[0], n * w, MPI_INT,
            (void*) &remote_contrib_sizes[0], n * w, MPI_INT,
            comm);
    displ recv(remote_contrib_sizes, n * w);
    row_t remote_coeffs(recv.total);

    // my_contrib_sizes[(peer * n + q) * w + r] =
    // my_contrib_sizes[(peer * n * w) + q * w + r] =
    // the number of coefficients I have in the r-th new row _I_ created,
    // from the q-th merge ordered by me, that will go to this peer.

    // remote_contrib_sizes[(p * n + q) * w + r] = the number of new
    // coefficients _I_ will have to store as elements of the r-th new
    // row among the ones in the q-th merge ordered by peer p

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
    for(size_t pq = 0 ; pq < n * comm_size ; pq++) {
        // my chunks of the q-th merge ordered by peer p. This will
        // entail reading values from peer p only.
        std::vector<size_t> const& col(all_row_indices[pq]);
        if (col.empty()) continue;
        for(size_t r = 0 ; r < col.size() - 1 ; r++) {
            // remote_contrib_sizes[(peer * n + q) * w + r] = the
            // number of new coefficients _I_ will have to store as
            // elements of the r-th new row among the ones in the
            // q-th merge ordered by peer p
            int dx = remote_contrib_sizes[pq * w + r];
            int peer = pq/n;
            row_t newrow(ptrs[peer], ptrs[peer] + dx);
            ptrs[peer] += dx;
            std::sort(newrow.begin(), newrow.end());
            res.push_back(std::move(newrow));
        }
    }
    return res;
}/*}}}*/
/*}}}*/

matrix<int> batch_compute_weights(collective_merge_operation::row_batch_t const & batch, size_t j)/*{{{*/
{
    /* compute the weight matrix for cancelling j with the rows given in
     * argument */

    /* 
     * The naive algorithm, for a level-w merge (w is the
     * weight of column j), and given rows with average weight L,
     * costs O(w^2*L) for creating the weight matrix.
     *
     * For the improved version here, first notice that if V_k is the
     * w-dimensional row vector giving coordinates for column k, then the
     * matrix giving at (i0,i1) the valuation of column k when row i1 is
     * cancelled from row i1 is given by trsp(V_k)*V_j-trsp(V_j)*V_k
     * (obviously antisymmetric).
     *
     * We need to do two things.
     *  - walk all w (sorted!) rows simultaneously, maintaining a
     *    priority queue for the "lowest next" column coordinate.
     *    By doing so, we'll be able to react on all column
     *    indices which happen to appear in several rows (because
     *    two consecutive pop()s from the queue will give the
     *    same index). So at select times, this gives a subset of
     *    w' rows, with that index appearing.  This step will
     *    cost O(w*L*log(w)).
     *  - Next, whenever we've identified a column index which
     *    appears multiple times (say w' times), build the w'*w'
     *    matrix trsp(V'_k)*V'_j-trsp(V'_j)*V'_k, where all
     *    vectors have length w'. Inspecting the 0 which appear
     *    in this matrix gives the places where we have
     *    exceptional cancellations. This can count as -1 in the
     *    total weight matrix.
     * Overall the cost becomes O(w*L*log(w)+w'*L'^2+w^2),
     * which should be somewhat faster.
     */

    size_t wj = batch.size();
    matrix<int> weights(wj,wj,0);

#ifdef FOR_DL
    typedef merge_matrix::row_value_t::exponent_type exponent_type;

    /* combination of rows i and j for i<j is done via:
     * row[i] * xx[i,j] - row[j] * xx[j,i]
     *
     * So in practive for i<j we expect xx[i,j] = ej[j] and xx[j,i] =
     * ej[i], with their gcd taken out.
     */
    matrix<exponent_type> xx(wj, wj, 0);
    std::vector<exponent_type> emax(wj);

    {
        /* We need to find ej in order to compute the matrix xx, but
         * we're happy to ditch it right thereafter */
        std::vector<exponent_type> ej(wj);
        for(size_t s = 0 ; s < wj ; s++) {
            for(auto const & v : batch[s]) {
                exponent_type e = v.exponent();
                if (std::abs(e) > emax[s])
                    emax[s] = std::abs(e);
                if (v.index() == j)
                    ej[s] = e;
            }
        }

        for(size_t s = 0 ; s < wj ; s++) {
            for(size_t t = 0 ; t < wj ; t++) {
                exponent_type d = gcd_int64(ej[s], ej[t]);
                xx(s,t) = ej[t] / d;
                xx(t,s) = ej[s] / d;
            }
        }
    }
#endif

    std::vector<size_t> ii(wj);

    typedef std::pair<size_t, size_t> pq_t;
    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> next;

    auto advance_and_push = [&ii, batch, &next](size_t r) {
        merge_matrix::shortlived_row_t const & R(batch[r]);
        if (ii[r] < R.size())
            next.push(std::make_pair(R[ii[r]].index(), r));
    };

    for(size_t r = 0 ; r < wj ; r++)
        advance_and_push(r);

    std::vector<size_t> local(wj);
#ifdef FOR_DL
    std::vector<exponent_type> ee(wj);
#endif

    for( ; !next.empty() ; ) {
        local.clear();
        pq_t me = next.top();
        /* By construction this loop will run at least once */
        for( ; !next.empty() && next.top().first == me.first ; ) {
            me = next.top();
            size_t r = me.second;
#ifdef FOR_DL
            ee[local.size()] = batch[r][ii[r]].exponent();
#endif
            local.push_back(r);
            next.pop();
            ii[r]++;
            advance_and_push(r);
        }

        if (local.size() == 1) {
            /* should occur pretty often */
            for(size_t s = 0 ; s < me.second ; s++)
                    weights(s, me.second)++;
            for(size_t s = me.second + 1 ; s < wj ; s++)
                    weights(me.second, s)++;
            continue;
        }

        if (me.first == j) {
            /* This one won't contribute to the weight */
            continue;
        }

        size_t zj = local.size();

        /* A priori we expect all of the values in the local[] array to
         * cause non-zero coefficients everywhere. */
        for(size_t s = 0 ; s < zj ; s++) {
            for(size_t t = 0 ; t < local[s] ; t++)
                weights(t, local[s])++;
            for(size_t t = local[s] + 1 ; t < wj ; t++)
                weights(local[s], t)++;
        }

        for(size_t s = 0 ; s < zj ; s++) {
            for(size_t t = s + 1 ; t < zj ; t++) {
                size_t i = local[s];
                size_t j = local[t];
#if FOR_DL
                /* what is the coefficient for this column when we
                 * combine rows local[s] and local[t] ? */
                exponent_type ei = ee[s];
                exponent_type ej = ee[t];
                exponent_type e = ei * xx(i,j) - ej * xx(j,i);
                if (e) {
                    /* see below -- at any rate counting two coefficients
                     * was too much ! */
                    weights(i, j)--;
                    weights(j, i)--;
                    continue;
                }
#endif
                /* ok, after all these do not have a coefficient
                 * here ! Two small catches: because our a priori update
                 * above added a contribution twice for this cell, we
                 * need to subtract two. Also, since we don't know
                 * whether i<j or j<i, we update both... */
                weights(i, j) -= 2;
                weights(j, i) -= 2;
            }
        }
    }

    /* make some final adjustments */
    for(size_t s = 0 ; s < wj ; s++) {
        for(size_t t = s + 1 ; t < wj ; t++) {
#if FOR_DL
            /* how about the max coefficent when we combine rows s and t ? */
            exponent_type e = emax[s] * std::abs(xx(s,t)) + emax[t] * std::abs(xx(t,s));
            if (e > 32)
                weights(s, t) = INT_MAX;
#endif
            /* make the matrix symmetric */
            weights(t,s) = weights(s,t);
        }
        weights(s,s) = INT_MAX;
    }

    return weights;
}/*}}}*/

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
void merge_matrix::parallel_pivot(collective_merge_operation & CM)
{
    std::vector<size_t> & my_col_indices(CM.my_col_indices);
    size_t n = CM.n;

    typedef collective_merge_operation::row_batch_t row_batch_t;
    typedef collective_merge_operation::row_batch_set_t row_batch_set_t;

    /* We'll get a vector of n vectors of up to w rows */
    row_batch_set_t row_batches = CM.merge_scatter_rows(*this);

    /* This has the same type, but with two fine subtleties:
     *  - all_merged_new_rows has n * comm_size fragments (all stored
     *    contiguously), because we receive fragments from everyone !
     *    (my_merged_new_rows has n entries).
     *  - each has wj-1 new rows, not wj.
     */
    row_batch_set_t my_merged_new_rows;

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
        row_batch_t const & batch(row_batches[q]);

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
            mst = minimum_spanning_tree(batch_compute_weights(batch, j));
        }

        /* Do the merge, but locally first */
        row_batch_t new_row_batch;
        for(auto const & edge : mst.second) {
            row_combiner C(batch[edge.first], batch[edge.second], j);
            new_row_batch.push_back(C.subtract_naively());
        }
        my_merged_new_rows.push_back(std::move(new_row_batch));
    }

    collective_merge_operation::row_batch_t all_merged_new_rows = CM.allgather_rows(my_merged_new_rows);

    std::vector<size_t> removed_row_indices;
    removed_row_indices.reserve(n * comm_size * cwmax);
    for(auto const & x : CM.all_row_indices)
        removed_row_indices.insert(removed_row_indices.end(), x.begin(), x.end());
    sort(removed_row_indices.begin(), removed_row_indices.end());

    std::vector<std::pair<size_t, row_weight_t>> inserted_row_indices;
    inserted_row_indices.reserve(all_merged_new_rows.size());
    /* This adds the rows, but skips the Rj update entirely */
    for(auto const & row : all_merged_new_rows) {
        size_t i = rows.size();
        rows.push_back(row);
        size_t newrow_size = row.size();
        MPI_Allreduce(MPI_IN_PLACE, &newrow_size, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
        inserted_row_indices.push_back(std::make_pair(i, newrow_size));
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
        }
    }
    for(auto i : removed_row_indices)
        remove_row_detached(i);
    MPI_Allreduce(&weight, &global_weight, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
    MPI_Allreduce(&ncols, &global_ncols, 1, MPI_MY_SIZE_T, MPI_SUM, comm);
}
/* }}} */

void merge_matrix::batch_Rj_update(std::vector<size_t> const & out, std::vector<std::pair<size_t, row_weight_t>> const & in)/*{{{*/
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

    auto advance_and_push = [this, &ii, in, out, &next](size_t r) {
        size_t i = r < out.size() ? out[r] : in[r-out.size()].first;
        auto R(rows[i]);
        for( ; ii[r] < R.second ; ii[r]++) {
            size_t j = R.first[ii[r]].index();
            ASSERT(is_my_col(j));
            size_t lj = j / comm_size;
            if (has_indirection(lj))
                break;
            col_weights[lj] += r < out.size() ? -1 : 1;
            /* don't do ncols-- just now, as we might increase the column
             * weight later on if j is both in and out */
        }
        if (ii[r] < R.second)
            next.push(std::make_pair(R.first[ii[r]].index(), r));
    };

    for(size_t r = 0 ; r < out.size() + in.size() ; r++)
        advance_and_push(r);

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
            ii[r]++;
            advance_and_push(r);
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
            if (!col_weights[lj]) ncols--;
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
        if ((ssize_t) w0 + sizediff > cwmax) {
            deactivate_indirection(lj, w1);
            /* don't forget to take it off the priority queue, or we'll
             * have a very bizarre situation */
            markowitz_table.remove(j);
            continue;
        }

        std::vector<size_t> Lout;
        std::vector<std::pair<size_t, row_weight_t>> Lin;
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

        const indirection_t * Rx0 = get_indirection(lj);
        temporary_indirection Rx1 = get_temporary_indirection(w1);

        /* Do the markowitz count while we're at it */
        size_t minrw = SIZE_MAX;

        // we rely on the fact that the Rj lists are sorted,
        // and so should be the lists of rows in and out.
        std::vector<std::pair<size_t, row_weight_t>>::const_iterator xin = Lin.begin();
        std::vector<size_t>::const_iterator xout = Lout.begin();

        const indirection_t * p = Rx0;
        indirection_t * q = Rx1;

        for( ; p != Rx0 + w0 && xout != Lout.end();) {
            for( ; p != Rx0 + w0 && p->first < *xout ; *q++ = *p++) {
                if (p->second < minrw)
                    minrw = p->second;
            }
            if (p == Rx0 + w0) break;
            /* here *p >= *xout */
            if (p->first == *xout)
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
        for( ; p != Rx0 + w0 ; *q++ = *p++) {
            if (p->second < minrw)
                minrw = p->second;
        }
        for( ; q != Rx1 + w1 ; *q++ = *xin++) {
            if (xin->second < minrw)
                minrw = xin->second;
        }
        ASSERT_ALWAYS(xin == Lin.end());

        replace_indirection(lj, Rx1, w1);
        if (w1 == 0) ncols--;

        /* last but not least ! */
        markowitz_table.update(std::make_pair(j, markowitz_count(j, minrw, w1)));
    }
}/*}}}*/

void merge_matrix::filter_merge_proposal(proposal_t & proposal, proposal_t & trash)/*{{{*/
{
    /* proposal is in decreasing order of priority. All are counted as
     * coefficient *decrease*, so that it's really a priority */
    // long maxprio = proposal.empty() ? LONG_MAX : proposal.front().second;
    long minprio = proposal.empty() ? LONG_MIN : proposal.back().second;
    long globalminprio = LONG_MAX;
    MPI_Allreduce(&minprio, &globalminprio, 1, MPI_LONG, MPI_MAX, comm);
    /*
    printf("rank %d has maxprio=%ld minprio=%ld global minprio=%ld\n",
            comm_rank, maxprio, minprio, globalminprio);
    */
    proposal_t ok;
    for(auto x : proposal) {
        if (x.second >= globalminprio)
            ok.push_back(x);
        else
            trash.push_back(x);
    }
    std::swap(ok, proposal);
}/*}}}*/

void merge_matrix::parallel_merge(size_t batch_size)/*{{{*/
{
    report_init();
    int spin = 0;
    report(true);
    for(cwmax = 2 ; cwmax <= maxlevel && spin < 10 ; cwmax++) {

        /* Note: we need the R table for remove_singletons. Note that
         * this function also updates markowitz_table */
        prepare_R_table();
        // expensive_check();

        remove_singletons_iterate();

        int wannabreak = spin && markowitz_table.empty();
        MPI_Allreduce(MPI_IN_PLACE, &wannabreak, 1, MPI_INT, MPI_LAND, comm);
        if (wannabreak) break;

        /* the control structure is global, as report() works on global
         * data */

        for( ; report() < target_density ; ) {
            int empty = markowitz_table.empty();
            MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_LAND, comm);
            if (empty) break;

            proposal_t merge_proposal;

            for( ; merge_proposal.size() < batch_size && !markowitz_table.empty() ; ) {
                markowitz_table_t::value_type q = markowitz_table.top();
                /* We're going to remove this column, so it's important that
                 * we remove it from the queue right now, otherwise our
                 * actions may cause a shuffling of the queue, and we could
                 * very well end up popping the wrong one if pop() comes too
                 * late.
                 */
                markowitz_table.pop();

                /* now that we have batch_Rj_update, this can be asserted */
                ASSERT_ALWAYS(col_weights[q.first / comm_size] <= (col_weight_t) cwmax);
                merge_proposal.push_back(q);
            }

            /* Make sure that the "best" merges really go first, and that
             * we don't hurry towards ordering sub-optimal merges. Some
             * fudge factor might be a good idea here.
             */

            proposal_t requeue;
            filter_merge_proposal(merge_proposal, requeue);
            /*
            printf("rank %d proposes %zu %d-merges, %zu to requeue\n", comm_rank, my_cols.size(), cwmax, requeue.size());
            */

            collective_merge_operation CM(comm, batch_size, cwmax);
            CM.announce_columns(merge_proposal);
            CM.deduce_and_share_rows(*this);
            proposal_t requeue2 = CM.discard_repeated_rows();
            discarded_merges[cwmax] += CM.my_discarded_merges;
            done_merges[cwmax] += CM.my_done_merges;

            /* requeue the merges which have been discarded. We need to
             * do it now, or otherwise we won't update their markowitz
             * count at all in batch_Rj_update */
            requeue.insert(requeue.end(), requeue2.begin(), requeue2.end());
            for(auto const & jp : requeue) {
                size_t j = jp.first;
                size_t lj = j / comm_size;
                if (!R_table[lj]) continue;
                markowitz_table.push(jp);
            }

            parallel_pivot(CM);
            // expensive_check();
            
        }
        // we remove the same set of rows on all nodes here.
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

    print_final_merge_stats();

    print_active_weight_count();

    if (!comm_rank)
        printf("# Output matrix has"
                " nrows=%zu ncols=%zu (excess=%zd) active_weight=%zu\n",
                nrows, global_ncols, nrows-global_ncols, global_weight);

}
/*}}}*/

int main(int argc, char *argv[])
{
    setbuf(stdout, NULL);
    MPI_Init(&argc, &argv);

    all_mpi_runtime_checks();

    char *argv0 = argv[0];

    int nthreads = 1;

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

    M.mpi_init(MPI_COMM_WORLD);

    int batch_size = (M.comm_size <= 8) ? (1 << (M.comm_size - 1)) : 128;
    param_list_parse_int(pl, "merge-batch-size", &batch_size);

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

    tt = seconds();

    M.read_rows(purgedname);

    if (!M.comm_rank)
        printf("# Time for filter_matrix_read: %2.2lfs\n", seconds() - tt);

    M.parallel_merge(batch_size);


#if 0
    /* resume from given history file if needed */
    if (resumename != NULL)
	resume(rep, mat, resumename);
#endif

    param_list_clear(pl);

    if (!M.comm_rank) {
        printf("Total merge time: %.2f seconds\n", seconds());
        print_timing_and_memory(stdout, wct0);
    }

    MPI_Finalize();
    return 0;
}
