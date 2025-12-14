#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mutex>
#include <vector>
#include <memory>
#include <algorithm>
#include <thread>
#include <list>
#include <functional>

#include <gmp.h>
#include "fmt/base.h"

#include "bit_vector.h"
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "filter_io.h"
#include "gmp_aux.h"
#include "gzip.h"
#include "macros.h"
#include "memalloc.h"
#include "mpz_poly.h"
#include "params.h"
#include "purgedfile.h"
#include "renumber.hpp"
#include "sm_utils.hpp"
#include "stats.h"
#include "timing.h"
#include "typedefs.h"
#include "verbose.h"
#include "utils_cxx.hpp"

#define DEBUG 0

static stats_data_t stats; /* struct for printing progress */

/*********************** mutex for multi threaded version ********************/
/* used as mutual exclusion lock for reading the status of logarithms */
static std::mutex lock;


// table of logs
struct logtab {// {{{
    uint64_t nprimes = 0;
    int nbsm = 0;
    cxx_mpz ell;
    cxx_cado_poly const & cpoly;
    std::vector<sm_side_info> const & sm_info;
    friend struct mpz_ro_accessor;
    friend struct mpz_rw_accessor;

  private:
    uint64_t nknown = 0;
    std::unique_ptr<mp_limb_t[]> data;

  public:
    uint64_t get_nknown() const { return nknown; }
    struct mpz_ro_accessor {
        /* This is a ***NON-OWNING*** mpz-like accessor for an integer
         * stored somewhere.
         *
         * A typical construct might be to use this as in:
         *
         * mpz_addmul(foo, log[i], bar);
         *
         * here log[i] would create an mpz_ro_accessor object, and then
         * call its mpz_srcptr converter (member function below) to
         * return an mpz_srcptr that is appropriate to pass to
         * mpz_addmul.
         *
         * The fine point is whether the dtor for the temporary
         * mpz_ro_accessor object is called before or after entering the
         * mpz_addmul function.
         *
         * The C++ standard says: after. [C++11 § 12.2.3] Namely:
         *
         *      When an implementation introduces a temporary object of a
         *      class that has a non-trivial constructor (12.1, 12.8), it
         *      shall ensure that a constructor is called for the
         *      temporary object. Similarly, the destructor shall be
         *      called for a temporary with a non-trivial destructor
         *      (12.4). Temporary objects are destroyed as the last step
         *      in evaluating the full-expression (1.9) that (lexically)
         *      contains the point where they were created.
         *
         * Bottom line: the construct mpz_addmul(foo, log[i], bar) is
         * safe. (However mpz_srcptr z = log[i]; followed by
         * mpz_addmul(foo, z, bar) is not !)
         */
      private:
        static int mpz_normalized_size(mp_limb_t * p, mp_size_t n)
        {
            for (; n && p[n - 1] == 0; n--)
                ;
            return n;
        }

      protected:
        mpz_t ugly;

      public:
        /* The following construct would be valid with gmp-6+
         * MPZ_ROINIT_N function, but unfortunately MPZ_ROINIT_N sets the
         * _mp_alloc field to zero, which gets in the way of our usage in
         * the rw_accessor below.
        mpz_ro_accessor(logtab const & log, mp_limb_t * place) : ugly
        MPZ_ROINIT_N(place, mpz_normalized_size(place, mpz_size(log.ell))) {}
         */
        mpz_ro_accessor(logtab const & log, mp_limb_t * place)
        {
            ugly->_mp_d = place;
            ugly->_mp_alloc = mpz_size(log.ell);
            ugly->_mp_size = mpz_normalized_size(place, mpz_size(log.ell));
        }
        operator mpz_srcptr() const { return ugly; }
        cxx_mpz as_cxx_mpz() const { return { static_cast<mpz_srcptr>(*this) }; }
    };
    struct mpz_rw_accessor : public mpz_ro_accessor {
      private:
        /* This is ***NON-OWNING*** as well, but has operator= which
         * carries over to the parent structure */
        uint64_t h;
        logtab & log;

      public:
        bool is_known() const { return log.is_known(h); }
        bool is_known_unlocked() const { return log.is_known_unlocked(h); }
        mpz_rw_accessor(logtab & log, uint64_t h, mp_limb_t * place)
            : mpz_ro_accessor(log, place)
            , h(h)
            , log(log)
        {
        }
        mpz_rw_accessor & operator=(mpz_srcptr v)
        {
            const std::lock_guard<std::mutex> dummy(lock);
            if (is_known_unlocked()) {
                cxx_mpz v0 = as_cxx_mpz();
                fmt::print(stderr,
                            "ERROR, inconsistent log for h = {} ;"
                            " we previously had {} in the database,"
                            " now we want to store {}\n",
                            h, v0, cxx_mpz(v));
                ASSERT_ALWAYS(mpz_cmp(v0, v) == 0);
                return *this;
            }
            if (mpz_cmp_ui(v, 0) < 0) {
                fmt::print(stderr,
                        "Warning, log is negative for h = {}\n", h);
                cxx_mpz vv;
                mpz_mod(vv, v, log.ell);
                (*this) = vv;
            } else if (mpz_cmp(v, log.ell) >= 0) {
                fmt::print(stderr, "Warning, log >= ell for h = {}\n", h);
                cxx_mpz vv;
                mpz_mod(vv, v, log.ell);
                (*this) = vv;
            } else {
                ASSERT_ALWAYS(mpz_size(v) <= mpz_size(log.ell));
                mpn_zero(ugly->_mp_d, mpz_size(log.ell));
                if (mpz_cmp_ui(v, 0) == 0) {
                    fmt::print(stderr,
                            "Warning, log is zero for h = {}\n", h);
                } else {
                    mpn_copyi(ugly->_mp_d, v->_mp_d, mpz_size(v));
                }
                // log of SM columns are not taken into account
                if (h < log.nprimes)
                    log.nknown++;
            }
            return *this;
        }
    };

  private:
    uint64_t smlog_index(int side, int idx_sm) const
    {
        uint64_t h = nprimes;
        for (int i = 0; i < side; i++) {
            h += sm_info[i].nsm;
        }
        return h + idx_sm;
    }

  public:
    /* The mpz_ro_accessor and mpz_rw_accessor structs are ugly hacks.
     * They're used in:
     *  - some mpz_addmul operations
     *  - parsing the table
     *  - printing the table
     */
    bool is_zero(uint64_t h) const
    {
        mpz_ro_accessor const z = (*this)[h];
        return mpz_cmp_ui((mpz_srcptr)z, 0) == 0;
    }
    bool is_known_unlocked(uint64_t h) const
    {
        mpz_ro_accessor const z = (*this)[h];
        return mpz_cmp((mpz_srcptr)z, ell) < 0;
    }
    bool is_known(uint64_t h) const
    {
        const std::lock_guard<std::mutex> dummy(lock);
        return is_known_unlocked(h);
    }
    mpz_ro_accessor operator[](uint64_t h) const
    {
        return { *this, data.get() + h * mpz_size(ell) };
    }
    mpz_rw_accessor operator[](uint64_t h)
    {
        return { *this, h, data.get() + h * mpz_size(ell) };
    }
    mpz_ro_accessor smlog(int side, int idx_sm) const
    {
        return (*this)[smlog_index(side, idx_sm)];
    }
    mpz_rw_accessor smlog(int side, int idx_sm)
    {
        return (*this)[smlog_index(side, idx_sm)];
    }
    void force_set(uint64_t h, mpz_srcptr v)
    {
        auto * p = data.get() + h * mpz_size(ell);
        ASSERT_ALWAYS(mpz_size(v) <= mpz_size(ell));
        mpn_zero(p, mpz_size(ell));
        mpn_copyi(p, v->_mp_d, mpz_size(v));
    }

    logtab(cxx_cado_poly const & cpoly, std::vector<sm_side_info> const & sm_info,
           uint64_t nprimes, mpz_srcptr ell)
        : nprimes(nprimes)
        , ell(ell)
        , cpoly(cpoly)
        , sm_info(sm_info)
    {
        for (int side = 0; side < cpoly->nb_polys; side++) {
            nbsm += sm_info[side].nsm;
        }
        data = std::make_unique<mp_limb_t[]>((nprimes + nbsm) * mpz_size(ell));
        /* set everything to the max value */
        std::fill_n(data.get(), (nprimes + nbsm) * mpz_size(ell), -1);
    }
    logtab(logtab const &) = delete;
};// }}}

/**** Relations structure used for computing the logarithms from the rels ****/
struct log_rel {// {{{
    weight_t nb_unknown = 0;
    ideal_merge_t * unknown = nullptr;
    cxx_mpz log_known_part;
    /* Return the number of unknown logarithms in the relation.
     * rels[i].nb_unknown may not be up-to-date (can only be greater than the actual
     * value) */
    weight_t nb_unknown_log(logtab & log)
    {
        weight_t j, k, len = nb_unknown;
        ideal_merge_t * p = unknown;
        int c = 0;
        for (j = 0, k = 0; k < len; k++) {
            bool known = log.is_known(p[k].id);

            if (!known) {
                if (j != k)
                    p[j] = p[k];
                j++;
            } else { // We know this log, add it to log_know_part
                mpz_addmul_si(log_known_part, log[p[k].id], p[k].e);
                c++;
            }
        }
        if (c)
            mpz_mod(log_known_part, log_known_part, log.ell);

        nb_unknown = j;
        return j;
    }

};// }}}


/************ Struct used for reading rels files with process_rels ***********/
struct read_data {// {{{
    std::vector<log_rel> rels;
    logtab & log;
    cxx_cado_poly & cpoly;
    std::vector<sm_side_info> const & sm_info;
    renumber_t & renum_tab;
    read_data(logtab & log, uint64_t nrels, cxx_cado_poly & cpoly,
              std::vector<sm_side_info> const & sm_info, renumber_t & renum_tab)
        : rels(nrels)
        , log(log)
        , cpoly(cpoly)
        , sm_info(sm_info)
        , renum_tab(renum_tab)
    {
        fflush(stdout);
    }
    read_data(read_data const &) = delete;

    private:
    /* In a relation with 1 missing logarithm of exponent e, compute its values,
     * i.e. compute   dest <- (-vlog / e) mod ell
     * Return 0 if dest was already known (i.e. computed between the call to
     * nb_unknown_log and the call to this functions), return 1 otherwise.
     */
    static unsigned int compute_missing_log(logtab::mpz_rw_accessor dest,
            mpz_t vlog, int32_t e, mpz_t ell)
    {
        unsigned int ret;
        cxx_mpz tmp = e;
        mpz_invert(tmp, tmp, ell);
        mpz_neg(vlog, vlog);
        mpz_mul(vlog, vlog, tmp);
        mpz_mod(tmp, vlog, ell);

        ret = !dest.is_known();
        if (ret)
            dest = tmp;

        return ret;
    }

    public:
    /* Compute all missing logarithms for relations in [start,end[.
     * Return the number of computed logarithms */
    uint64_t chunk(bit_vector not_used,
            uint64_t start, uint64_t end)
    {
        uint64_t computed = 0;

        for (size_t i = start; i < end; i++) {
            if (bit_vector_getbit(not_used, i)) {
                weight_t const nb = rels[i].nb_unknown_log(log);
                if (nb <= 1) {
                    bit_vector_clearbit(not_used, i);
                    mpz_ptr vlog = rels[i].log_known_part;
                    if (nb == 0 && mpz_cmp_ui(vlog, 0) != 0) {
                        fmt::print(stderr,
                                "Error, no unknown log in rel {} and sum"
                                " of log is not zero, sum is: {}\n",
                                i, cxx_mpz(vlog));
                        exit(EXIT_FAILURE);
                    } else if (nb == 1) {
                        ideal_merge_t const ideal = rels[i].unknown[0];
                        computed += compute_missing_log(log[ideal.id], vlog,
                                ideal.e, log.ell);
                    }
                }
            }
        }

        return computed;
    }

};// }}}


/************************** Dependency graph *********************************/
struct node_dep {
    enum { UNKNOWN, KNOWN_FROM_LOGFILE, RECONSTRUCTED } state = UNKNOWN;
    uint64_t i = 0;
};

struct graph_dep;

/***** Light relations structure (for constructing the dependency graph) *****/
struct light_rel {
    weight_t len = 0;
    index_t * needed = nullptr;

    /* Return the number of unknown logarithms in the relation and put in
     * h the last unknown logarithm of the relations (useful when the
     * number of unknown logarithms is 1) */
    weight_t dep_nb_unknown_log(graph_dep const & G, index_t & h) const;
};


struct graph_dep : public std::vector<node_dep> {// {{{
    explicit graph_dep(size_t size)
        : std::vector<node_dep>(size)
    {}
/* Set G[h].state accordingly to log[h] values */
    void set_log_already_known(logtab const & log)
    {
        for (uint64_t h = 0; h < log.nprimes; h++) {
            if (log.is_known(h))
                (*this)[h].state = node_dep::KNOWN_FROM_LOGFILE;
        }
    }

    uint64_t needed_rels_from_index(index_t h,
            std::vector<light_rel> & rels,
            bit_vector needed_rels)
    {
        if ((*this)[h].state == node_dep::UNKNOWN) {
            fmt::print(stderr,
                    "Error: logarithms of {:x} cannot be reconstructed "
                    "from this set of relations. Abort!\n",
                    h);
            abort();
        } else if ((*this)[h].state == node_dep::KNOWN_FROM_LOGFILE) {
            /* We know the wanted logarithm from linear algebra, no new relation is
             * necessary */
#if DEBUG >= 1
            fmt::print(stderr, "DEBUG: h = {:x} is known from logfile\n", h);
#endif
            return 0;
        } else {
            uint64_t const relnum = (*this)[h].i;
            bit_vector_setbit(needed_rels, relnum);
            uint64_t nadded = 1;
            weight_t const nb_needed = rels[relnum].len;

#if DEBUG >= 1
            fmt::print(stderr, "DEBUG: h = {:x} can be reconstructed\n", h);
            fmt::print(stderr, "DEBUG:     relation {} added\n", relnum);
            fmt::print(stderr, "DEBUG:     depends of %u others logs\n",
                    nb_needed - 1);
#endif

            for (weight_t j = 0; j < nb_needed; j++) {
                index_t const hh = rels[relnum].needed[j];
                if (!bit_vector_getbit(needed_rels, (*this)[hh].i)) {
                    nadded += needed_rels_from_index(hh, rels, needed_rels);
                }
            }
            return nadded;
        }
    }
};// }}}

inline weight_t light_rel::dep_nb_unknown_log(graph_dep const & G, index_t & h) const// {{{
{
    weight_t nb = 0;
    for (weight_t k = 0; k < len; k++) {
        int unknown;
        {
            const std::lock_guard<std::mutex> dummy(lock);
            unknown = G[needed[k]].state == node_dep::UNKNOWN;
        }

        if (unknown) { // we do not know the log if this ideal
            nb++;
            h = needed[k];
        }
    }
    return nb;
}// }}}



/************************ Handling of the SMs *******************************/
/* number of SM that must be used. */

/* Callback function called by filter_rels in compute_log_from_rels */
static void * thread_sm(void * context_data, earlyparsed_relation_ptr rel)
{
    read_data & data = *(read_data *)context_data;
    log_rel & lrel = data.rels[rel->num];

    mpz_ptr l = lrel.log_known_part;
    int64_t const a = rel->a;
    uint64_t const b = rel->b;

    uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */
    if (data.cpoly->nb_polys > 2) {
        for (weight_t i = 0; i < rel->nb; i++) {
            index_t const h = rel->primes[i].h;
            int const side = data.renum_tab.p_r_from_index(h).side;
            nonvoidside |= ((uint64_t)1) << side;
        }
        /* nonvoidside must *not* be a power of two. If it is, then we
         * have a nasty problem similar to bug 21707: in a sense, we have
         * true gem of a relation that yields a trivial norm on one side,
         * but it's really too bad that we have no effective way to check
         * for it. */
        ASSERT_ALWAYS(nonvoidside & (nonvoidside - 1));
        /* one thing we might do at this point is recompute the norm from
         * a, b, and data.cpoly->pols[side], and see if we get \pm1.
         */
    } else {
        nonvoidside = 3;
    }

    if (rel->sm_size) {
        /* use the SM values which are already present in the input file,
         * because some goodwill computed them for us.
         */
        int c = 0;
        for (int side = 0; side < data.cpoly->nb_polys; side++) {
            sm_side_info const & S = data.sm_info[side];
            if (S.nsm > 0 && (nonvoidside & (((uint64_t)1) << side))) {
#define xxxDOUBLECHECK_SM
#ifdef DOUBLECHECK_SM
                /* I doubt that this is really compatible with our
                 * changes in the SM mode.
                 */
                mpz_poly u;
                mpz_poly_init(u, MAX(1, S.f->deg - 1));
                mpz_poly_setcoeff_int64(u, 0, a);
                mpz_poly_setcoeff_int64(u, 1, -b);
                compute_sm_piecewise(u, u, S);
                ASSERT_ALWAYS(u->deg < S.f->deg);
                ASSERT_ALWAYS(u->deg == S.f->deg - 1);
                for (int i = 0; i < S.nsm; i++) {
                    if (S.mode == SM_MODE_LEGACY_PRE2018)
                        ASSERT_ALWAYS(mpz_cmp(u->coeff[S.f->deg - 1 - i],
                                              rel->sm[c + i]) == 0);
                    else
                        ASSERT_ALWAYS(mpz_cmp(u->coeff[i], rel->sm[c + i]) ==
                                      0);
                }

#endif
                ASSERT_ALWAYS(c + S.nsm <= rel->sm_size);
                for (int i = 0; i < S.nsm; i++, c++) {
                    mpz_addmul(l, data.log.smlog(side, i), rel->sm[c]);
                }
                mpz_mod(l, l, data.log.ell);
#ifdef DOUBLECHECK_SM
                mpz_poly_clear(u);
#endif
            }
        }
    } else {
        mpz_srcptr ell = data.log.ell;
        for (int side = 0; side < data.cpoly->nb_polys; side++) {
            sm_side_info const & S = data.sm_info[side];
            if (S.nsm > 0 && (nonvoidside & (((uint64_t)1) << side))) {
                cxx_mpz_poly u;
                mpz_poly_set_ab(u, a, b);
                S.compute_piecewise(u, u);
                ASSERT_ALWAYS(u->deg < S.f->deg);
                if (S.mode == SM_MODE_LEGACY_PRE2018) {
                    for (int i = S.f->deg - 1 - u->deg; i < S.nsm; i++) {
                        mpz_addmul(l, data.log.smlog(side, i),
                                   mpz_poly_coeff_const(u, S.f->deg - 1 - i));
                    }
                } else {
                    for (int i = 0; i < S.nsm; i++) {
                        mpz_addmul(l, data.log.smlog(side, i),
                                   mpz_poly_coeff_const(u, i));
                    }
                }
                mpz_mod(l, l, ell);
            }
        }
    }

    return NULL;
}

/****************** Computation of missing logarithms ************************/
/* Callback function called by filter_rels in compute_log_from_rels */
static void * thread_insert(void * context_data, earlyparsed_relation_ptr rel)
{
    read_data & data = *(read_data *)context_data;
    log_rel & lrel = data.rels[rel->num];
    unsigned int next = 0;
    ideal_merge_t buf[REL_MAX_SIZE];
    int c = 0;
    for (unsigned int i = 0; i < rel->nb; i++) {
        index_t const h = rel->primes[i].h;
        exponent_t const e = rel->primes[i].e;

        if (data.log.is_known(h)) {
            mpz_addmul_si(lrel.log_known_part, data.log[h], e);
            c++;
        } else
            buf[next++] = (ideal_merge_t) {.id = h, .e = e};
    }
    if (c)
        mpz_mod(lrel.log_known_part, lrel.log_known_part, data.log.ell);

    lrel.unknown = ideal_merge_my_malloc(next);
    lrel.nb_unknown = next;
    memcpy(lrel.unknown, buf, next * sizeof(ideal_merge_t));

    return nullptr;
}



/* Trivial structure containing the data necessary for reading rels for dep.
 * graph */
struct dep_read_data {// {{{
    std::vector<light_rel> & rels;
    graph_dep & G;
    dep_read_data(std::vector<light_rel> & rels, graph_dep & G)
        : rels(rels)
        , G(G)
    {
    }

    /* Compute all dependencies for relations in [start,end[.
     * Return the number of dependencies found */
    uint64_t chunk(bit_vector not_used, uint64_t start, uint64_t end)
    {
        uint64_t computed = 0;

        for (size_t i = start; i < end; i++) {
            if (bit_vector_getbit(not_used, i)) {
                index_t h = 0; // Placate gcc
                weight_t const nb = rels[i].dep_nb_unknown_log(G, h);
                if (nb <= 1) {
                    bit_vector_clearbit(not_used, i);
                    if (nb == 1) {
                        G[h].state = node_dep::RECONSTRUCTED;
                        G[h].i = i;
                        computed++;
                    }
                }
            }
        }

        return computed;
    }
};// }}}

/* Callback function called by filter_rels in compute_needed_rels */
static void * dep_thread_insert(void * context_data,
                                earlyparsed_relation_ptr rel)
{
    dep_read_data const & data = *(dep_read_data *)context_data;
    light_rel & lrel = data.rels[rel->num];
    unsigned int next = 0;
    index_t buf[REL_MAX_SIZE];

    for (unsigned int i = 0; i < rel->nb; i++) {
        index_t const h = rel->primes[i].h;
        if (data.G[h].state == node_dep::UNKNOWN)
            buf[next++] = h;
    }

    lrel.needed = index_my_malloc(next);
    lrel.len = next;
    memcpy(lrel.needed, buf, next * sizeof(ideal_merge_t));

    return NULL;
}

/******************** Code for multi thread version **************************/
#define SIZE_BLOCK 1024

template<typename D>
static uint64_t do_one_iter_mt(D & data,
                               bit_vector not_used, size_t nt,
                               uint64_t nrels)
{
    if (nt <= 1)
        return data.chunk(not_used, 0, nrels);

    uint64_t computed = 0;

    std::list<std::thread> threads;
    size_t alive = 0;

    // Main loop
    for (uint64_t i = 0 ; i < nrels || !threads.empty() ; ) {
        // Start / restart as many threads as allowed
        if (alive < nt && i < nrels) {
            const uint64_t nb = MIN(SIZE_BLOCK, nrels - i);
            threads.emplace_back([&](auto start, auto end) {
                    auto c = data.chunk(not_used, start, end);
                    const std::lock_guard<std::mutex> dummy(lock);
                    computed += c; }, i, i + nb);
            i += nb;
            alive++;
            continue;
        }
        // Wait for the oldest thread to finish in order to print result.
        threads.front().join();
        threads.pop_front();
        alive--;
    }

    return computed;
}

/***************** Important functions called by main ************************/
/* Read the logarithms computed by the linear algebra */
static void read_log_format_LA(logtab & log, char const * logfile,
                               char const * idealsfile,
                               std::vector<sm_side_info> const & sm_info,
                               int nb_polys)
{
    uint64_t i, ncols, col;
    index_t h;
    mpz_t tmp_log;
    FILE *flog = NULL, *fid = NULL;

    fmt::print("# Reading logarithms in LA format from {}\n", logfile);
    fmt::print("# Reading links between matrix columns and ideals from {}\n",
           idealsfile);
    fflush(stdout);
    flog = fopen_maybe_compressed(logfile, "r");
    FATAL_ERROR_CHECK(flog == NULL, "Cannot open file for reading logarithms");
    fid = fopen_maybe_compressed(idealsfile, "r");
    FATAL_ERROR_CHECK(fid == NULL, "Cannot open ideals file");

    if (fscanf(fid, "# %" SCNu64 "\n", &ncols) != 1) {
        fmt::print(stderr, "Error while reading first line of {}\n", idealsfile);
        abort();
    }

    mpz_init(tmp_log);
    i = 0;
    stats_init(stats, stdout, &i, nbits(ncols) - 5, "Read", "logarithms", "",
               "logs");
    while (fscanf(fid, "%" SCNu64 " %" SCNid "\n", &col, &h) == 2) {
        FATAL_ERROR_CHECK(col >= ncols, "Too big value of column number");
        FATAL_ERROR_CHECK(h >= log.nprimes, "Too big value of index");

        int const ret = gmp_fscanf(flog, "%Zd\n", tmp_log);
        FATAL_ERROR_CHECK(ret != 1,
                          "Error in file containing logarithms values");

        ASSERT_ALWAYS(col == i);
        log[h] = tmp_log;
        i++;
        if (stats_test_progress(stats))
            stats_print_progress(stats, i, 0, 0, 0);
    }
    stats_print_progress(stats, i, 0, 0, 1);
    ASSERT_ALWAYS(feof(fid));
    ASSERT_ALWAYS(i == ncols);

    for (int side = 0; side < nb_polys; side++) {
        for (int ism = 0; ism < sm_info[side].nsm; ism++) {
            int const ret = gmp_fscanf(flog, "%Zd\n", tmp_log);
            FATAL_ERROR_CHECK(ret != 1,
                              "Error in file containing logarithms values");
            log.smlog(side, ism) = tmp_log;
        }
    }
    /* If we are not at the end of the file, it means that it remains some
     * values and we do not know to what "ideals" they correspond. Probably an
     * error somewhere, it is better to abort. */
    ASSERT_ALWAYS(feof(flog));

    for (int side = 0; side < nb_polys; side++) {
        if (sm_info[side].nsm)
            fmt::print("# Logarithms for {} SM columns on side {}"
                    " were also read\n",
                    sm_info[side].nsm, side);
    }
    mpz_clear(tmp_log);
    fclose_maybe_compressed(flog, logfile);
    fclose_maybe_compressed(fid, idealsfile);
}

/* Read the logarithms in output format of reconstructlog */
static void read_log_format_reconstruct(logtab & log,
                                        MAYBE_UNUSED renumber_t const & renumb,
                                        char const * filename)
{
    uint64_t nread = 0;
    index_t h;
    mpz_t tmp_log;
    FILE * f = NULL;
    int ret;

    fmt::print("# Reading logarithms in reconstruct format from {}\n", filename);
    fflush(stdout);
    f = fopen_maybe_compressed(filename, "r");
    FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

    mpz_init(tmp_log);
    stats_init(stats, stdout, &nread, nbits(renumb.size()) - 5, "Read",
               "logarithms", "", "logs");
    for (index_t i = 0; i < renumb.number_of_additional_columns(); i++) {
        ret = gmp_fscanf(f, "%" SCNid " added column %Zd\n", &h, tmp_log);
        ASSERT_ALWAYS(ret == 2);
        ASSERT_ALWAYS(renumb.is_additional_column(h));
        nread++;
        log[h] = tmp_log;
    }
    for (index_t i = 0; i < renumb.number_of_bad_ideals(); i++) {
        ret = gmp_fscanf(f, "%" SCNid " bad ideals %Zd\n", &h, tmp_log);
        ASSERT_ALWAYS(ret == 2);
        ASSERT_ALWAYS(renumb.is_bad(h));
        nread++;
        log[h] = tmp_log;
    }
    while (gmp_fscanf(f, "%" SCNid " %*" SCNpr " %*d %*s %Zd\n", &h, tmp_log) ==
           2) {
        nread++;
        log[h] = tmp_log;
        if (stats_test_progress(stats))
            stats_print_progress(stats, nread, 0, 0, 0);
    }
    stats_print_progress(stats, nread, 0, 0, 1);

    for (int nsm = 0; nsm < log.nbsm; nsm++) {
        unsigned int n, side;
        if (nsm == 0) /* h was already read by previous gmp_fscanf */
        {
            ret = gmp_fscanf(f, "SM %u %u %Zd\n", &side, &n, tmp_log);
            ASSERT_ALWAYS(ret == 3);
        } else {
            ret = gmp_fscanf(f, "%" SCNid " SM %u %u %Zd\n", &h, &side, &n,
                             tmp_log);
            ASSERT_ALWAYS(ret == 4);
        }
        //    ASSERT_ALWAYS (n == nsm); // obsolete with new coding
        ASSERT_ALWAYS(h == (index_t)nsm + log.nprimes);
        log[h] = tmp_log;
    }
    ASSERT_ALWAYS(feof(f));

    mpz_clear(tmp_log);
    fclose_maybe_compressed(f, filename);
}

/* Write values of the known logarithms. */
static void write_log(char const * filename, logtab & log,
                      renumber_t const & tab,
                      std::vector<sm_side_info> const & sm_info)
{
    uint64_t i;
    FILE * f = NULL;

    fmt::print("# Opening {} for writing logarithms\n", filename);
    fflush(stdout);
    f = fopen_maybe_compressed(filename, "w");
    FATAL_ERROR_CHECK(f == NULL, "Cannot open file for writing");

    /* Divide all known logs by 'base' so that the first known non-zero
     * logarithm is equal to 1.
     * TODO: make a command line argument to choose this 'base'.
     */
    int base_already_set = 0;
    cxx_mpz base, scaled;
    for (i = 0; i < log.nprimes + log.nbsm; i++) {
        if (!log.is_known(i))
            continue;
        if (log.is_zero(i))
            continue;

        if (!base_already_set) {
            base_already_set = 1;
            /* base = 1/log[i] mod ell */
            int const ret = mpz_invert(base, log[i], log.ell);
            ASSERT_ALWAYS(ret != 0);
            mpz_set_ui(scaled, 1);
            log.force_set(i, scaled);
        } else {
            mpz_mul(scaled, log[i], base);
            mpz_mod(scaled, scaled, log.ell);
            log.force_set(i, scaled);
        }
    }

    uint64_t nknown = 0;
    stats_init(stats, stdout, &nknown, nbits(tab.size()) - 5, "Wrote",
               "known logarithms", "ideals", "logs");
    i = 0;
    for (auto it = tab.begin() ; it != tab.end() ; ++it, ++i) {
        if (!log.is_known(i))
            continue;
        nknown++;
        if (tab.is_additional_column(i)) {
            fmt::print(f, "{:x} added column {}\n", i, log[i].as_cxx_mpz());
        } else if (tab.is_bad(i)) {
            fmt::print(f, "{:x} bad ideals {}\n", i, log[i].as_cxx_mpz());
        } else {
            renumber_t::p_r_side const x = *it;

            if (x.side != tab.get_rational_side())
                fmt::print(f, "{:x} {:x} {} {:x} {}\n", i, x.p, x.side, x.r, log[i].as_cxx_mpz());
            else
                fmt::print(f, "{:x} {:x} {} rat {}\n", i, x.p, x.side, log[i].as_cxx_mpz());
        }
        if (stats_test_progress(stats))
            stats_print_progress(stats, nknown, i + 1, 0, 0);
    }
    stats_print_progress(stats, nknown, tab.size(), 0, 1);
    for (int nsm = 0 ; nsm < log.nbsm; nsm++) {
        auto const i = tab.size();
        // compute side
        int side, nsm_tot = sm_info[0].nsm, jnsm = nsm;
        for (side = 0; nsm >= nsm_tot; side++) {
            nsm_tot += sm_info[side + 1].nsm;
            jnsm -= sm_info[side].nsm;
        }
        ASSERT_ALWAYS((jnsm >= 0) && (jnsm < sm_info[side].nsm));
        if (log.is_zero(i + nsm)) {
            fmt::print("# Note: on side {}, log of SM number {} is zero\n", side,
                   jnsm);
        } else {
            ASSERT_ALWAYS(log.is_known(i + nsm));
        }
        fmt::print(f, "{:x} SM {} {} {}\n", i + nsm, side, jnsm, log[i + nsm].as_cxx_mpz());
    }

    uint64_t const missing = tab.size() - nknown;
    fmt::print("# factor base contains {} elements\n"
            "# logarithms of {} elements are known ({:.1f}%)\n"
            "# logarithms of {} elements are missing ({:.1f}%)\n",
            tab.size(),
            nknown, 100.0 * double_ratio(nknown, tab.size()),
            missing, 100.0 * double_ratio(missing, tab.size()));
    fclose_maybe_compressed(f, filename);
    ASSERT_ALWAYS(log.get_nknown() == nknown);
}

/* Given a filename, compute all the possible logarithms of ideals appearing in
 * the file. Return the number of computed logarithms.
 * Modify its first argument bit_vector needed_rels */
static uint64_t compute_log_from_rels(bit_vector needed_rels,
                                      char const * relspfilename,
                                      uint64_t nrels_purged,
                                      char const * relsdfilename,
                                      uint64_t nrels_del, uint64_t nrels_needed,
                                      int nt, read_data & data)
{
    double wct_tt0, wct_tt;
    uint64_t total_computed = 0, iter = 0, computed;
    uint64_t const nrels = nrels_purged + nrels_del;
    ASSERT_ALWAYS(nrels_needed > 0);

    /* Reading all relations */
    fmt::print("# Reading relations from {} and {}\n", relspfilename,
           relsdfilename);
    if (nrels_needed != nrels)
        fmt::print("# Parsing only {} needed relations out of {}\n",
               nrels_needed, nrels);
#if DEBUG >= 1
    fmt::print("# DEBUG: Using {} thread(s) for thread_sm\n", nt);
#endif
    fflush(stdout);
    char const * fic[3] = {(char const *)relspfilename,
                           (char const *)relsdfilename, nullptr};

    /* When purged.gz and relsdel.gz both have SM info included, we may
     * have an advantage in having more threads for thread_insert. Note
     * though that we'll most probably be limited by gzip throughput */

    int const ni = 1;
    int ns = nt;

    if (!filename_matches_one_compression_format(relspfilename) &&
        !filename_matches_one_compression_format(relsdfilename)) {
        fmt::print(
            "# Files {} and {} are uncompressed, limiting consumer threads\n",
            relspfilename, relsdfilename);
        ns = 4;
    }
    struct filter_rels_description desc[3] = {
        {.f = thread_insert, .arg = (void *)&data, .n = ni},
        {.f = thread_sm, .arg = (void *)&data, .n = ns},
        {.f = NULL, .arg = 0, .n = 0}};
    filter_rels2(fic, desc,
                 EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX |
                     EARLYPARSE_NEED_SM, /* It's fine (albeit slow) if we
                                            recompute them */
                 needed_rels, NULL);

    /* computing missing log */
    fmt::print("# Starting to compute missing logarithms from rels\n");

    /* adjust the number of threads based on the number of needed relations */
    double const ntm = ceil((nrels_needed + 0.0) / SIZE_BLOCK);
    if (nt > ntm)
        nt = (int)ntm;

    if (nt > 1)
        fmt::print("# Using multithread version with {} threads\n", nt);
    else
        fmt::print("# Using monothread version\n");

    wct_tt0 = wct_seconds();
    do {
        fmt::print("# Iteration {}: starting... [{} known logs]\n",
               iter, data.log.get_nknown());
        fflush(stdout);
        wct_tt = wct_seconds();

        /* Compute all missing logarithms possible.
         * Run through all the relations once.
         * Return the number of computed logarithms */
        computed = do_one_iter_mt(data, needed_rels, nt, nrels);
        total_computed += computed;

        fmt::print("# Iteration {}: {} new logarithms computed\n",
               iter, computed);
        fmt::print("# Iteration {} took {:.1f}s (wall-clock time).\n", iter,
               wct_seconds() - wct_tt);

        iter++;
    } while (computed);

    fmt::print("# Computing {} new logarithms took {:.1f}s (wall-clock "
           "time)\n",
           total_computed, wct_seconds() - wct_tt0);

    size_t const c = bit_vector_popcount(needed_rels);
    if (c != 0)
        fmt::print(stderr, "### Warning, {} relations were not used\n", c);

    return total_computed;
}

/* Given a filename, compute all the relations needed to compute the logarithms
 * appearing in the file.
 * needed_rels should be initialized before calling this function. Its size
 * must be nrels_purged + nrels_del.
 * Output:
 *    bit_vector needed_rels, where bits of needed rels are set.
 *    Return the number of needed_rels.*/
static uint64_t compute_needed_rels(bit_vector needed_rels,
                                    char const * relspfilename,
                                    uint64_t nrels_purged,
                                    char const * relsdfilename,
                                    uint64_t nrels_del, logtab & log,
                                    char const * wanted_filename, int nt)
{
    double wct_tt0, wct_tt;
    // uint64_t total_computed = 0;
    uint64_t iter = 0, computed;
    uint64_t const nrels = nrels_purged + nrels_del;
    graph_dep dep_graph(log.nprimes);
    std::vector<light_rel> rels(nrels);

    dep_graph.set_log_already_known(log);

    dep_read_data data(rels, dep_graph);

    /* Init bit_vector to remember which relations were already used */
    bit_vector_set(needed_rels, 1);

    /* Reading all relations */
    fmt::print("# Reading relations from {} and {}\n", relspfilename,
           relsdfilename);
    fflush(stdout);
    char const * fic[3] = {(char const *)relspfilename,
                           (char const *)relsdfilename, NULL};
    filter_rels(fic, (filter_rels_callback_t)&dep_thread_insert, (void *)&data,
                EARLYPARSE_NEED_INDEX, NULL, NULL);

    /* computing dependencies */
    fmt::print("# Starting to compute dependencies from rels\n");

    /* adjust the number of threads based on the number of relations */
    double const ntm = ceil((nrels + 0.0) / SIZE_BLOCK);
    if (nt > ntm)
        nt = (int)ntm;

    if (nt > 1)
        fmt::print("# Using multithread version with {} threads\n", nt);
    else
        fmt::print("# Using monothread version\n");

    wct_tt0 = wct_seconds();
    do {
        fmt::print("# Iteration {}: starting...\n", iter);
        fflush(stdout);
        wct_tt = wct_seconds();

        /* Compute all missing logarithms possible.
         * Run through all the relations once.
         * Return the number of computed logarithms */
        computed = do_one_iter_mt(data, needed_rels, nt, nrels);
        // total_computed += computed;

        fmt::print("# Iteration {}: {} new dependencies computed\n",
               iter, computed);
        fmt::print("# Iteration {} took {:.1f}s (wall-clock time).\n", iter,
               wct_seconds() - wct_tt);

        iter++;
    } while (computed);

    fmt::print("# Computing dependencies took {:.1f}s (wall-clock time)\n",
           wct_seconds() - wct_tt0);

    fmt::print("# Reading wanted logarithms from {}\n", wanted_filename);
    fflush(stdout);
    FILE * f = fopen_maybe_compressed(wanted_filename, "r");
    FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

    bit_vector_set(needed_rels, 0);
    index_t h;
    uint64_t nrels_necessary = 0, nwanted_log = 0;
    wct_tt = wct_seconds();
    while (fscanf(f, "%" SCNid "\n", &h) == 1) {
        FATAL_ERROR_CHECK(h >= log.nprimes, "Too big value of index");
        fmt::print("# Computing rels necessary for wanted log {:x}\n", h);
        fflush(stdout);
        auto nadded = dep_graph.needed_rels_from_index(h, rels, needed_rels);
        nrels_necessary += nadded;
        fmt::print("-> {} needed relations were added ({} so far)\n",
               nadded, nrels_necessary);
        nwanted_log++;
    }

    fclose_maybe_compressed(f, wanted_filename);
    fmt::print("# Reading {} wanted logarithms took {:.1f}s\n", nwanted_log,
           wct_seconds() - wct_tt);
    fmt::print("# {} relations are needed to compute these logarithms\n",
           nrels_necessary);
    ASSERT_ALWAYS(nrels_necessary == bit_vector_popcount(needed_rels));

    return nrels_necessary;
}

/********************* usage functions and main ******************************/
static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "log", "input file containing known logarithms");
    param_list_decl_usage(pl, "logformat",
                          "format of input log file: 'LA' or "
                          "'reconstruct' (default is 'LA')");
    param_list_decl_usage(pl, "ell", "group order (see sm -ell parameter)");
    param_list_decl_usage(pl, "out", "output file for logarithms");
    param_list_decl_usage(pl, "renumber", "input file for renumbering table");
    param_list_decl_usage(pl, "poly", "input polynomial file");
    param_list_decl_usage(pl, "ideals",
                          "link between matrix cols and ideals "
                          "(see replay -ideals parameter)");
    param_list_decl_usage(pl, "purged",
                          "file with purged relations "
                          "(see purge -out parameter)");
    param_list_decl_usage(pl, "relsdel",
                          "file with relations deleted by purge "
                          "(see purge -outdel parameter)");
    param_list_decl_usage(pl, "nrels",
                          "number of relations (same as purge "
                          "-nrels parameter)");
    param_list_decl_usage(pl, "partial",
                          "do not reconstruct everything "
                          "that can be reconstructed");
    param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
    param_list_decl_usage(pl, "nsm", "number of SM's to add on side 0,1,...");
    param_list_decl_usage(pl, "mt", "number of threads (default 1)");
    param_list_decl_usage(pl, "wanted", "file containing list of wanted logs");
    param_list_decl_usage(pl, "force-posix-threads",
                          "force the use of posix threads, do not rely on "
                          "platform memory semantics");
    param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
    verbose_decl_usage(pl);
}

static void usage(param_list pl, char const * argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    char const * argv0 = argv[0];

    uint64_t nrels_tot = 0, nrels_purged, nrels_del, nrels_needed;
    uint64_t nprimes;
    int mt = 1;
    int partial = 0;

    cxx_mpz ell;
    cxx_cado_poly cpoly;

    cxx_param_list pl;
    declare_usage(pl);
    argv++, argc--;

    param_list_configure_switch(pl, "partial", &partial);
    param_list_configure_switch(pl, "force-posix-threads",
                                &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY; /* Binary open for all files */
#endif

    if (argc == 0)
        usage(pl, argv0);

    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        fmt::print(stderr, "Unknown option: {}\n", argv[0]);
        usage(pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    char const * logfilename = param_list_lookup_string(pl, "log");
    char const * logformat = param_list_lookup_string(pl, "logformat");
    char const * idealsfilename = param_list_lookup_string(pl, "ideals");
    char const * relsdfilename = param_list_lookup_string(pl, "relsdel");
    char const * relspfilename = param_list_lookup_string(pl, "purged");
    char const * outfilename = param_list_lookup_string(pl, "out");
    char const * renumberfilename = param_list_lookup_string(pl, "renumber");
    char const * polyfilename = param_list_lookup_string(pl, "poly");
    char const * wantedfilename = param_list_lookup_string(pl, "wanted");
    param_list_parse_int(pl, "mt", &mt);
    char const * path_antebuffer =
        param_list_lookup_string(pl, "path_antebuffer");

    /* Some checks on command line arguments */
    if (!param_list_parse_mpz(pl, "ell", ell) || mpz_cmp_ui(ell, 0) <= 0) {
        fmt::print(stderr, "Error, missing -ell command line argument "
                        "(or ell <= 0)\n");
        usage(pl, argv0);
    }
    if (!param_list_parse_uint64(pl, "nrels", &nrels_tot) || nrels_tot == 0) {
        fmt::print(stderr, "Error, missing -nrels command line argument "
                        "(or nrels = 0)\n");
        usage(pl, argv0);
    }
    if (logfilename == NULL) {
        fmt::print(stderr, "Error, missing -log command line argument\n");
        usage(pl, argv0);
    }
    if (relspfilename == NULL) {
        fmt::print(stderr, "Error, missing -purged command line argument\n");
        usage(pl, argv0);
    }
    if (relsdfilename == NULL) {
        fmt::print(stderr, "Error, missing -relsdel command line argument\n");
        usage(pl, argv0);
    }
    if (outfilename == NULL) {
        fmt::print(stderr, "Error, missing -out command line argument\n");
        usage(pl, argv0);
    }
    if (renumberfilename == NULL) {
        fmt::print(stderr, "Error, missing -renumber command line argument\n");
        usage(pl, argv0);
    }
    if (polyfilename == NULL) {
        fmt::print(stderr, "Error, missing -poly command line argument\n");
        usage(pl, argv0);
    }
    if (mt < 1) {
        fmt::print(stderr, "Error: parameter mt must be at least 1\n");
        usage(pl, argv0);
    }

    if (logformat != NULL) {
        if (strcmp(logformat, "LA") != 0 &&
            strcmp(logformat, "reconstruct") != 0) {
            fmt::print(stderr,
                    "Error, unknown -formatlog argument. Must be 'LA' or "
                    "'reconstruct'\n");
            usage(pl, argv0);
        }
    }
    if ((logformat == NULL || strcmp(logformat, "LA") == 0) &&
        idealsfilename == NULL) {
        fmt::print(stderr, "Error, missing -ideals command line argument\n");
        usage(pl, argv0);
    }

    if (wantedfilename != NULL && !partial) {
        fmt::print(stderr, "Warning, -wanted command line argument is ignored if "
                        "-partial is not set\n");
    }

    if (!cado_poly_read(cpoly, polyfilename)) {
        fmt::print(stderr, "Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }

    /* Read number of sm to be printed from command line */
    std::vector<int> nsm_arg(cpoly->nb_polys, -1);
    param_list_parse_int_args_per_side(pl, "nsm", nsm_arg.data(),
                                       cpoly->nb_polys,
                                       ARGS_PER_SIDE_DEFAULT_AS_IS);

    for (int side = 0; side < cpoly->nb_polys; side++) {
        if (nsm_arg[side] < 0)
            continue;
        if (nsm_arg[side] > cpoly->pols[side]->deg) {
            fmt::print(stderr, "Error: nsm{}={} can not exceed the degree={}\n",
                    side, nsm_arg[side], cpoly->pols[side]->deg);
            exit(EXIT_FAILURE);
        }
    }

    char const * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

    if (param_list_warn_unused(pl)) {
        fmt::print(stderr, "Error, unused parameters are given\n");
        usage(pl, argv0);
    }

    set_antebuffer_path(argv0, path_antebuffer);

    /* Init data for computation of the SMs. */
    std::vector<sm_side_info> sm_info;
    for (int side = 0; side < cpoly->nb_polys; side++) {
        sm_info.emplace_back(cpoly->pols[side], ell, 0);
        sm_info[side].set_mode(sm_mode_string);
        fmt::print(stdout, "\n# Polynomial on side {}:\n# F[{}] = ", side, side);
        mpz_poly_fprintf(stdout, cpoly->pols[side]);
        fmt::print("# SM info on side {}:\n", side);
        sm_info[side].print(stdout);
        if (nsm_arg[side] >= 0)
            sm_info[side].nsm = nsm_arg[side]; /* command line wins */
        fmt::print("# Will use {} SMs on side {}\n", sm_info[side].nsm, side);

        /* do some consistency checks */
        if (sm_info[side].unit_rank != sm_info[side].nsm) {
            fmt::print(stderr,
                    "# On side {}, unit rank is {}, computing {} SMs ; "
                    "weird.\n",
                    side, sm_info[side].unit_rank, sm_info[side].nsm);
            /* for the 0 case, we haven't computed anything: prevent the
             * user from asking SM data anyway */
            ASSERT_ALWAYS(sm_info[side].unit_rank != 0);
        }
    }
    fflush(stdout);

    /* Reading renumber file */
    /* XXX legacy format insists on getting the badidealinfo file */
    fmt::print("\n###### Reading renumber file ######\n");
    renumber_t renumber_table(cpoly);
    renumber_table.read_from_file(renumberfilename, 1);
    nprimes = renumber_table.size();

    /* Read number of rows and cols on first line of purged file */
    {
        uint64_t nideals_purged;
        purgedfile_read_firstline(relspfilename, &nrels_purged,
                                  &nideals_purged);
        nrels_del = nrels_tot - nrels_purged;
    }

    /* Malloc'ing log tab and reading values of log */
    fmt::print("\n###### Reading known logarithms ######\n");
    fflush(stdout);

    logtab log(cpoly, sm_info, nprimes, ell);

    if (logformat == NULL || strcmp(logformat, "LA") == 0)
        read_log_format_LA(log, logfilename, idealsfilename, sm_info,
                           cpoly->nb_polys);
    else
        read_log_format_reconstruct(log, renumber_table, logfilename);

    /* Init bit_vector of rels that must be process by compute_log_from_rels */
    bit_vector rels_to_process;
    bit_vector_init(rels_to_process, nrels_tot);

    if (partial) {
        if (wantedfilename == NULL) {
            bit_vector_set(rels_to_process, 0);
            nrels_needed = 0;
        } else /* We compute needed rels for logarithms in wantedfilename */
        {
            fmt::print("\n###### Computing needed rels ######\n");
            nrels_needed = compute_needed_rels(
                rels_to_process, relspfilename, nrels_purged, relsdfilename,
                nrels_del, log, wantedfilename, mt);
        }
    } else {
        bit_vector_set(rels_to_process, 1);
        nrels_needed = nrels_tot;
    }

    /* Computing logs using rels in purged file */
    fmt::print("\n###### Computing logarithms using rels ######\n");
    if (nrels_needed > 0) {
        read_data data(log, nrels_purged + nrels_del, cpoly, sm_info,
                       renumber_table);

        compute_log_from_rels(rels_to_process, relspfilename, nrels_purged,
                              relsdfilename, nrels_del, nrels_needed, mt, data);
        fmt::print("# {} logarithms are known.\n", log.get_nknown());
    } else
        fmt::print(
            "# All wanted logarithms are already known, skipping this step\n");
    fflush(stdout);

    /* Writing all the logs in outfile */
    fmt::print("\n###### Writing logarithms in a file ######\n");
    write_log(outfilename, log, renumber_table, sm_info);

    /* freeing and closing */
    bit_vector_clear(rels_to_process);
    return EXIT_SUCCESS;
}
