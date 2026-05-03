#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <atomic>
#include <mutex>
#include <vector>
#include <algorithm>
#include <functional>
#include <thread>
#include <list>
#include <utility>
#include <string>
#include <limits>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "cado_poly.hpp"
#include "cxx_mpz.hpp"
#include "filter_io.hpp"
#include "gmp_aux.h"
#include "gzip.h"
#include "macros.h"
#include "memalloc.hpp"
#include "mpz_poly.h"
#include "params.hpp"
#include "renumber.hpp"
#include "sm_utils.hpp"
#include "stats.h"
#include "timing.h"
#include "typedefs.h"
#include "verbose.hpp"
#include "utils_cxx.hpp"
#include "fstream_maybe_compressed.hpp"

#define DEBUG 0

static stats_data_t stats; /* struct for printing progress */

/*********************** mutex for multi threaded version ********************/

// table of logs
struct logtab { /* {{{ */

    parameter_mandatory<std::string,
        "log",
        "input file containing known logarithms">
            logfilename;
    parameter_with_default<std::string,
        "logformat",
        "format of input log file: 'LA' or 'reconstruct'",
        "LA">
            logformat;
    parameter<std::string,
        "ideals",
        "link between matrix cols and ideals (for logformat=='LA')">
            idealsfilename;

    size_t nprimes = 0;
    std::vector<int> sm_per_side;
    int nbsm = 0;

    cxx_mpz const & ell;
    cxx_cado_poly const & cpoly;
    friend struct mpz_ro_accessor;
    friend struct mpz_rw_accessor;

    private:
    uint64_t nknown = 0;
    std::vector<mp_limb_t> data;
    mutable std::mutex lock;

    public:

    static void configure(cxx_param_list & pl) {
        pl.declare_usage_section("Know Logarithms provided as input");
        decltype(logfilename)::configure(pl);
        decltype(logformat)::configure(pl);
        decltype(idealsfilename)::configure(pl);
    }

    explicit logtab(cxx_param_list & pl,
            cxx_cado_poly const & cpoly,
            cxx_mpz const & ell)
        : logfilename(pl)
        , logformat(pl)
        , idealsfilename(pl)
        , ell(ell)
        , cpoly(cpoly)
    {
        if (logformat() != "LA" && logformat() != "reconstruct")
            pl.fail("-logformat must be 'LA' or 'reconstruct'");
        if (logformat() == "LA" && idealsfilename.is_default())
            pl.fail("missing -ideals");
        else if (logformat() == "reconstruct" && idealsfilename.is_provided())
            fmt::print(stderr, "# argument -ideals is ignored\n");
    }


    public:
    uint64_t get_nknown() const { return nknown; }
    struct mpz_ro_accessor { /* {{{ */
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
            static mp_size_t mpz_normalized_size(const mp_limb_t * p, mp_size_t n)
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
            mpz_ro_accessor(logtab const & log, const mp_limb_t * place)
            {
                ugly->_mp_d = const_cast<mp_limb_t *>(place);
                ugly->_mp_alloc = mpz_size(log.ell);
                ugly->_mp_size = mpz_normalized_size(place, mpz_size(log.ell));
            }
            operator mpz_srcptr() const { return ugly; }
            // cxx_mpz as_cxx_mpz() const { return { static_cast<mpz_srcptr>(*this) }; }
    }; /* }}} */
    struct mpz_rw_accessor : public mpz_ro_accessor { /* {{{ */
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
            mpz_rw_accessor & operator=(cxx_mpz const & v)
            {
                const std::scoped_lock dummy(log.lock);
                if (is_known_unlocked()) {
                    cxx_mpz v0 { *this };
                    fmt::print(stderr,
                            "ERROR, inconsistent log for h = {} ;"
                            " we previously had {} in the database,"
                            " now we want to store {}\n",
                            h, v0, v);
                    ASSERT_ALWAYS(v == v0);
                    return *this;
                }
                if (v < 0) {
                    fmt::print(stderr,
                            "Warning, log is negative for h = {}\n", h);
                    cxx_mpz vv;
                    mpz_mod(vv, v, log.ell);
                    (*this) = vv;
                } else if (v >= log.ell) {
                    fmt::print(stderr, "Warning, log >= ell for h = {}\n", h);
                    cxx_mpz vv;
                    mpz_mod(vv, v, log.ell);
                    (*this) = vv;
                } else {
                    ASSERT_ALWAYS(mpz_size(v) <= mpz_size(log.ell));
                    mpn_zero(ugly->_mp_d, mpz_size(log.ell));
                    if (v == 0) {
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
            operator mpz_ptr() { return ugly; }
    }; /* }}} */

    private:
    uint64_t smlog_index(int side, int idx_sm) const
    {
        uint64_t h = nprimes;
        for (int i = 0; i < side; i++)
            h += sm_per_side[i];
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
        const cxx_mpz z { (*this)[h] };
        return z == 0;
    }
    bool is_known_unlocked(uint64_t h) const
    {
        const cxx_mpz z { (*this)[h] };
        return z < ell;
    }
    bool is_known(uint64_t h) const
    {
        const std::scoped_lock dummy(lock);
        return is_known_unlocked(h);
    }
    mpz_ro_accessor operator[](uint64_t h) const
    {
        return { *this, data.data() + h * mpz_size(ell) };
    }
    mpz_rw_accessor operator[](uint64_t h)
    {
        return { *this, h, data.data() + h * mpz_size(ell) };
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
        auto * p = data.data() + h * mpz_size(ell);
        ASSERT_ALWAYS(mpz_size(v) <= mpz_size(ell));
        mpn_zero(p, mpz_size(ell));
        mpn_copyi(p, v->_mp_d, mpz_size(v));
    }

    /* This is a bit of a "late ctor". We just prepare the table with
     * placeholder values.
     */
    void prepare(size_t nprimes, std::vector<sm_side_info> const & sm_info)
    {
        this->nprimes = nprimes;
        this->nbsm = 0;
        sm_per_side.clear();
        for (auto const & s : sm_info) {
            sm_per_side.push_back(s.nsm);
            this->nbsm += s.nsm;
        }
        /* set everything to the max value */
        data.assign((nprimes + nbsm) * mpz_size(ell),
                std::numeric_limits<mp_limb_t>::max());
    }

    private:
    /* Read the logarithms computed by the linear algebra */
    void read_from_LA(std::vector<sm_side_info> const & sm_info)
    {
        uint64_t col;
        index_t h;
        cxx_mpz tmp_log;

        fmt::print("# Reading logarithms in LA format from {}\n",
                logfilename());
        fmt::print("# Reading links between matrix columns and ideals from {}\n",
                idealsfilename());
        fflush(stdout);

        ifstream_maybe_compressed flog(logfilename());
        if (!flog.good())
            throw cado::error("cannot open {}", logfilename());

        ifstream_maybe_compressed fideals(idealsfilename());
        if (!fideals.good())
            throw cado::error("cannot open {}", idealsfilename());

        size_t ncols;
        if (!(fideals >> expect("# ") >> ncols))
            throw cado::error("Error while reading first line of {}",
                    idealsfilename());

        uint64_t i = 0;
        stats_init(stats, stdout, &i, nbits(ncols) - 5,
                "Read", "logarithms", "", "logs");
        for( ; (fideals >> std::dec >> col >> std::hex >> h) && (flog >> tmp_log); i++) {
            FATAL_ERROR_CHECK(col >= ncols, "Column number too large");
            FATAL_ERROR_CHECK(h >= nprimes, "Index too large");
            ASSERT_ALWAYS(col == i);
            (*this)[h] = tmp_log;
            if (stats_test_progress(stats))
                stats_print_progress(stats, i, 0, 0, 0);
        }
        stats_print_progress(stats, i, 0, 0, 1);
        ASSERT_ALWAYS(fideals.eof());
        ASSERT_ALWAYS(i == ncols);

        for (int side = 0; side < cpoly.nsides(); side++) {
            for (int ism = 0; ism < sm_info[side].nsm; ism++) {
                if (!(flog >> tmp_log))
                    throw cado::error("error in {}", logfilename());
                smlog(side, ism) = tmp_log;
            }
        }
        /* If we are not at the end of the file, it means that it remains
         * some values and we do not know to what "ideals" they
         * correspond. Probably an error somewhere, it is better to
         * abort. */
        {
            /* technically this would abort on comment lines at the end
             * of the file... Not sure it's much to worry about.
             */
            std::string x;
            ASSERT_ALWAYS(!flog.eof() && !(flog>>x) && flog.eof());
        }

        for (int side = 0; side < cpoly.nsides(); side++) {
            if (sm_info[side].nsm)
                fmt::print("# Logarithms for {} SM columns on side {}"
                        " were also read\n",
                        sm_info[side].nsm, side);
        }
    }

    /* Read the logarithms in output format of reconstructlog */
    void read_from_reconstruct(renumber_t const & renumb)
    {
        cxx_mpz tmp_log;

        fmt::print("# Reading logarithms in reconstruct format from {}\n",
                logfilename());
        fflush(stdout);
        ifstream_maybe_compressed flog(logfilename());
        if (!flog.good())
            throw cado::error("cannot open {}", logfilename());

        /* stats_init wants to follow an uint64_t, but really it should
         * be a size_t instead... */
        uint64_t i = 0;
        stats_init(stats, stdout, &i, nbits(renumb.size()) - 5, "Read",
                "logarithms", "", "logs");

        index_t h;
        for (index_t k = 0; k < renumb.number_of_additional_columns(); k++, i++) {
            if (!(flog >> h >> expect("added column") >> tmp_log))
                throw cado::error("error in {}", logfilename());
            ASSERT_ALWAYS(renumb.is_additional_column(h));
            (*this)[h] = tmp_log;
        }
        for (index_t k = 0; k < renumb.number_of_bad_ideals(); k++, i++) {
            if (!(flog >> h >> expect("bad ideals") >> tmp_log))
                throw cado::error("error in {}", logfilename());
            ASSERT_ALWAYS(renumb.is_bad(h));
            (*this)[h] = tmp_log;
        }
        p_r_values_t p, r;
        int side;
        for(; flog >> std::hex >> h >> p >> side >> r >> std::dec >> tmp_log ; i++) {
            (*this)[h] = tmp_log;
            if (i % 1000 == 0) {
                const auto x = renumb.p_r_from_index(i);
                const renumber_t::p_r_side y { .p=p, .r=r, .side=side };
                ASSERT_ALWAYS(x == y);
            }
            if (stats_test_progress(stats))
                stats_print_progress(stats, i, 0, 0, 0);
        }
        stats_print_progress(stats, i, 0, 0, 1);

        /* This code is unfortunately not seen by coverage, but I think
         * there was a parsing bug here at sime point, with data being
         * read differently for the first SM and the following. (or maybe
         * it was a relic of old code).
         */
        for (int nsm = 0; nsm < nbsm; nsm++) {
            unsigned int n, side;
            if (!(flog >> h >> expect("SM") >> side >> n >> tmp_log))
                throw cado::error("error in {}", logfilename());
            ASSERT_ALWAYS(h == (index_t)nsm + nprimes);
            (*this)[h] = tmp_log;
        }
        ASSERT_ALWAYS(flog.eof());
    }

    public:

    void read(renumber_t const & tab, std::vector<sm_side_info> const & sm_info)
    {
        prepare(tab.size(), sm_info);
        if (logformat() == "LA")
            read_from_LA(sm_info);
        else if (logformat() == "reconstruct")
            read_from_reconstruct(tab);
        else
            ASSERT_ALWAYS(0);
    }

    /* Write values of the known logarithms. */
    void write(std::string const & filename,
            renumber_t const & tab,
            std::vector<sm_side_info> const & sm_info)
    {
        fmt::print("# Opening {} for writing logarithms\n", filename);
        fflush(stdout);
        ofstream_maybe_compressed f(filename);
        FATAL_ERROR_CHECK(!f.good(), "Cannot open file for writing");

        /* Divide all known logs by 'base' so that the first known non-zero
         * logarithm is equal to 1.
         * TODO: make a command line argument to choose this 'base'.
         */
        int base_already_set = 0;
        cxx_mpz base, scaled;
        for (size_t i = 0; i < nprimes + nbsm; i++) {
            if (!is_known(i))
                continue;
            if (is_zero(i))
                continue;

            if (!base_already_set) {
                base_already_set = 1;
                /* base = 1/log[i] mod ell */
                int const ret = mpz_invert(base, (*this)[i], ell);
                ASSERT_ALWAYS(ret != 0);
                mpz_set_ui(scaled, 1);
                force_set(i, scaled);
            } else {
                mpz_mul(scaled, (*this)[i], base);
                mpz_mod(scaled, scaled, ell);
                force_set(i, scaled);
            }
        }

        uint64_t nwritten = 0;
        stats_init(stats, stdout, &nwritten, nbits(tab.size()) - 5, "Wrote",
                "known logarithms", "ideals", "logs");
        size_t i = 0;
        for (auto it = tab.begin() ; it != tab.end() ; ++it, ++i) {
            if (!is_known(i))
                continue;
            nwritten++;
            if (tab.is_additional_column(i)) {
                fmt::print(f, "{:x} added column {}\n", i, cxx_mpz((*this)[i]));
            } else if (tab.is_bad(i)) {
                fmt::print(f, "{:x} bad ideals {}\n", i, cxx_mpz((*this)[i]));
            } else {
                auto const [ p, r, side ] = *it;

                if (side != tab.get_rational_side())
                    fmt::print(f, "{:x} {:x} {} {:x} {}\n", i, p, side, r, cxx_mpz((*this)[i]));
                else
                    fmt::print(f, "{:x} {:x} {} rat {}\n", i, p, side, cxx_mpz((*this)[i]));
            }
            if (stats_test_progress(stats))
                stats_print_progress(stats, nwritten, i + 1, 0, 0);
        }
        stats_print_progress(stats, nwritten, tab.size(), 0, 1);
        for (int nsm = 0 ; nsm < nbsm; nsm++) {
            auto const i = tab.size() + nsm;
            // compute side
            int side, nsm_tot = sm_info[0].nsm, jnsm = nsm;
            for (side = 0; nsm >= nsm_tot; side++) {
                nsm_tot += sm_info[side + 1].nsm;
                jnsm -= sm_info[side].nsm;
            }
            ASSERT_ALWAYS((jnsm >= 0) && (jnsm < sm_info[side].nsm));
            if (is_zero(i)) {
                fmt::print("# Note: on side {}, log of SM number {} is zero\n", side,
                        jnsm);
            } else {
                ASSERT_ALWAYS(is_known(i));
            }
            fmt::print(f, "{:x} SM {} {} {}\n", i, side, jnsm, cxx_mpz((*this)[i]));
        }

        uint64_t const missing = tab.size() - nwritten;
        fmt::print("# factor base contains {} elements\n"
                "# logarithms of {} elements are known ({:.1f}%)\n"
                "# logarithms of {} elements are missing ({:.1f}%)\n",
                tab.size(),
                nwritten, 100.0 * double_ratio(nwritten, tab.size()),
                missing, 100.0 * double_ratio(missing, tab.size()));
        ASSERT_ALWAYS(get_nknown() == nwritten);
    }
}; /* }}} */

template <typename T>
concept container_with_h = requires {
    T().begin();
    T().end();
    T()[0].h;
};

/* this is inspired by purge_matrix.hpp */
template<typename T>
struct reconstruct_matrix : simple_minded_chunk_allocator<T> /* {{{ */ {
    std::vector<std::pair<size_t, T *>> rows;
    protected:
    mutable std::mutex m;
    public:

    /* Insert row i into the matrix, containing the primes in the given
     * list that are within the column index range [col0, col1).
     *
     * remaining_rows increases by one
     * remaining_columns increases if needed.
     *
     * the matrix is reallocated if needed, with lock protection.
     */
    template <container_with_h C,
             typename Predicate,
             typename Proj = std::identity>
    void new_row(size_t i, C & primes, Predicate const & f, Proj const & proj = Proj())
    {
        auto * wprimes = primes.begin();

        for (auto const & c: primes) {
            if (f(c))
                *wprimes++ = c;
        }

        size_t active_weight = wprimes - primes.begin();

        auto * p = simple_minded_chunk_allocator<T>::alloc(active_weight);

        {
            auto * pp = p;
            for(const auto * q = primes.begin() ; q != wprimes ; q++)
                *pp++ = proj(*q);
        }

        std::scoped_lock const dummy(m);
        if (i >= rows.size())
            rows.insert(rows.end(), i + 1 - rows.size(), {0, nullptr});
        rows[i] = {active_weight, p};
    }

#define SIZE_BLOCK 1024
    template<typename Derived>
    size_t do_one_iter_mt(std::vector<bool> & not_yet_used, size_t nt)
    {
        auto & me = dynamic_cast<Derived &>(*this);

        if (nt <= 1)
            return me.chunk(not_yet_used, 0, rows.size());

        std::atomic<size_t> computed = 0;

        std::list<std::thread> threads;
        size_t alive = 0;

        // Main loop
        for (size_t i = 0 ; i < rows.size() || !threads.empty() ; ) {
            // Start / restart as many threads as allowed
            if (alive < nt && i < rows.size()) {
                const size_t nb = MIN(SIZE_BLOCK, rows.size() - i);
                threads.emplace_back([&](auto start, auto end) {
                        computed += me.chunk(not_yet_used, start, end);
                        }, i, i + nb);
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

    /* we need this in order to make the type polymorphic, and allow
     * do_one_iter_mt above to access the derived type safely
     */
    virtual ~reconstruct_matrix() = default;
}; /* }}} */


struct input_relation_pool /* {{{ */ {
    parameter_mandatory<std::string,
        "purged",
        "file with purged relations (see purge -out parameter)">
            relspfilename;
    parameter_mandatory<std::string,
        "relsdel",
        "file with relations deleted by purge (see purge -outdel parameter)">
            relsdfilename;
    /* TODO: we want to make this optional, as we did in purge.
    parameter_with_default<size_t,
        "nrels",
        "number of relations (same as purge -nrels parameter",
        "0">
            nrels;
     */
    parameter_switch<"force-posix-threads",
        "force the use of posix threads,"
            " do not rely on platform memory semantics">
            force_posix_threads;
    parameter<std::string,
        "path_antebuffer",
        "path to antebuffer program">
            path_antebuffer;
    static void configure(cxx_param_list & pl) {
        pl.declare_usage_section("Input specification (files)");
        decltype(relspfilename)::configure(pl);
        decltype(relsdfilename)::configure(pl);
        // decltype(nrels)::configure(pl);
        decltype(force_posix_threads)::configure(pl);
        decltype(path_antebuffer)::configure(pl);
    }
    explicit input_relation_pool(cxx_param_list & pl)
        : relspfilename(pl)
        , relsdfilename(pl)
        // , nrels(pl)
        , force_posix_threads(pl)
        , path_antebuffer(pl)
    {
        if (path_antebuffer.is_provided())
            set_antebuffer_path(pl.binary_name().c_str(),
                    path_antebuffer.parameter_value().c_str());
    }

    template<typename relation_type, typename... Args>
    size_t filter(std::vector<bool> const * only_some_rels, Args && ...args) const
    {
        const std::vector<std::string> files { relspfilename, relsdfilename };
        if (force_posix_threads) {
            using locking_type = cado::filter_io_details::ifb_locking_posix;
            return filter_rels<locking_type, relation_type>(files,
                    only_some_rels, nullptr,
                    std::forward<Args>(args)...);
        } else {
            using locking_type = cado::filter_io_details::ifb_locking_lightweight;
            return filter_rels<locking_type, relation_type>(files,
                    only_some_rels, nullptr,
                    std::forward<Args>(args)...);
        }
    }
    bool both_files_uncompressed() const {
        for(auto const & s : { relspfilename(), relsdfilename() })
            if (filename_matches_one_compression_format(s.c_str()))
                return false;
        return true;
    }
};
namespace fmt {
    template<>
    struct formatter<input_relation_pool> : public formatter<string_view> {
        // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
        auto format(input_relation_pool const & R, format_context & ctx) const
        {
            return fmt::format_to(ctx.out(), "input relation pool consisting of {} and {}", R.relspfilename(), R.relsdfilename());
        }
    };
} /* namespace fmt */
/* }}} */

/* This structure is the main device of reconstructlog. We do multiple
 * passes over the input, and work on a table that collapses the known
 * logs as we encounter them.
 */
struct compute_log_from_rels    /* {{{ */
    : public reconstruct_matrix<prime_type_for_indexed_relations>
{

    using relation_type =
        cado::relation_building_blocks::sm_block<
        cado::relation_building_blocks::primes_block<
            prime_type_for_indexed_relations,
        cado::relation_building_blocks::ab_block<uint64_t, 16>>>;

    std::vector<cxx_mpz> known_log_part_of_rows;

    logtab & log;
    cxx_cado_poly const & cpoly;
    std::vector<sm_side_info> const & sm_info;
    renumber_t & renum_tab;

    /* compare the ideals stored in relation i with the log table. Adjust
     * the known part according to what we have in the table, and discard
     * the rest. Return the new number of unknown logs.
     */
    void adjust_number_of_unknowns(size_t i) { /* {{{ */
        auto [ n, p ] = rows[i];
        auto & L = known_log_part_of_rows[i];

        auto * q = p;
        for(size_t i = 0 ; i < n ; i++) {
            auto const & [ h, e ] = p[i];
            if (log.is_known(h))
                mpz_addmul_si(L, log[h], e);
            else
                *q++ = p[i];
        }
        if (q != p + n)
            mpz_mod(L, L, log.ell);
        rows[i] = { q-p, p };
    } /* }}} */

    /* Compute all missing logarithms for relations in [start,end[.
     * Return the number of computed logarithms */
    size_t chunk(std::vector<bool> & not_yet_used, uint64_t start, uint64_t end) { /* {{{ */
        size_t computed = 0;

        for (size_t i = start; i < end; i++) {
            if (!not_yet_used[i])
                continue;

            adjust_number_of_unknowns(i);
            auto [ nb, row ] = rows[i];

            if (nb >= 2)
                continue;

            not_yet_used[i] = false;

            cxx_mpz const & vlog = known_log_part_of_rows[i];

            if (nb == 0 && vlog != 0)
                throw cado::error(
                        "no unknown log in rel {} and sum"
                        " of log is not zero, sum is: {}",
                        i, cxx_mpz(vlog));

            if (nb == 0)
                continue;

            auto [ h, e ] = row[0];

            if (log.is_known(h))
                continue;

            /* compute log[h] <- (-vlog / e) mod ell */
            cxx_mpz tmp = e;
            mpz_invert(tmp, tmp, log.ell);
            mpz_neg(tmp, tmp);
            mpz_mul(tmp, vlog, tmp);
            mpz_mod(tmp, tmp, log.ell);
            /* log[h] is a tricky beast here. Only use it via its
             * operator= overload, and nothing else.
             */
            log[h] = tmp;
            computed++;

            /* TODO: should we release memory at this point?  Perhaps.
             * Well, the trick is that the simple-minded allocator
             * doesn't support freeing anyway...  */
            rows[i] = { 0, nullptr };
        }
        return computed;
    } /* }}} */

    compute_log_from_rels(logtab & log, cxx_cado_poly & cpoly,
              std::vector<sm_side_info> const & sm_info,
              renumber_t & renum_tab)
        : log(log)
        , cpoly(cpoly)
        , sm_info(sm_info)
        , renum_tab(renum_tab)
    {
        fflush(stdout);
    }
    compute_log_from_rels(compute_log_from_rels const &) = delete;

    /* Callback function called by filter_rels in compute_log_from_rels */
    void thread_insert(relation_type & rel) /* {{{ */
    {
        int c = 0;
        {
            std::scoped_lock const dummy(reconstruct_matrix::m);
            auto const i = rel.num;
            auto & v = known_log_part_of_rows;
            if (i >= v.size())
                v.insert(v.end(), i + 1 - v.size(), 0);
        }
        cxx_mpz & vlog = known_log_part_of_rows[rel.num];
        auto pred = [&](prime_type_for_indexed_relations const & p) {
            auto const [ h, e ] = p;
            const bool b = log.is_known(h);
            if (b) {
                mpz_addmul_si(vlog, log[h], e);
                c++;
            }
            return !b;
        };
        new_row(rel.num, rel.primes, pred);
        if (c)
            mpz_mod(vlog, vlog, log.ell);
    } /* }}} */

    /* Callback function called by filter_rels in compute_log_from_rels */
    void thread_sm(relation_type const & rel) /* {{{ */
    {
        auto const & a = rel.a;
        auto const & b = rel.b;

        int nonvoidside = 0; /* bit vector of which sides appear in the rel */

        if (cpoly.nsides() > 2) {
            for (auto const & [ h, e ] : rel.primes) {
                int const side = renum_tab.p_r_from_index(h).side;
                nonvoidside |= 1 << side;
            }
            /* nonvoidside must *not* be a power of two. If it is, then we
             * have a nasty problem similar to bug 21707: in a sense, we have
             * true gem of a relation that yields a trivial norm on one side,
             * but it's really too bad that we have no effective way to check
             * for it. */
            ASSERT_ALWAYS(nonvoidside & (nonvoidside - 1));
            /* one thing we might do at this point is recompute the norm from
             * a, b, and data.cpoly[side], and see if we get \pm1.
             */
        } else {
            nonvoidside = 3;
        }

        cxx_mpz sum_sm_contribs = 0;
        if (!rel.sm.empty()) {
            /* use the SM values which are already present in the input file,
             * because some goodwill computed them for us.
             */
            size_t c = 0;
            for (int side = 0; side < cpoly.nsides(); side++) {
                sm_side_info const & S = sm_info[side];
                if (S.nsm > 0 && (nonvoidside & (1 << side))) {
#define xxxDOUBLECHECK_SM
#ifdef DOUBLECHECK_SM
                    /* I doubt that this is really compatible with our
                     * changes in the SM mode.
                     */
                    cxx_mpz_poly u;
                    mpz_poly_setcoeff_int64(u, 0, a);
                    mpz_poly_setcoeff_int64(u, 1, -b);
                    compute_sm_piecewise(u, u, S);
                    ASSERT_ALWAYS(u->deg < S.f->deg);
                    ASSERT_ALWAYS(u->deg == S.f->deg - 1);
                    for (int i = 0; i < S.nsm; i++) {
                        if (S.mode == SM_MODE_LEGACY_PRE2018)
                            ASSERT_ALWAYS(mpz_cmp(u->coeff[S.f->deg - 1 - i],
                                        rel.sm[c + i]) == 0);
                        else
                            ASSERT_ALWAYS(mpz_cmp(u->coeff[i], rel.sm[c + i]) ==
                                    0);
                    }

#endif
                    ASSERT_ALWAYS(c + S.nsm <= rel.sm.size());
                    for (int i = 0; i < S.nsm; i++, c++)
                        mpz_addmul(sum_sm_contribs, log.smlog(side, i),
                                rel.sm[c]);
                }
            }
        } else {
            for (int side = 0; side < cpoly.nsides(); side++) {
                sm_side_info const & S = sm_info[side];
                if (S.nsm > 0 && (nonvoidside & (1 << side))) {
                    cxx_mpz_poly u;
                    mpz_poly_set_ab(u, a, b);
                    S.compute_piecewise(u, u);
                    ASSERT_ALWAYS(u->deg < S.f->deg);
                    if (S.mode == SM_MODE_LEGACY_PRE2018) {
                        for (int i = S.f->deg - 1 - u->deg; i < S.nsm; i++)
                            mpz_addmul(sum_sm_contribs, log.smlog(side, i),
                                    mpz_poly_coeff_const(u, S.f->deg - 1 - i));
                    } else {
                        for (int i = 0; i < S.nsm; i++)
                            mpz_addmul(sum_sm_contribs, log.smlog(side, i),
                                    mpz_poly_coeff_const(u, i));
                    }
                }
            }
        }
        /* thread_insert may still insert relations at this
         * point, which might trigger reallocation of the
         * known_log_part_of_rows vector */
        std::scoped_lock const dummy(reconstruct_matrix::m);
        cxx_mpz & vlog = known_log_part_of_rows[rel.num];
        mpz_add(vlog, vlog, sum_sm_contribs);
        mpz_mod(vlog, vlog, log.ell);
    } /* }}} */

    /* Compute all the possible logarithms of ideals given the input.
     * Return the number of computed logarithms, and
     * modify needed_rels */
    size_t process_input(std::vector<bool> * needed_rels, /* {{{ */
            input_relation_pool const & input,
            int nt)
    {
        fmt::print("\n###### Computing logarithms using rels ######\n");

        double wct_tt0, wct_tt;
        uint64_t total_computed = 0, iter = 0, computed;
        // uint64_t const nrels = nrels_purged + nrels_del;
        // ASSERT_ALWAYS(nrels_needed > 0);

        /* Reading all relations */
        fmt::print("# Reading relations from {}\n", input);
#if DEBUG >= 1
        fmt::print("# DEBUG: Using {} thread(s) for thread_sm\n", nt);
#endif
        fflush(stdout);

        /* When purged.gz and relsdel.gz both have SM info included, we may
         * have an advantage in having more threads for thread_insert. Note
         * though that we'll most probably be limited by gzip throughput */

        int const ni = 1;
        int ns = nt;

        if (input.both_files_uncompressed() && ns > 4) {
            fmt::print("# {} is uncompressed, limiting consumer threads\n",
                    input);
            ns = 4;
        }

        /* if needed_rels != nullptr, nrels is the popcount of
         * *needed_rels.
         */
        size_t nrels = input.filter<relation_type>(needed_rels,
                cado::filter_io_details::multithreaded_call(ni,
                    [&](relation_type & rel) { thread_insert(rel); }),
                cado::filter_io_details::multithreaded_call(ns,
                    [&](relation_type const & rel) { thread_sm(rel); })
                );

        /* computing missing log */
        fmt::print("# Starting to compute missing logarithms from rels\n");

        std::vector<bool> not_yet_used;
        if (needed_rels) {
            not_yet_used = std::move(*needed_rels);
        } else {
            not_yet_used.assign(nrels, true);
        }
        /* TODO: adjust the number of threads based on the number of
         * needed relations. If we don't have many in the bitmap, then
         * there's little point in firing many threads...
         */
        // nt = std::min(nt, static_cast<int>(double_ratio(nrels_needed, SIZE_BLOCK)));
        if (nt > 1)
            fmt::print("# Using multithread version with {} threads\n", nt);
        else
            fmt::print("# Using monothread version\n");

        wct_tt0 = wct_seconds();
        do {
            fmt::print("# Iteration {}: starting... [{} known logs]\n",
                    iter, log.get_nknown());
            fflush(stdout);
            wct_tt = wct_seconds();

            /* Compute all missing logarithms possible.
             * Run through all the relations once.
             * Return the number of computed logarithms */
            computed = do_one_iter_mt<compute_log_from_rels>(not_yet_used, nt);
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

        size_t const c = std::ranges::count(not_yet_used, true);
        if (c != 0)
            fmt::print(stderr, "### Warning, {} relations were not used\n", c);

        fmt::print("# {} logarithms are known.\n", log.get_nknown());

        return total_computed;
    } /* }}} */
}; /* }}} */


/************************** Dependency graph *********************************/

struct compute_needed_rels    /* {{{ */
: public reconstruct_matrix<index_t>
{
    enum class state { UNKNOWN, KNOWN_FROM_LOGFILE, RECONSTRUCTED };

    using node_dep = std::pair<state, size_t>;

    using relation_type =
        cado::relation_building_blocks::primes_block<
        prime_type_for_indexed_relations,
        cado::relation_building_blocks::ab_ignore<16>>;

    /* This dependency graph maps the set of column indices to a
     * node_dep, which says if the vlog for this column can be computed,
     * and where from (which row index).
     *
     * The reconstruct_matrix contains no column index for which the
     * state is KNOWN_FROM_LOGFILE;
     *
     * If the graph says that the log of column index h can be
     * reconstructed from row i, then the same holds for all column
     * indices found in rows[i] (which isn't necessarily of weight 1).
     */
    class graph_dep : private std::vector<node_dep> { /* {{{ */
        using super = std::vector<node_dep>;
        mutable std::mutex lock;
        public:
        using super::size;

        explicit graph_dep(size_t size)
            : super(size)
        {}

        /* Set G[h].state accordingly to log[h] values */
        void set_log_already_known(logtab const & log)
        {
            const std::scoped_lock dummy(lock);
            for (uint64_t h = 0; h < log.nprimes; h++) {
                if (log.is_known(h))
                    super::operator[](h) = { state::KNOWN_FROM_LOGFILE, 0 };
            }
        }

        bool is_unknown(size_t n) const {
            const std::scoped_lock dummy(lock);
            return super::operator[](n).first == state::UNKNOWN;
        }

        bool is_known(size_t n) const {
            const std::scoped_lock dummy(lock);
            return super::operator[](n).first != state::UNKNOWN;
        }

        bool is_known_from_logfile(size_t n) const {
            const std::scoped_lock dummy(lock);
            return super::operator[](n).first == state::KNOWN_FROM_LOGFILE;
        }
        bool is_known_reconstructed(size_t n) const {
            const std::scoped_lock dummy(lock);
            return super::operator[](n).first == state::RECONSTRUCTED;
        }

        void mark_node_reconstructed(index_t h, size_t i) {
            const std::scoped_lock dummy(lock);
            super::operator[](h) = { state::RECONSTRUCTED, i };
        }

        node_dep operator[](size_t n) const {
            const std::scoped_lock dummy(lock);
            return super::operator[](n);
        }
    }; /* }}} */

    graph_dep G;

    /* Given the dependency graph G, determine the set of relations that
     * are needed to compute the vlog of column h. Store the bitmask of
     * the needed relations to needed_rels (in addition to the ones that
     * are already marked), and return the number of newly added
     * relations in thie bitmap.
     */
    size_t needed_rels_from_index(index_t h, std::vector<bool> & needed_rels) /* {{{ */
    {
        auto [ s, relnum ] = G[h];

        if (s == state::UNKNOWN)
            throw cado::error(
                    "Error: logarithms of {:x} cannot be reconstructed "
                    "from this set of relations. Abort!\n", h);

        if (s == state::KNOWN_FROM_LOGFILE) {
            /* We know the wanted logarithm from linear algebra, no
             * new relation is necessary */
#if DEBUG >= 1
            fmt::print(stderr, "DEBUG: h = {:x} is known from logfile\n", h);
#endif
            return 0;
        }

        /* this h can be reconstructed */
        size_t nadded = 1;
        needed_rels[relnum] = true;
        auto [ rowsize, row ] = rows[relnum];

#if DEBUG >= 1
        fmt::print(stderr, "DEBUG: h = {:x} can be reconstructed\n", h);
        fmt::print(stderr, "DEBUG:     relation {} added\n", relnum);
        fmt::print(stderr, "DEBUG:     depends of {} others logs\n",
                nb_needed - 1);
#endif

        for (size_t j = 0 ; j < rowsize ; j++) {
            index_t const hh = row[j];
            auto [ nextstate, next_i ] = G[hh];
            ASSERT(nextstate == state::RECONSTRUCTED);
            if (!needed_rels[next_i])
                nadded += needed_rels_from_index(hh, needed_rels);
        }
        return nadded;
    } /* }}} */

    /* Examine relations in [start,end[. Find which logs can be computed
     * from them.
     */
    size_t chunk(std::vector<bool> & not_yet_used, uint64_t start, uint64_t end) /* {{{ */
    {
        size_t computed = 0;

        for (size_t i = start; i < end; i++) {
            if (!not_yet_used[i])
                continue;

            size_t nb = 0;
            index_t h = 0;
            auto [ rowsize, row ] = rows[i];
            for(size_t j = 0 ; j < rowsize ; j++) {
                index_t const hh = row[j];
                if (G[hh].first == state::UNKNOWN) {
                    nb++;
                    h = hh;
                }
            }
            if (nb <= 1) {
                not_yet_used[i] = false;
                if (nb == 1) {
                    G.mark_node_reconstructed(h, i);
                    computed++;
                }
            }
        }
        return computed;
    } /* }}} */

    /* Callback function called by filter_rels in compute_log_from_rels */
    void thread_insert(relation_type & rel) /* {{{ */
    {
        auto pred = [&](prime_type_for_indexed_relations const & p) {
            auto const [ h, e ] = p;
            return G[h].first == state::UNKNOWN;
        };
        auto proj = [&](prime_type_for_indexed_relations const & p) {
            auto const [ h, e ] = p;
            return h;
        };
        new_row(rel.num, rel.primes, pred, proj);
    } /* }}} */

    compute_needed_rels(logtab & log)
        : G(log.nprimes)
    {
        G.set_log_already_known(log);
    }

    /* Given a filename, compute all the relations needed to compute the
     * logarithms appearing in the file.  needed_rels should be
     * initialized before calling this function. Its size must be
     * nrels_purged + nrels_del.
     *
     * Output:
     *    needed_rels, where bits of needed rels are set.
     */
    std::vector<bool> process_input(
            input_relation_pool const & input,
            int nt,
            std::istream& wanted)
    {
        /* Reading all relations */
        fmt::print("# Reading relations from {}\n", input);
        fflush(stdout);
        size_t nrels = input.filter<relation_type>(nullptr,
                [&](relation_type & rel) { thread_insert(rel); });
        fmt::print("# Total: {} relations found\n", nrels);

        std::vector<bool> not_yet_used(nrels, true);

        /* computing dependencies */
        fmt::print("# Starting to compute dependencies from rels\n");

        if (nt > 1)
            fmt::print("# Using multithread version with {} threads\n", nt);
        else
            fmt::print("# Using monothread version\n");

        const double wct_tt0 = wct_seconds();
        for(int iter = 0 ; ; iter++) {
            fmt::print("# Iteration {}: starting...\n", iter);
            fflush(stdout);
            const double wct_tt = wct_seconds();

            /* Compute all missing logarithms possible.
             * Run through all the relations once.
             * Return the number of computed logarithms */ // E: use of undeclared …
            size_t computed = do_one_iter_mt<compute_needed_rels>(not_yet_used, nt);
            fmt::print("# Iteration {}: {} new dependencies computed\n",
                    iter, computed);
            fmt::print("# Iteration {} took {:.1f}s (wall-clock time).\n", iter,
                    wct_seconds() - wct_tt);

            if (!computed)
                break;
        }
        not_yet_used.clear();
        not_yet_used.shrink_to_fit();

        fmt::print("# Computing dependencies took {:.1f}s (wall-clock time)\n",
               wct_seconds() - wct_tt0);

        fmt::print("# Reading wanted logarithms\n");
        fflush(stdout);

        std::vector<bool> needed_rels(nrels, false);
        size_t nwanted = 0;
        size_t nrels_necessary = 0;
        const double wct_tt = wct_seconds();
        wanted >> std::hex;
        for(index_t h ; wanted >> h ; nwanted++) {
            FATAL_ERROR_CHECK(h >= G.size(), "index too large");
            fmt::print("# Computing rels necessary for wanted log {:x}\n", h);
            fflush(stdout);
            auto nadded = needed_rels_from_index(h, needed_rels);
            nrels_necessary += nadded;
            fmt::print("-> {} needed relations were added ({} so far)\n",
                   nadded, nrels_necessary);
        }

        fmt::print("# Reading {} wanted logarithms took {:.1f}s\n", nwanted,
               wct_seconds() - wct_tt);
        fmt::print("# {} relations (out of {}) are needed to compute these logarithms\n",
               nrels_necessary, nrels);
        ASSERT_ALWAYS(nrels_necessary == static_cast<size_t>(std::ranges::count(needed_rels, true)));

        return needed_rels;
    } /* }}} */
};

/********************* usage functions and main ******************************/
static void declare_usage(cxx_param_list & pl)
{
    verbose_decl_usage(pl);
}

struct sm_computer { /* {{{ */
    /* this wraps around the sm_info structure, in a way that exposes
     * some command-line arguments. This functionality should perhaps
     * rather be embedded in the sm_info structure itself.
     */
    cxx_cado_poly const & cpoly;
    cxx_mpz const & ell;
    parameter<std::vector<int>,
        "nsm",
        "number of SM's to add on side 0,1,...">
            nsm;
    parameter<std::string,
        "sm-mode",
        "SM mode (see sm-portability.h)">
            sm_mode;

    static void configure(cxx_param_list & pl) {
        decltype(nsm)::configure(pl);
        decltype(sm_mode)::configure(pl);
    }

    sm_computer(cxx_param_list & pl,
            cxx_cado_poly const & cpoly,
            cxx_mpz const & ell)
        : cpoly(cpoly)
        , ell(ell)
        , nsm(pl)
        , sm_mode(pl)
    {
        /* TODO: parse_per_side has no equivalent in the params
         * declarator model, we have to do it by hand */
        pl.parse_per_side("nsm", nsm.parameter_value(), cpoly.nsides(), -1);
        for (int side = 0; side < cpoly.nsides(); side++) {
            if (nsm.parameter_value()[side] < 0)
                continue;
            if (nsm.parameter_value()[side] > cpoly[side]->deg)
                pl.fail("nsm{}={} can not exceed the degree={}\n",
                        side, nsm.parameter_value()[side], cpoly[side]->deg);
        }
    }

    std::vector<sm_side_info> compute_sm_info() const
    {
        /* Init data for computation of the SMs. */
        std::vector<sm_side_info> sm_info;
        for (int side = 0; side < cpoly.nsides(); side++) {
            sm_info.emplace_back(cpoly[side], ell, 0);
            sm_info[side].set_mode(sm_mode.parameter_value().c_str());
            fmt::print("\n# Polynomial on side {}:\n# F[{}] = {}\n",
                    side, side, cpoly[side]);
            fmt::print("# SM info on side {}:\n", side);
            sm_info[side].print(stdout);
            /* command line wins */
            if (nsm.parameter_value()[side] >= 0)
                sm_info[side].nsm = nsm.parameter_value()[side];
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
        return sm_info;
    }
}; /* }}} */

struct reconstructlog_process { /* {{{ */

    parameter_mandatory<cxx_mpz,
        "ell",
        "group order (see sm -ell parameter)">
            ell;

    parameter_switch<"partial",
        "do not reconstruct everything that can be reconstructed">
            partial;

    parameter<std::string,
        "wanted",
        "file containing list of wanted logs">
            wantedfilename;

    parameter<std::string,
        "renumber",
        "input file for renumbering table">
            renumberfilename;

    parameter_with_default<int,
        "mt",
        "number of threads",
        "1">
            mt;

    parameter_mandatory<std::string,
        "out",
        "output file for logarithms">
            outfilename;

    cxx_cado_poly cpoly;
    renumber_t renumber_table;
    const sm_computer SM;
    const input_relation_pool input;
    logtab log;

    std::vector<sm_side_info> sm_info;

    static void configure(cxx_param_list & pl) {
        pl.declare_usage_section("General operational flags");
        decltype(ell)::configure(pl);
        decltype(partial)::configure(pl);
        decltype(wantedfilename)::configure(pl);
        decltype(renumberfilename)::configure(pl);
        decltype(mt)::configure(pl);
        decltype(outfilename)::configure(pl);

        decltype(SM)::configure(pl);
        decltype(cpoly)::declare_usage(pl);

        decltype(input)::configure(pl);

        decltype(log)::configure(pl);
    }

    explicit reconstructlog_process(cxx_param_list & pl)
        : ell(pl)
        , partial(pl)
        , wantedfilename(pl)
        , renumberfilename(pl)
        , mt(pl)
        , outfilename(pl)
        , cpoly(pl)
        , renumber_table(cpoly)
        , SM(pl, cpoly, ell())
        , input(pl)
        , log(pl, cpoly, ell())
    {
        if (ell() <= 0)
            pl.fail("ell must be >0");
        if (wantedfilename.is_provided() && !partial.is_provided())
            pl.fail("-wanted implies -partial");
        if (mt < 1)
            pl.fail("Error: parameter mt must be at least 1\n");
    }

    void prepare() {
        sm_info = SM.compute_sm_info();

        fmt::print("\n###### Reading renumber file ######\n");
        renumber_table.read_from_file(renumberfilename, true);
        // nprimes = renumber_table.size();

        fmt::print("\n###### Reading known logarithms ######\n");
        fflush(stdout);
        log.read(renumber_table, sm_info);
    }

    void do_main_stuff_partial() {

        ifstream_maybe_compressed wanted(wantedfilename());
        if (!wanted.good())
            throw cado::error("Cannot read {}", wantedfilename());

        auto needed_rels = compute_needed_rels(log)
            .process_input(input, mt, wanted);

        if (std::ranges::none_of(needed_rels, std::identity())) {
            fmt::print("# All wanted logarithms are already known\n");
            return;
        }
        /* Computing logs using rels in purged file */
        compute_log_from_rels(log, cpoly, sm_info, renumber_table)
            .process_input(&needed_rels, input, mt);
    }

    void do_main_stuff_complete() {
        compute_log_from_rels(log, cpoly, sm_info, renumber_table)
            .process_input(nullptr, input, mt);
    }

    void do_main_stuff() {
        if (partial && wantedfilename.is_provided()) {
            do_main_stuff_partial();
        } else if (partial) {
            /* -partial without -wanted can make sense if we want to just
             * write the logs */
            fmt::print("\n# -partial without -wanted: nothing to do\n");
        } else {
            do_main_stuff_complete();
        }

        /* Writing all the logs in outfile */
        fmt::print("\n###### Writing logarithms in a file ######\n");
        log.write(outfilename(), renumber_table, sm_info);

    }
}; /* }}} */

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    declare_usage(pl);
    reconstructlog_process::configure(pl);

    param_list_process_command_line(pl, &argc, &argv, false);

    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    reconstructlog_process R(pl);

    if (pl.warn_unused())
        pl.fail("Unused parameters are given\n");


    R.prepare();

    R.do_main_stuff();
    return EXIT_SUCCESS;
}
