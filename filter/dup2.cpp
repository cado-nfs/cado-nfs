/* dup2: 2nd pass

   {{{ some (very certainly outdated) documentation:

   Usage: dup2 -nrels <nrels> -renumber <renumberfile> [ -outdir <dir> ]
               [ -outfmt <fmt> ] [ -dl ]
               [ -filelist <fl> [ -basepath <dir> ] | file1 ... filen ]

   Input files can be given on command line, or via a filelist file.

   By default dup2 overwrites input files. To modify this behaviour an output
   directory can be given with '-outdir'.
   In case a filelist is given, the -basepath option enables to tell in which
   directory those files are.

   By default, the output will be bzipped/gzipped according to the status
   of the input file. The output format can be chosen with '-outfmt'.

   Allocates a hash-table of size next_prime(100 + 1.2*nrels).

   The relations are renumbered according to the file given via the
   '-renumberfile' argument. Input relations are of the following format:

       a,b[@s1,s2]:p1,p2,...,pj:q1,q2,...,qk                         (*)

   where p1,p2,...,pj are ideals on side s1 (possibly duplicate), and
   q1,q2,...,qk are ideals on side s2 (possibly duplicate). The part '@s1,s2'
   is optional and will default to s1=0 and s2=1 if not present. Output is:

       a,b:r1,r2,...,rm                                              (**)

   where each index r1,r2,...,rm refering to ideals via the renumbering table.
   By default valuations are reduced modulo 2. This can be changed by using the
   '-dl' argument.

   The format of each file is recognized by counting the number of ':' in the
   first line: if two we have the raw format (*), if only one we have the
   renumbered format (**). It is assumed that all files in renumbered format
   come first.

   Algorithm: for each (a,b) pair, we compute h(a,b) = (CA*a+CB*b) % 2^64.

   h has 64 bits. We know bits 2..6 are identical for a given slice, thus
   we remove them, and we obtain a 59-bit value h'.

   Let h' = j * 2^k + i.

   We store j % 2^32 at the next empty cell after index i in the hash-table.
  
   }}}
*/

#include "cado.h"     // IWYU pragma: keep
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <ostream>

#include <memory>

#include <gmp.h>
#include "fmt/base.h"

#include "cado_poly.hpp"
#include "filter_config.h"
#include "filter_io.hpp"
#include "gmp_aux.h"
#include "gzip.h"
#include "macros.h"
#include "filelist.hpp"
#include "params.hpp"
#include "portability.h"
#include "relation-tools.h"
#include "renumber.hpp"
#include "misc.h"
#include "typedefs.h"
#include "verbose.hpp"
#include "cxx_mpz.hpp"
#include "utils_cxx.hpp"
#include "fstream_maybe_compressed.hpp"

#define xxxDEBUG

/* For debugging */
// #define TRACE_HASH_TABLE
#ifdef TRACE_HASH_TABLE
#define TRACE_I 42
#define TRACE_J 17
#endif

#include <variant>


struct missing_parameter : public parameter_error { /* {{{ */
    explicit missing_parameter(std::string const & what)
        : parameter_error("missing parameter " + what)
    {}
}; /* }}} */

/* {{{ check_whether_file_is_renumbered
 *
 * Probe the input file to see if it looks like a file that has
 * already undergone dup2 renumbering. It's always a bit shaggy. Here
 * we count the colons in the first few lines, but we can certainly
 * imagine scenarios where this would fail (what if we have just one
 * side, for example?).
 */
bool check_whether_file_is_renumbered(std::string const & filename)
{
    unsigned int count = 0;
    char s[1024];
    FILE * f_tmp = fopen_maybe_compressed(filename.c_str(), "rb");

    if (!f_tmp)
        throw cado::error("{}: {}", filename, strerror(errno));

    if (feof(f_tmp)) /* file is empty */
    {
        fclose_maybe_compressed(f_tmp, filename.c_str());
        return true; /* an empty file might be considered as renumbered */
    }

    /* Look for first non-comment line */
    for(;;) {
        const char * ret = fgets(s, 1024, f_tmp);
        if (ret == NULL) {
            /* fgets returns NULL when the end of file occurs or when there is
               an error */
            if (feof(f_tmp)) {
                fclose_maybe_compressed(f_tmp, filename.c_str());
                return true;
            }
            throw cado::error("Error while reading {}", filename);
        }
        if (strlen(s) >= 1023)
            throw cado::error("Too long line while reading {}", filename);
        size_t i = 0;
        while (s[i] == ' ')
            i++;
        if (s[i] != '#')
            break;
    }
    for (unsigned int i = 0; i < strlen(s); i++)
        count += s[i] == ':';
    fclose_maybe_compressed(f_tmp, filename.c_str());

    if (count == 1)
        return true;
    else if (count == 2)
        return false;
    else
        throw cado::error("Invalid line in {}"
                    " (line has {} colons but 1 or 2 where expected:\n{}",
                    filename, count, s);
}
/* }}} */

struct output_specification { /* {{{ */
    std::string outfmt;
    std::string outdir;

    static void declare_usage(cxx_param_list & pl) {
        pl.declare_usage("outdir", "by default, input files are overwritten");
        pl.declare_usage("outfmt", "format of output file (default same as input)");
    }

    output_specification() = default;

    explicit output_specification(cxx_param_list & pl)
    {
        pl.parse("outdir", outdir);
        pl.parse("outfmt", outfmt);

        if (!outfmt.empty())
            if (!is_supported_compression_format(outfmt.c_str()))
                throw parameter_error("output compression format not supported");
    }
    /* return in {oname, oname_tmp} two file names for writing the output
     * of processing the given input file f. Both files are placed in the
     * directory outdir if not NULL, otherwise in the current directory.
     * The parameter outfmt specifies the output file extension and
     * format (semantics are as for fopen_maybe_compressed).
     *
     * proper use requires that data be first written to the file whose name
     * is *oname_tmp, and later on upon successful completion, that file must
     * be renamed to *oname. Otherwise disaster may occur, as there is a slim
     * possibility that *oname == infilename on return.
     */
    std::pair<std::string, std::string>
    get_outfilename_from_infilename(std::string const & f) const
    {
        const std::string suffix_in = get_suffix(f, "");
        const std::string suffix_out = outfmt.empty() ? suffix_in : outfmt;
        std::string newname(f.begin(), f.end() - suffix_in.size());

        std::pair<std::string, std::string> res;

        if (!outdir.empty()) {
            const std::string basename(path_basename(newname.c_str()));
            newname = outdir + "/" + basename;
        }
        res = {
            newname + suffix_out,
            newname + ".tmp" + suffix_out,
        };

#ifdef DEBUG
        fmt::print(stderr,
                "DEBUG: Input file name: {}\n"
                "DEBUG: temporary output file name: {}\n"
                "DEBUG: final output file name: {}\n",
                f, res.second, res.first);
#endif

        return res;
    }
}; /* }}} */

struct dup2_process { /* {{{ */
    output_specification output;

    std::string polyfilename;
    cxx_cado_poly cpoly;

    /* Renumbering table to convert from (p,r) to an index */
    std::string renumberfilename;
    renumber_t renumber_tab;

    bool is_for_dl = false; /* By default we do dup2 for factorization */
    bool largeab = false; /* By default, do not use mpz for a,b */

    explicit dup2_process(cxx_param_list & pl)
        : output(pl)
        , polyfilename(pl.parse_mandatory<std::string>("poly")) // E: Member in…
        , renumberfilename(pl.parse_mandatory<std::string>("renumber"))
        , nrels_expected(pl.parse_mandatory<size_t>("nrels"))
    {
        /* those are switches, but declared with a nullptr target, so
         * that they can be parsed later. Which means now.
         */
        pl.parse("dl", is_for_dl);
        pl.parse("large-ab", largeab);

    }

    static void declare_usage(cxx_param_list & pl) {
        pl.declare_usage("poly", "input polynomial file");
        pl.declare_usage("renumber", "input file for renumbering table");
        pl.declare_usage("nrels", "number of expected relations");
        pl.declare_usage("dl", "do not reduce exponents modulo 2");
        pl.declare_usage("large-ab", "enable support for a and b beyond 64 bits");
    }

    static void configure_switches(cxx_param_list & pl) {
        pl.configure_switch("dl");
        pl.configure_switch("large-ab");
    }

    void read()
    {
        /* This is sort of a "late" ctor.  */

        if (!cpoly.read(polyfilename))
            throw cado::error("cannot read {}", polyfilename); // E: Use of undeclared identifier 'polyfilename'

        renumber_tab = renumber_t(cpoly);
        renumber_tab.read_from_file(renumberfilename, is_for_dl);

        /* determine the expected size of the hash table, and allocate */

        K = 100 + nrels_expected + (nrels_expected / 5);
        K = uint64_nextprime(K);

        H = std::make_unique<uint32_t[]>(K);
        std::fill_n(H.get(), K, 0);
        fmt::print(stderr, "Allocated hash table of {} entries ({})\n",
                K, size_disp(K * sizeof(uint32_t)));

        /* sanity check: since we allocate two 64-bit words for each, instead of
           one 32-bit word for the hash table, taking K/100 will use 2.5% extra
           memory */
        sanity_ab.resize(1U + (K / 100U));
        fmt::print(stderr, "[checking true duplicates on sample of {} cells]\n",
                sanity_ab.size());

    }


    std::unique_ptr<uint32_t[]> H;   /* H contains the hash table */
    size_t K = 0; /* Size of the hash table */
    unsigned long nrels_expected = 0;

    double cost = 0.0; /* Cost to insert all rels in the hash table */

    /* Number of duplicates and rels on the current file */
    size_t ndup = 0;
    size_t nrels = 0;

    /* Number of duplicates and rels on all read files */
    size_t ndup_tot = 0;
    size_t nrels_tot = 0;
    size_t nrels_already_renumbered = 0;

    /* {{{ sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
       and check for hash collisions */
    std::vector<std::pair<cxx_mpz, cxx_mpz>> sanity_ab;
    unsigned long sanity_checked = 0;
    unsigned long sanity_collisions = 0;
    template <typename Ta, typename Tb>
    void sanity_check(uint64_t i, Ta a, Tb b)
    {
        if (i >= sanity_ab.size())
            return;
        sanity_checked++;
        if (sanity_ab[i].first == 0) {
            sanity_ab[i].first = a;
            sanity_ab[i].second = b;
        } else if (sanity_ab[i].first != a) {
            sanity_collisions++;
            fmt::print(stderr,
                    "Collision between ({}, {}) and ({}, {})\n",
                    sanity_ab[i].first, sanity_ab[i].second,
                    cxx_mpz(a), cxx_mpz(b));
        }
    }
    /* }}} end sanity check */

    /* {{{ periodically check for overfilling */
    double factor = 1;
    void check_overfilling()
    {
        if (cost >= factor * (double)(nrels_tot - ndup_tot)) {
            const uint64_t nodup = nrels_tot - ndup_tot;
            const double full_table = 100.0 * double_ratio(nodup, K);
            fmt::print(stderr,
                    "Warning, hash table is {:1.0f}%"
                    " full (avg cost {:1.2f})\n",
                    full_table, double_ratio(cost, nrels_tot));
            if (full_table >= 99) {
                fprintf(stderr, "Error, hash table is full\n");
                exit(1);
            }
            factor += 1.0;
        }
    }
    /* }}} */

    /* The number of threads used to compute roots mod p is hard-coded to 3,
       since it was observed that using a large value gives some contention.
       The value of 3 is optimal for a c130 on a 64-core node. */
    int nthreads_for_roots = 3;

    /* it isn't totally clear if we need (and even if we can) increase
     * this */
    int nthreads_hash = 1;

    /* {{{ compute_hash (two overloads) */
    static uint64_t compute_hash(int64_t a, uint64_t b)
    {
        return CA_DUP2 * (uint64_t)a + CB_DUP2 * b;
    }

    static uint64_t compute_hash(cxx_mpz const & a, cxx_mpz const & b)
    {
        cxx_mpz r = a * CA_DUP2 + b * CB_DUP2;
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }
    /* }}} */

    /* {{{ insert_relation_in_dup_hashtable
     *
     * if duplicate is_dup = 1, else is_dup = 0
     * return i for sanity check, which is the place in the hash table where
     * the relation is stored: more precisely we store in H[i] the value
     * of floor(h/2^32), where h(a,b) is a 64-bit value.
     */
    template<typename relation_type>
    size_t insert_relation_in_dup_hashtable(relation_type & rel, bool & is_dup)
    {
        auto h = compute_hash(rel.a, rel.b);

        /* We put j = floor(h/2^32) in cell H[i] where i = h % K (or the first
           available cell after i if H[i] is already occupied). We have extraneous
           duplicates when:
           (a) either we have a collision on h: this should happen with probability
           2^(-64)
           (b) j = j' and |i-i'| is small, in particular i=i'. If K is at least
           twice the number of relations, then |i-i'| is bounded by say 2 on
           average, thus this happens when h' = h + n*K or h' = h + n*K + 1. This
           happens with probability 2/K. Moreover we should have j = j', which
           happens with probability 2^(-32). Since K is an odd prime, the global
           probability is about 2^(-31)/K. */

        uint64_t i = h % K;
        auto j = (uint32_t)(h >> 32);

#ifdef TRACE_HASH_TABLE
        uint64_t old_i = i;
#endif

        double local_cost = 0;
        while (H[i] != 0 && H[i] != j) {
            i++;
            if (UNLIKELY(i == K))
                i = 0;
            local_cost++;
        }

        /* The largest block for a linear probing table of size m with l entries
         * behaves like (log(m)-1.5*log(log(m)))/(beta-1-log(beta)) where beta =
         * l/m. See B. Pittel, Linear probing: the probable largest search time
         * grows logarithmically with the number of records, Journal of Algorithms
         * 8, number 2, 236-249, 1987. Here we have m = K and l <= nrels_expected.
         * Since we take K >= 6/5 * nrels_expected, we have beta <= 5/6,
         * the factor 1/(beta-1-log(beta)) is bounded by 64.
         * For m=2e8, the largest block can have up to length 938,
         * thus it is not surprising to have local_cost > 100. */

        cost += local_cost;

        /* Note: since we use 0 for uninitialized entries, entries with
         * j=0 will get always marked as 'duplicate' and be lost. */

        is_dup = H[i] == j;
        H[i] = j;

#ifdef TRACE_HASH_TABLE
        if (i == TRACE_I && j == TRACE_J) {
            fmt::print(stderr,
                    "TRACE: a = {}\n"
                    "TRACE: b = {}\n"
                    "TRACE: i = {}\n"
                    "TRACE: j = {}\n"
                    "TRACE: initial value of i was {}\n"
                    "TRACE: h = {}\n"
                    "TRACE: is_dup = {}\n",
                    rel.a, rel.b, i, j, old_i, h, is_dup);
        }
#endif
        return i;
    }
    /* }}} */

    template<typename relation_type>
    void hash_renumbered_relation(relation_type & rel)
    {
        bool is_dup;

        nrels++;
        nrels_tot++;
        uint64_t i = insert_relation_in_dup_hashtable(rel, is_dup);

        static unsigned long count = 0;

        // They should be no duplicate in already renumbered file
        if (is_dup && count++ < 10) {
            fmt::print(stderr,
                    "Warning, duplicate relation in already renumbered files:"
                    "\na = {}\nb = {}\ni = {}\nj = {}\n"
                    "This warning may be due to a collision on the hash"
                    " function or to an actual duplicate\n"
                    "relation. If it appears often you should check the"
                    " input set of relations.\n\n",
                    rel.a, rel.b, i, H[i]);
        }

        sanity_check(i, rel.a, rel.b);
        check_overfilling();
    }

    template<typename ab_type> void filter_old(std::vector<std::string> const & files) /* {{{ */
    {
        /* for "old" relation, (a,b) are always in hex.
         *
         * However we probably need to templatize based on the exponent
         * type or things of that sort. */
        using prime_type = prime_type_for_indexed_relations;

        using relation_type = 
            cado::relation_building_blocks::primes_block<prime_type,
            cado::relation_building_blocks::ab_block<ab_type, 16>>;

        filter_rels<relation_type>(files, nullptr, nullptr,
                [this](relation_type & rel) { hash_renumbered_relation(rel); });
    } /* }}} */

    void filter_already_renumbered_rels(std::vector<std::string> const & files) /* {{{ */
    {
        fmt::print(stderr, "Reading {} files already renumbered:\n", files.size());
        if (largeab)
            filter_old<cxx_mpz>(files);
        else
            filter_old<uint64_t>(files);
        nrels_already_renumbered = nrels_tot;
    } /* }}} */

    /* compute_index_rel {{{
     *
     * modify in place the relation rel to take into account:
     *  - the renumbering
     *  - the bad ideals
     * 
     * It is a non-trivial operation. We need to go from
     *  (prime number, side, exponent)
     * to
     *  (prime number, root, side, exponent) (using the information on a,b)
     * and finally to
     *  (prime ideal index, exponent) (using the renumber table)
     *
     */

    template<typename relation_type>
    void compute_index_rel(relation_type & rel)
    requires requires
        {
            /* we want relation_type to be a variant, one with a sieve
             * relation, one with an indexed relation */
            std::variant_size_v<relation_type>;
        }
    {
        /* We used to reduce exponents early, here. However, it's not a
         * good idea if we want to get the valuations at bad ideals.
         * (maybe we could reduce exponents for primes that we know are
         * above the bad ideal bound...).
         *
         * Anyway. Reduction is now done at the _end_ of
         * compute_index_rel.
         */
        using sieve_rel_type = std::variant_alternative_t<0, relation_type>;
        using indexed_rel_type = std::variant_alternative_t<1, relation_type>;

        ASSERT_ALWAYS(std::holds_alternative<sieve_rel_type>(rel));

        auto & srel = std::get<0>(rel);

        uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */

        indexed_rel_type irel;

        using prime_type = indexed_rel_type::prime_type;

        for(auto const & pse : srel.primes) {
            auto p = pse.p;
            auto side = pse.side;
            auto e = pse.e;

            p_r_values_t r = 0;

            if (e == 0) continue;

            if (side != renumber_tab.get_rational_side()) {
#ifdef DEBUG
                // Check for this bug : [#15897] [las] output
                // "ideals" that are not prime
                /* p_r_values_t may be 32-bit... */
                const unsigned long p_ul = p;
                if (!modul_isprime(&p_ul)) {
                    throw std::runtime_error(
                            fmt::format(
                                "Error, relation with a={} b={} contains {}"
                                " which is not prime.\n"
                                "Remove this relation from the file and re-run dup2.\n",
                                srel.a, srel.b, p));
                }
#endif
                r = relation_compute_r(srel.a, srel.b, p);
            } else {
                r = 0; // on the rational side we need not compute r,
                       // which is m mod p.
            }

            renumber_t::p_r_side const prside { p, r, side };

            nonvoidside |= ((uint64_t)1) << side;

            if (renumber_tab.is_bad(prside)) {
                auto [ first_index, exps ] = renumber_tab.indices_from_p_a_b(
                        prside, e, srel.a, srel.b);

                for(size_t i = 0 ; i < exps.size() ; i++) {
                    irel.primes.push_back(prime_type {
                            .h = index_t(first_index + i),
                            .e = exps[i]});
                }
            } else {
                auto h = renumber_tab.index_from_p_r(prside);
                int const f = renumber_tab.inertia_from_p_r(prside);
                if (f > 1) {
                    /* XXX there's a catch here. A non-bad ideal can
                     * still have non-trivial inertia (say f=2), in which
                     * case we must divide the valuation (which comes
                     * from the norm) by the inertia in order to obtain
                     * the valuation at the prime ideal
                     */
                    ASSERT_ALWAYS(e % f == 0);
                    e /= f;
                }
                irel.primes.push_back(prime_type { .h = h, .e = e });
            }
        }

        /* If needed, print the additional columns (they are always at
         * the beginning of the renumbering table).
         * if naddcols == 0:
         *    do nothing.
         * else:
         *    we add the columns i if and only if the polynomial on side
         *    i is non monic and the relation contains at least one prime
         *    on side i.
         *
         * nb_polys==2 is special (see also  renumber_t::index_from_p_r)
         * ; in that case, we only use a single combined additional
         * column. (except for free relations, of course).
         */

        if (renumber_tab.number_of_additional_columns()) {
            size_t n = renumber_tab.get_nb_polys();
            if (n == 2) {
                /* Possible cases when we have two sides:
                 *  - both are monic, there is no J ideal
                 *  - only one is monic, there is only one J ideal
                 *  - neither is monic, we have two J ideals, but since n==2
                 *  it's okay to have only one additional column to count both
                 *  of them together.
                 * Either way, we only have one additional column (if we have
                 * any, of course), and it's number zero
                 */
                irel.primes.emplace_back(prime_type { .h = 0, .e = 1 });
            } else {
                auto sides = renumber_tab.get_sides_of_additional_columns();
                for (size_t idx = 0; idx < sides.size(); idx++) {
                    int side = sides[idx];
                    if ((nonvoidside & (((uint64_t)1) << side))) {
                        irel.primes.emplace_back(
                                prime_type { .h = static_cast<index_t>(idx),
                                             .e = 1 });
                    }
                }
            }
        }

        if (!is_for_dl) { /* Do we reduce mod 2 */
            auto jt = irel.primes.begin();
            for(auto & pe : irel.primes) {
                pe.e &= 1;
                if (pe.e)
                    *jt++ = pe;
            }
            irel.primes.erase(jt, irel.primes.end());
        }

        /* replace the sieved relation by the indexed relation */
        irel.a = std::move(srel.a);
        irel.b = std::move(srel.b);
        rel = std::move(irel);
    } /* }}} */

    template<typename relation_type> void thread_dup2(std::ostream& out, relation_type & rel)
    {
        using indexed_rel_type = std::variant_alternative_t<1, relation_type>;
        ASSERT_ALWAYS(std::holds_alternative<indexed_rel_type>(rel));
        auto & irel = std::get<1>(rel);

        bool is_dup;
        nrels++;
        nrels_tot++;
        const uint64_t i = insert_relation_in_dup_hashtable(irel, is_dup);

        if (!is_dup) {
            sanity_check(i, irel.a, irel.b);
            check_overfilling();
            fmt::print(out, "{}\n", irel);
        } else {
            ndup++;
            ndup_tot++;
        }
    }

    void dup_print_stat(char const * s) const
    {
        fmt::print(stderr, "{}: nrels={} dup={} ({:.2f}%) rem={}\n",
                s, nrels, ndup,
                100.0 * double_ratio(ndup, nrels),
                nrels - ndup);
        fmt::print(stderr, "Total so far: nrels={} dup={} ({:.2f}%) rem={}\n",
                nrels_tot, ndup_tot,
                100.0 * double_ratio(ndup_tot, nrels_tot),
                nrels_tot - ndup_tot);
    }

    template<typename ab_type>
    void filter_new(std::vector<std::string> const & files)
    {
        /* TODO: tnfs and such! */

        static constexpr int base = 10; /* FIXME for abhexa */
        using relation_type = 
            std::variant<
                cado::relation_building_blocks::primes_block<
                    prime_type_for_sieve_relations,
                    cado::relation_building_blocks::ab_block<ab_type, base>
                >,
                cado::relation_building_blocks::primes_block<
                    prime_type_for_indexed_relations,
                    cado::relation_building_blocks::ab_block<ab_type, 16>
                >
            >;
        fmt::print(stderr, "Reading new files (using {} auxiliary threads for "
                "roots mod p):\n", nthreads_for_roots);
        for(auto const & f : files) {
            auto [ oname, oname_tmp ] = output.get_outfilename_from_infilename(f);
            ofstream_maybe_compressed out(oname_tmp);
            nrels = ndup = 0;

            using cado::filter_io_details::multithreaded_call;

            filter_rels<relation_type>(f, nullptr, nullptr,
                    multithreaded_call(nthreads_for_roots,
                        [this](relation_type & rel) {
                            compute_index_rel(rel);
                        }),
                    multithreaded_call(nthreads_hash,
                        [&, this](relation_type & rel) {
                            thread_dup2(out, rel);
                        })
             );

            out.close();
            int rc;
            rc = rename(oname_tmp.c_str(), oname.c_str());
            if (rc < 0)
                throw cado::error("rename({} -> {}): {}",
                            oname_tmp, oname, strerror(errno));

            // stat for the current file
            // + for all the files already read
            dup_print_stat(path_basename(f.c_str()));
        }
    }

    void filter_new_rels(std::vector<std::string> const & files) /* {{{ */
    {
        fmt::print(stderr, "Reading {} new files:\n", files.size());
        if (largeab)
            filter_new<cxx_mpz>(files);
        else
            filter_new<uint64_t>(files);
    } /* }}} */

    void final_stats() { /* {{{ */
        fmt::print(stderr, "At the end: {} remaining relations\n",
                nrels_tot - ndup_tot);

        fmt::print(stderr,
                "At the end: hash table is {:1.2f}% full\n"
                "            hash table cost: {:1.2f} per relation\n",
                100.0 * double_ratio(nrels_tot - ndup_tot, K),
                1.0 + double_ratio(cost, nrels_tot));
        fmt::print(stderr,
                "  [found {} true duplicates on sample of {} relations]\n",
                sanity_collisions, sanity_checked);

        if (nrels_already_renumbered == 0) {
            if (nrels_tot != nrels_expected) {
                fmt::print(stderr,
                        "Warning: number of relations read ({}) does not match"
                        " the number of relations expected ({})\n",
                        nrels_tot, nrels_expected);
            }
        } else {
            /* when we have renumbered files, we know that we won't have the
             * total number of relations... */
            if (nrels_tot > nrels_expected) {
                fmt::print(stderr,
                        "Warning: number of relations read ({}) exceeds"
                        " the number of relations expected ({})\n",
                        nrels_tot, nrels_expected);
            }
        }
    } /* }}} */
}; /* }}} */

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    verbose_decl_usage(pl);
    dup2_process::declare_usage(pl);
    filelist::configure(pl);
    output_specification::declare_usage(pl);
    cado::filter_io_details::configure(pl);

    dup2_process::configure_switches(pl);

    param_list_process_command_line(pl, &argc, &argv, true);

    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    cado::filter_io_details::interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    const filelist input(pl, argc, argv);

    dup2_process D(pl);

    /* we're done looking up parameters */
    if (param_list_warn_unused(pl))
        pl.fail("Error, unused parameters are given\n");

    D.read();
    auto [ oldlist, newlist ] = input.separate_file_list(check_whether_file_is_renumbered);
    D.filter_already_renumbered_rels(oldlist);
    D.filter_new_rels(newlist);
    D.final_stats();


    return 0;
}
