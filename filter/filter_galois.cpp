#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <algorithm>
#include <atomic>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <mutex>
#include <utility>

#ifdef HAVE_MINGW
#include <fcntl.h>
#endif

#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "cado_poly.hpp"
#include "filter_config.h"
#include "filter_io.hpp"
#include "fmt/base.h"
#include "galois_action.hpp"
#include "gzip.h"
#include "fstream_maybe_compressed.hpp"
#include "macros.h"
#include "params.hpp"
#include "portability.h"
#include "relation-tools.h"
#include "renumber.hpp"
#include "timing.h"
#include "typedefs.h"
#include "verbose.hpp"
#include "filelist.hpp"

/* exactly the same as in dup2... */
struct output_specification {
    parameter<std::string,
        "outdir",
        "store output files in another place instead of overwriting inputs">
            outdir;

    parameter<std::string,
        "outfmt",
        "format of output file (default same as input)">
            outfmt;

    static void configure(cxx_param_list & pl) {
        decltype(outdir)::configure(pl);
        decltype(outfmt)::configure(pl);
    }

    output_specification() = default;

    explicit output_specification(cxx_param_list & pl)
        : outdir(pl)
        , outfmt(pl)
    {
        if (!outfmt().empty())
            if (!is_supported_compression_format(outfmt().c_str()))
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
        const std::string suffix_out = outfmt.is_provided() ? outfmt() : suffix_in;
        std::string newname(f.begin(), f.end() - suffix_in.size()); // W: Narro…

        std::pair<std::string, std::string> res;

        if (outdir.is_provided()) {
            const std::string basename(path_basename(newname.c_str()));
            newname = outdir() + "/" + basename;
        }
        res = {
            newname + suffix_out,
            newname + ".tmp" + suffix_out,
        };

#if DEBUG >= 1
        fmt::print(stderr,
                "DEBUG: Input file name: {}\n"
                "DEBUG: temporary output file name: {}\n"
                "DEBUG: final output file name: {}\n",
                f, res.second, res.first);
#endif

        return res;
    }

};


struct filter_galois_process {
    output_specification output;

    parameter_mandatory<std::string,
        "poly",
        "input polynomial file">
            polyfilename;

    parameter_mandatory<std::string,
        "renumber",
        "input file for renumbering table">
            renumberfilename;

    parameter_mandatory<size_t,
        "nrels",
        "(approximate) number of input relations">
            nrels_expected;

    parameter_mandatory<galois_action,
        "galois",
        "Galois action among 1/y or _y">
            action;

    parameter_with_default<int,
        "t",
        "number of threads",
        "1">
            nthreads;

    parameter_switch<
        "dl",
        "for DL (untested)">
            is_for_dl;

    parameter_switch<
        "large-ab",
        "enable support for a and b beyond 64 bits">
            largeab;

    cxx_cado_poly cpoly;

    /* Renumbering table to convert from (p,r) to an index */
    renumber_t renumber_tab;


    /* created by galois_action_on_ideals() */
    std::vector<index_t> ga_id_cache;


    /* see alloc_hashtable */
    std::unique_ptr<uint32_t[]> H;   /* H contains the hash table */
    size_t K = 0; /* Size of the hash table */

    mutable std::mutex io_lock;

    std::atomic<size_t> ndups = 0;

    static void configure(cxx_param_list & pl) {
        decltype(output)::configure(pl);
        decltype(polyfilename)::configure(pl);
        decltype(renumberfilename)::configure(pl);
        decltype(nrels_expected)::configure(pl);
        decltype(action)::configure(pl);
        decltype(nthreads)::configure(pl);
        decltype(is_for_dl)::configure(pl);
        decltype(largeab)::configure(pl);
    }

    explicit filter_galois_process(cxx_param_list & pl)
        : output(pl)
        , polyfilename(pl)
        , renumberfilename(pl)
        , nrels_expected(pl)
        , action(pl)
        , nthreads(pl)
        , is_for_dl(pl)
        , largeab(pl)
    {
        if (!cpoly.read(polyfilename))
            throw cado::error("cannot read {}", polyfilename());
    }

    void read() {
        fmt::print("# Using {}\n", action());
        renumber_tab = renumber_t(cpoly);
        renumber_tab.read_from_file(renumberfilename, is_for_dl);

        if (action().get_order() > 2) {
            std::cerr << "Error, Galois action of order > 2 are not supported "
                         "yet. The missing piece of code is the one that takes "
                         "care of computing the new valuation of the the "
                         "rewritten ideals.\n";
            exit(EXIT_FAILURE);
        }


        if (renumber_tab.number_of_bad_ideals() > 0 && action().get_order() > 1)
            std::cout << "\n/!\\/!\\/!\\/!\\\nWARNING, bad ideals will be left "
                         "unchanged, the output may not be usable depending on "
                         "your use case\nSee comments in utils/galois_action.cpp "
                         "for more info\n/!\\/!\\/!\\/!\\\n\n";
    }

    void alloc_hashtable() {
        K = 100 + 1.2 * double(nrels_expected);
        H = std::unique_ptr<uint32_t[]>(new uint32_t[K]);
        ASSERT_ALWAYS(H);
        std::fill(H.get(), H.get() + K, 0);
    }


    void galois_action_on_ideals()
    {
        std::cout << "Computing Galois action on ideals\n";
        size_t norb = action().compute_action_on_index(ga_id_cache, renumber_tab);
        fmt::print("Found {} orbits of length {},"
                " {} columns were left unchanged "
                "(among which {} additional column(s) and {} column(s) "
                "corresponding to badideals)\n",
                norb, action().get_order(),
                ga_id_cache.size() - norb * action().get_order(),
                renumber_tab.number_of_additional_columns(),
                renumber_tab.number_of_bad_ideals());
    }

    template <typename relation_type>
    size_t insert_relation_in_dup_hashtable(relation_type & rel, bool & is_dup)
    {
        auto h = action().hash_ab(rel.a, rel.b, CA_DUP2, CB_DUP2);

        uint64_t i = h % K;
        auto j = (uint32_t)(h >> 32);
        while (H[i] != 0 && H[i] != j) {
            i++;
            if (UNLIKELY(i == K))
                i = 0;
        }

        is_dup = H[i] == j;
        H[i] = j;
        return i;
    }

    template<typename relation_type>
    void thread_galois(std::ostream & out, relation_type & rel)
    {
        bool is_dup;
        std::vector<index_t> const & sigma = ga_id_cache;
        insert_relation_in_dup_hashtable(rel, is_dup);

        if (is_dup) {
            ndups++;
            return;
        }

        /* transform primes via the Galois action if needed */
        for(auto & [h, e] : rel.primes) {
            auto hrep = sigma[h];
            /* since we have only order-two actions, the action on the
             * exponent is fairly easy to write.
             */
            if (hrep != h && is_for_dl)
                e = -e;
            h = hrep;
        }
        fmt::print(out, "{}\n", rel);
    }

    /* This is the main entry point. It returns the number of
     * non-duplicate relations found in the input file set.
     */

    size_t filter(std::vector<std::string> const & files)
    {
        if (!largeab) {
            using R = cado::relation_building_blocks::primes_block<
                prime_type_for_indexed_relations,
                cado::relation_building_blocks::ab_block<uint64_t, 16>>;
            return filter<R>(files);
        } else {
            using R = cado::relation_building_blocks::primes_block<
                prime_type_for_indexed_relations,
                cado::relation_building_blocks::ab_block<cxx_mpz, 16>>;
            return filter<R>(files);
        }
    }
    private:

    template<typename relation_type>
    size_t filter(std::vector<std::string> const & files) {
        using R = relation_type;
        size_t n = 0;
        for(auto const & f : files) {
            auto [ oname, oname_tmp ] = output.get_outfilename_from_infilename(f);
            ofstream_maybe_compressed out(oname_tmp);

            n += filter_rels<R>(f, nullptr, nullptr,
                    cado::filter_io_details::multithreaded_call(nthreads,
                        [&](R & rel) {
                    thread_galois(out, rel);
                    }));

            out.close();
            int rc;
            rc = rename(oname_tmp.c_str(), oname.c_str());
            if (rc < 0)
                throw cado::error("rename({} -> {}): {}",
                            oname_tmp, oname, strerror(errno));

        }
        return n - ndups;
    }
};

static void declare_usage(cxx_param_list & pl)
{
    verbose_decl_usage(pl);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    declare_usage(pl);
    filelist::configure(pl);
    filter_galois_process::configure(pl);
    cado::filter_io_details::configure(pl);

    param_list_process_command_line(pl, &argc, &argv, true);

    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    cado::filter_io_details::interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    const filelist input(pl, argc, argv);

    filter_galois_process fg(pl);

    if (param_list_warn_unused(pl))
        pl.fail("Error, unused parameters are given\n");

    /* Renumbering table to convert from (p,r) to an index */
    fg.read();
    fg.alloc_hashtable();

    fg.galois_action_on_ideals();

    std::cout << "Rewriting relations files\n";


    std::cout << "Reading files (using 1 auxiliary thread):\n";
    timingstats_dict_t stats;
    timingstats_dict_init(stats);

    size_t noutrels = fg.filter(input.create_file_list());

    /* XXX This printout is parsed by cado-nfs.py */
    fmt::print(stderr, "Number of output relations: {}\n", noutrels);

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);
    return 0;
}
