/* dup1: 1st duplicate pass, split relation files into 'nslices'
         slices (adapted from check).

   Usage:
   dup1 [-bz] [-n nslices_log] -out <dir> file1 ... filen
   by default nslices_log = 1 (nslices = 2).

   Files file1 ... filen are split into 'nslices' slices in
   <dir>/0/filej ... <dir>/31/filej.

   If option -bz is given, then the output is compressed with bzip2
   instead of gzip.
   Input can be in gzipped or bzipped format.
*/

#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>

#include <string>
#include <memory>
#include <utility>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "filter_config.h"
#include "filter_io.hpp"
#include "gzip.h"
#include "macros.h"
#include "portability.h" // strdup // IWYU pragma: keep
#include "filelist.hpp"
#include "params.hpp"
#include "timing.h"
#include "verbose.hpp"
#include "fstream_maybe_compressed.hpp"
#include "cxx_mpz.hpp"

/* {{{ set_of_files
 *
 * a set of files is given by a filename pattern, and a max number of lines per
 * file. Each write is done by the call to write_line(p, q). The newline
 * is added automatically.
 */
struct set_of_files {
    std::string pattern;
    std::string message_pattern;
    std::string filename;
    std::unique_ptr<ofstream_maybe_compressed> file;
    unsigned int next_idx = 0;
    size_t lines_per_file;
    size_t lines_left = 0;

    set_of_files(std::string pattern,
            const size_t lines_per_file,
            std::string message_pattern = {})
        : pattern(std::move(pattern))
        , message_pattern(std::move(message_pattern))
        , lines_per_file(lines_per_file)
    {
    }

    private:
    /* Closes the currently open file, if any, and opens the next one */
    void open_next_file()
    {
        filename = fmt::format(fmt::runtime(pattern), next_idx++);
        file = std::make_unique<ofstream_maybe_compressed>(filename);

        if (!message_pattern.empty())
            fmt::print(stderr, fmt::runtime(message_pattern), filename);
        lines_left = lines_per_file;
    }

    public:
    void write_line(const char * p, const char * q)
    {
        if (lines_left == 0)
            open_next_file();
        file->write(p, q-p);
        file->put('\n');
        if (!file->good())
            throw cado::error("Error writing relation to {}", filename);
        lines_left--;
    }
};// }}}

template<typename ab_type, int base>
using dup1_relation =
    cado::relation_building_blocks::line_block<
    cado::relation_building_blocks::ab_block<ab_type, base>>;

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage_section("output specification");
    pl.declare_usage("out", "output directory");
    pl.declare_usage("prefix", "prefix for output files");
    pl.declare_usage("outfmt", "format of output file (default same as input)");
}

struct dup1_process {
    parameter_switch<"abhexa", "read a and b as hexa not decimal"> abhexa;
    parameter_switch<"ab", "only print a and b in the output"> only_ab;
    parameter_switch<"large-ab", "enable support for a,b beyond 64 bits">
        largeab;

    parameter_with_default<unsigned int,
                        "n",
                        "log of number of slices",
                        "1">
                        nslices_log;
    parameter_with_default<unsigned int,
                        "lognrels",
                        "log of number of rels per output file",
                        "25">
                        log_max_nrels_per_files;

    parameter_with_default<int,
                        "only",
                        "do only slice i (-1 means all)",
                        "-1">
                        only_slice;


    /* nslices controls the two arrays below */
    unsigned int nslices = 0;
    std::vector<size_t> nr_rels_tot;
    std::vector<bool> do_slice;
    void update_nslices_log() {
        nslices = 1 << nslices_log;
        nr_rels_tot.assign(1 << nslices_log, 0);
        do_slice.assign(1 << nslices_log, true);
    }

    std::vector<set_of_files> S;

    explicit dup1_process(cxx_param_list & pl)
        : abhexa(pl)
        , only_ab(pl)
        , largeab(pl)
        , nslices_log(pl)
        , log_max_nrels_per_files(pl)
        , only_slice(pl)
    {
        update_nslices_log();

        if (only_slice >= 0) {
            do_slice.assign(nslices, false);
            do_slice[only_slice] = true;
        }
    }

    static void configure(cxx_param_list & pl)
    {
        pl.declare_usage_section("general operational flags");
        decltype(dup1_process::only_slice)::configure(pl);
        decltype(dup1_process::nslices_log)::configure(pl);
        decltype(dup1_process::log_max_nrels_per_files)::configure(pl);
        decltype(dup1_process::only_ab)::configure(pl);
        decltype(dup1_process::abhexa)::configure(pl);
        decltype(dup1_process::largeab)::configure(pl);
    }


    void prepare_filesets(
            std::string const & outdir,
            std::string const & prefix_files,
            std::string const & outfmt)
    {
        S.clear();
        S.reserve(nslices);
        for(unsigned int i = 0 ; i < nslices ; i++) {
            const std::string pattern = fmt::format("{}/{}/{}.{{:04d}}{}{}",
                    outdir, i, prefix_files,
                    only_ab ? ".ab" : "",
                    outfmt);
            /* XXX this comment is significant, and parsed by
             * scripts/cadofactor/cadotask.py
             */
            const std::string message_pattern = fmt::format(
                    "# Opening output file for slice {}: {{}}\n", i);
            S.emplace_back(pattern, 1UL<<log_max_nrels_per_files, message_pattern);
        }
    }

    /* Must be called only when nslices_log > 0 */
    unsigned int compute_slice (int64_t a, uint64_t b) const
    {
        const uint64_t h = CA_DUP1 * static_cast<uint64_t>(a) + CB_DUP1 * b;
        /* Using the low bit of h is not a good idea, since then
           odd values of i are twice more likely. The second low bit
           also gives a small bias with RSA768 (but not for random
           coprime a, b). We use here the nslices_log high bits.
           */
        return static_cast<unsigned int>(h >> (64 - nslices_log));
    }

    unsigned int compute_slice (cxx_mpz const & a, cxx_mpz const & b) const
    {
        const cxx_mpz t = a * CA_DUP1 + b * CA_DUP1;
        mp_limb_t h = 0U;
        ASSERT_ALWAYS(mpz_size(t) <= std::numeric_limits<mp_size_t>::max());
        for (size_t i = 0; i < mpz_size(t); ++i)
            h ^= mpz_getlimbn(t, static_cast<mp_size_t>(i));
        return static_cast<unsigned int>(h >> (GMP_NUMB_BITS - nslices_log));
    }


    template<bool slice0_only, typename relation_type>
    void process(relation_type & rel)
    {
        unsigned int slice = 0;
        if constexpr (!slice0_only) 
            slice = compute_slice (rel.a, rel.b);
        if (!do_slice[slice])
            return;

        const char * p = rel.line.data();
        const char * q = rel.line.data() + rel.line.size();

        /* note that rel.line does not contain the trailing newline. */
        if (only_ab)
            q = std::ranges::find(p, q, ':');
        
        S[slice].write_line(p, q);
        nr_rels_tot[slice]++;
    }

    template<typename relation_type>
    void filter(std::vector<std::string> const & files)
    {
        if (nslices == 1) {
            filter_rels<relation_type>(files, 
                    nullptr, nullptr,
                    [this](relation_type & rel) { process<true>(rel); });
        } else {
            filter_rels<relation_type>(files, 
                    nullptr, nullptr,
                    [this](relation_type & rel) { process<false>(rel); });
        }
    }

    template<typename ab_type>
    void filter0(std::vector<std::string> const & files)
    {
        if (abhexa) {
            filter<dup1_relation<ab_type, 16>>(files);
        } else {
            filter<dup1_relation<ab_type, 10>>(files);
        }
    }

    void filter(std::vector<std::string> const & files)
    {
        if (largeab) {
            filter0<cxx_mpz>(files);
        } else {
            filter0<uint64_t>(files);
        }
        for (unsigned int i = 0; i < nslices; i++) {
            fmt::print (stderr, "# slice {} received {} relations\n", i,
                    nr_rels_tot[i]);
        }
    }

};

int
main (int argc, char const * argv[])
{
    cxx_param_list pl;

    declare_usage(pl);

    filelist::configure(pl);
    dup1_process::configure(pl);
    verbose_decl_usage(pl);
    cado::filter_io_details::configure(pl);


    /* used for counting time in different processes */
    timingstats_dict_t stats;

    param_list_process_command_line(pl, &argc, &argv, true);

    dup1_process D(pl);
    verbose_interpret_parameters(pl);
    cado::filter_io_details::interpret_parameters(pl);

    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    std::string outfmt;
    pl.parse("outfmt", outfmt);

    const char * outdir = param_list_lookup_string(pl, "out");
    const char * prefix_files = param_list_lookup_string(pl, "prefix");

    const filelist input(pl, argc, argv);

    if (param_list_warn_unused(pl))
        pl.fail("Error, unused parameters are given\n");
    if (!prefix_files)
        pl.fail("Error, missing -prefix command line argument\n");
    if (!outdir)
        pl.fail("Error, missing -out command line argument\n");
    if (!outfmt.empty() && !is_supported_compression_format(outfmt.c_str()))
        pl.fail("Error, output compression format unsupported\n");

    const auto input_files = input.create_file_list();

    // If not output suffix is specified, use suffix of first input file
    if (outfmt.empty() && !input_files.empty()) {
        const char * tmp;
        get_suffix_from_filename (input_files[0].c_str(), &tmp);
        if (tmp)
            outfmt = std::string(tmp);
    }

    D.prepare_filesets(outdir, prefix_files, outfmt);

    timingstats_dict_init(stats);

    D.filter(input_files);

    // double thread_times[2];
    // thread_seconds_user_sys(thread_times);
    timingstats_dict_add_mythread(stats, "main");
    // fmt::print(stderr, "Main thread ends after having spent {:.2f}s+{:.2f}s on cpu \n", thread_times[0], thread_times[1]);
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return 0;
}
