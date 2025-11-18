#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <climits>
#include <cstdint>
#include <cstdlib>

#include <istream>
#include <fstream>
#include <ios>
#include <string>
#include <utility>
#include <regex>
#include <sstream>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include <unistd.h>
#include <gmp.h>

#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "lingen_average_matsize.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_checkpoints.hpp"
#include "lingen_io_matpoly.hpp"
#include "lingen_matpoly_select.hpp"
#include "logline.hpp"
#include "macros.h"
#include "misc.h"
#include "params.h"
#include "runtime_numeric_cast.hpp"
#include "select_mpi.h"
#include "tree_stats.hpp"
#include "utils_cxx.hpp"

/* Checkpoints */

template<bool is_binary>
std::string lingen_checkpoint<is_binary>::default_directory;
template<bool is_binary>
unsigned int lingen_checkpoint<is_binary>::threshold = UINT_MAX;
template<bool is_binary>
int lingen_checkpoint<is_binary>::save_gathered = 0;

/* There's much copy-paste here */

template<bool is_binary>
void lingen_checkpoint<is_binary>::decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "checkpoint-directory",
            "where to save checkpoints");
    param_list_decl_usage(pl, "checkpoint-threshold",
            "threshold for saving checkpoints");
    param_list_decl_usage(pl, "lingen_checkpoint_save_gathered",
            "save global checkpoints files, instead of per-job files");
}

template<bool is_binary>
void lingen_checkpoint<is_binary>::lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "checkpoint-directory");
    param_list_lookup_string(pl, "checkpoint-threshold");
    param_list_lookup_string(pl, "lingen_checkpoint_save_gathered");
}

template<bool is_binary>
void lingen_checkpoint<is_binary>::interpret_parameters(cxx_param_list & pl)
{
    param_list_parse(pl, "checkpoint-directory", default_directory);
    param_list_parse(pl, "checkpoint-threshold", threshold);
    param_list_parse(pl, "lingen_checkpoint_save_gathered", save_gathered);


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!save_gathered && !default_directory.empty()) {
        int ok = access(default_directory.c_str(), X_OK) == 0;
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (!ok) {
            if (!rank)
                fmt::print("# Checkpoint directory {}/ does not exist at all ranks, falling back to gathered checkpoints\n", default_directory);
            lingen_checkpoint<is_binary>::save_gathered = 1;
        }
    }
    if (save_gathered) {
        if (!default_directory.empty()) {
            int ok = access(default_directory.c_str(), X_OK) == 0;
            MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (!ok) {
                if (!rank)
                    fmt::print("# Checkpoint directory {}/ does not exist"
                            " at rank zero, checkpoint settings ignored\n",
                            default_directory);
                default_directory.clear();
                threshold = UINT_MAX;
            }
        }
    }
}

static bool remove_suffix(std::string & S, std::string const & suffix)
{
    if (has_suffix(S.c_str(), suffix.c_str())) {
        S = S.substr(0, S.size() - suffix.size());
        return true;
    } else {
        return false;
    }
}

static void remove_suffix_regexp(std::string & S, std::string const & suffix)
{
    S = regex_replace(S, std::regex(suffix + "$"), "");
}


template<bool is_binary>
lingen_checkpoint<is_binary>::lingen_checkpoint(bmstatus<is_binary> & bm, unsigned int t0, unsigned int t1, int mpi, std::string base)
    : bm(bm)
    , level(bm.depth)
    , t0(t0)
    , t1(t1)
    , target_t(t0)
    , mpi(mpi)
{
    if (mpi)
        MPI_Comm_rank(bm.com[0], &(rank));
    else
        rank = 0;

    /* If the filename ends with any of the suffixes we generally attach
     * to checkpoint file names, then derive a new basename from that.
     */
    remove_suffix(base, ".aux");
    remove_suffix(base, ".single.data");
    remove_suffix_regexp(base, R"(\.\d+\.data)");

    auxfile = base + ".aux";
    gdatafile = base + ".single.data";
    sdatafile = base + fmt::format(".{}.data", rank);

    datafile = mpi ? sdatafile : gdatafile;
}

template<bool is_binary>
std::string lingen_checkpoint<is_binary>::get_cp_basename(bmstatus<is_binary> & bm, cp_which which, unsigned int t0, unsigned int t1)
{
    int const level = bm.depth;
    std::string base;
    switch(which) {
        case LINGEN_CHECKPOINT_E:
            base = fmt::format("{}.{}.{}.E", level, t0, t1);
            break;
        case LINGEN_CHECKPOINT_PI:
            base = fmt::format("{}.{}.{}.pi", level, t0, t1);
            break;
    }
    if (!default_directory.empty())
        base = fmt::format("{}/{}", default_directory, base);
    return base;
}


template<bool is_binary>
bool lingen_checkpoint<is_binary>::save_aux_file(size_t Xsize) const /*{{{*/
{
    auto const& d = bm.d;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    if (rank) return true;
    std::ofstream os(auxfile);
    os << "format " << format << "\n";
    os << m << " " << n << "\n";
    os << level << " " << t0 << " " << t1 << " " << bm.t << "\n";
    os << Xsize << "\n";
    os << bm.d.ab.characteristic() << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.delta[i];
    os << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.lucky[i];
    os << "\n";
    os << bm.done;
    os << "\n";
    os << logline_serialize();
    os << "\n";
    try {
        os << bm.hints;
    } catch (std::runtime_error const & e) {
        std::string s = e.what();
        s += " when writing checkpoint";
        throw std::runtime_error(s);
    }
    os << "\n";
    os << bm.stats;
    os << "\n";
    bool const ok = os.good();
    if (!ok)
        unlink(auxfile.c_str());
    return ok;
}/*}}}*/

/* See load_mpi_checkpoint_file_* for the way we reconcile the aux file
 * that is read at rank 0 and the in-memory structure at all ranks.
 */

template<bool is_binary>
bool lingen_checkpoint<is_binary>::checkpoint_already_present() const/*{{{*/
{
    lingen_checkpoint test = *this;
    size_t Xsize;
    int ok=0;
    int exc=0;

    /* Do we have an aux file at rank 0 ? */
    try {
        ok = test.load_aux_file(Xsize);     /* rank>0 always returns 1 */
    } catch (typename lingen_checkpoint<is_binary>::invalid_aux_file const & inv) {
        if (!rank)
            fmt::print(stderr, "Overwriting bogus checkpoint file {} [{}]\n",
                    datafile, inv.what());
        exc = 1;
    }
    if (mpi) MPI_Allreduce(MPI_IN_PLACE, &exc, 1, MPI_INT, MPI_MAX, bm.com[0]);
    if (exc) return false;

    if (mpi) MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (!ok) return false;

    /* Do we have a scattered file set (preferred) ? */
    ok = access(sdatafile.c_str(), R_OK) == 0;
    if (mpi) MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (ok) return true;

    /* Do we have a gathered file ? */
    ok = rank || access(gdatafile.c_str(), R_OK) == 0;
    if (mpi) MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (ok) return true;

    return false;
}/*}}}*/

/* This function can be used to parse a checkpoint header _and_ use the
 * context directly from the header.
 */
template<bool is_binary>
std::istream& lingen_checkpoint<is_binary>::header_info::input(std::istream & is)
{
    std::string hfstring;
    unsigned long hformat;
    is >> hfstring >> hformat;
    if (hfstring != "format")
        throw invalid_aux_file("version < 1");
    if (hformat != lingen_checkpoint<is_binary>::format)
        throw invalid_aux_file(
                fmt::format("version {} < {}",
                    hformat, lingen_checkpoint<is_binary>::format));
    if (!(is >> m >> n))
        throw invalid_aux_file("parse error");
    if (!(is >> level >> t0 >> t1 >> t))
        throw invalid_aux_file("parse error");
    if (t < t0 || t > t1)
        throw invalid_aux_file("t is out of bounds");
    if (!(is >> ncoeffs >> p))
        throw invalid_aux_file("parse error");

    is >> read_container(delta, m + n);
    is >> read_container(lucky, m + n);
    is >> done;
    return is;
}

template<bool is_binary>
bool lingen_checkpoint<is_binary>::load_aux_file(size_t & Xsize)/*{{{*/
{
    if (rank) return true;
    std::ifstream is(auxfile);
    if (!is.good()) return false;

    bmstatus<is_binary> nbm = bm;

    header_info u;
    try {
        if (!(is >> u))
            return false;
        if (u.m != bm.d.m || u.n != bm.d.n)
            throw invalid_aux_file(fmt::format("made for (m,n)=({},{})",
                        u.m, u.n));
        if (level != u.level || t0 != u.t0 || t1 != u.t1)
            throw invalid_aux_file(fmt::format("made for depth={} t0={} t1={}",
                        u.level, u.t0, u.t1));
        if (mpz_cmp(u.p, bm.d.ab.characteristic()) != 0)
            throw invalid_aux_file("made for wrong p");

        target_t = u.t;
        nbm.delta = u.delta;
        nbm.lucky = u.lucky;
        nbm.done = u.done;
        Xsize = runtime_numeric_cast<size_t>(u.ncoeffs);

        if (u.t0 > u.t || u.t > u.t1)
            throw invalid_aux_file("inconsistent values for t");
        if (!(u.t0 + u.ncoeffs <= u.t1))
            throw invalid_aux_file("inconsistent values for t and/or ncoeffs");
    } catch (invalid_aux_file const& e) {
        throw invalid_aux_file(
                fmt::format(
                    "checkpoint file {} cannot be used ({}",
                    auxfile, e.what()));
    }

    double tt;
    is >> tt;
    logline_unserialize(tt);
    if (!(is >> nbm.hints)) {
        throw invalid_aux_file("parse error while reading hints file from checkpoint");
    }
    for(auto const & x : nbm.hints) {
        if (!x.second.check()) {
            is.setstate(std::ios::failbit);
            throw invalid_aux_file("checkpoint contains invalid schedule information");
        }
    }

    if (bm.hints != nbm.hints) {
        is.setstate(std::ios::failbit);
        std::stringstream os;
        os << bm.hints;
        throw invalid_aux_file(
                    "checkpoint file cannot be used since"
                    " it was made for another set of schedules (stats"
                    " would be incorrect)\n"
                    "textual description of the schedule set that we"
                    " expect to find:\n"
                + os.str());
    }

    if (!(is >> nbm.stats))
        throw invalid_aux_file("stats in checkpoint file cannot be merged with existing stats");
   
    bm = std::move(nbm);

    return is.good();
}/*}}}*/

/* TODO: adapt for GF(2) */
template<bool is_binary>
int lingen_checkpoint<is_binary>::load_data_file(matpoly<is_binary> & X)/*{{{*/
{
    auto & d = bm.d;
    auto & ab = d.ab;
    FILE * data = fopen(datafile.c_str(), "rb");
    int rc;
    if (data == nullptr) {
        fmt::print(stderr, "Warning: cannot open {}\n", datafile);
        return 0;
    }
    rc = lingen_io_matpoly<is_binary>::read(&ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
/* I think we always have Xsize == X.size, the only questionable
 * situation is when we're saving part of a big matrix */
template<bool is_binary>
int lingen_checkpoint<is_binary>::save_data_file(matpoly<is_binary> const & X)/*{{{*/
{
    int ok = 0;
    auto unlink_aux = call_dtor([&]() {if (!ok) unlink(auxfile.c_str());});
    auto & ab = bm.d.ab;
    std::ofstream data(datafile, std::ios_base::out | std::ios_base::binary);
    int rc;
    if (!data) {
        fmt::print(stderr, "Warning: cannot open {}\n", datafile);
        return 0;
    }
    auto unlink_data = call_dtor([&]() {if (!ok) unlink(datafile.c_str());});
    rc = lingen_io_matpoly<is_binary>::write(&ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size())
        return 0;
    ok = data.good();;
    return ok;
}/*}}}*/

template<bool is_binary>
int load_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, matpoly<is_binary> & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    auto & d = bm.d;
    auto * ab = & d.ab;
    unsigned int const m = d.m;
    unsigned int const n = d.n;

    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;

    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 0);

    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = 1;
    try {
        ok = cp.load_aux_file(Xsize);
        if (!ok) return false;
    } catch (typename lingen_checkpoint<is_binary>::invalid_aux_file const & inv) {
        fmt::print(stderr, "Error with checkpoint file {}:\n{}\n", cp.datafile, inv.what());
        abort();
    }
    logline_begin(stdout, SIZE_MAX, "Reading %s", cp.datafile.c_str());
    switch (which) {
        case LINGEN_CHECKPOINT_E:
            X = matpoly<is_binary>(ab, m, m+n, runtime_numeric_cast<int>(Xsize));
            break;
        case LINGEN_CHECKPOINT_PI:
            X = matpoly<is_binary>(ab, m+n, m+n, runtime_numeric_cast<int>(Xsize));
            break;
    }
    X.set_size(Xsize);
    ok = cp.load_data_file(X);
    X.clear_high_word();
    logline_end(&bm.t_cp_io,"");
    if (!ok)
        fmt::print(stderr, "Warning: I/O error while reading {}\n", cp.datafile);
    bm.t = cp.target_t;
    return ok;
}/*}}}*/

template<bool is_binary>
int save_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, matpoly<is_binary> const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    // coverity[tainted_data_transitive]
    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 0);
    if (cp.checkpoint_already_present())
        return 1;
    logline_begin(stdout, SIZE_MAX, "Saving %s%s",
            cp.datafile.c_str(),
            cp.mpi ? " (MPI, scattered)" : "");
    int ok = cp.save_aux_file(X.get_size());
    if (ok) ok = cp.save_data_file(X);
    logline_end(&bm.t_cp_io,"");
    if (!ok && !cp.rank)
        fmt::print(stderr, "Warning: I/O error while saving {}\n", cp.datafile);
    return ok;
}/*}}}*/

/* At some point the aux file, which is read at rank 0 only, must bcast
 * its data to other ranks. We list here the fields of interest.
 *
 * Maybe a function that is dedicated to sharing this structure info
 * would be a good idea.
    target_t
    Xsize
    bm.delta[0]@(m+n)
    bm.lucky[0]@(m+n)
    bm.done
    bm.hints -- should match the expected stuff. no need to share !
    bm.stats -- not needed at rank>0
 */

template<bool is_binary>
static int load_mpi_checkpoint_file_scattered(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    auto const& d = bm.d;
    auto * ab = & bm.d.ab;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 1);
    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    int ok = 1, error = 0;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (typename lingen_checkpoint<is_binary>::invalid_aux_file const & inv) {
        if (!rank)
            fmt::print(stderr, "Error with checkpoint file {}:\n{}\n", cp.datafile, inv.what());
        error = 1;
    }
    MPI_Bcast(&error, 1, MPI_INT, 0, bm.com[0]);
    ASSERT_ALWAYS(!error);
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (!ok) {
        return false;
    }
    MPI_Bcast(&Xsize, 1, CADO_MPI_SIZE_T, 0, bm.com[0]);
    int const b = runtime_numeric_cast<int>(m + n);
    MPI_Bcast(bm.delta.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(bm.lucky.data(), b, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    int commsize;
    MPI_Comm_size(bm.com[0], &commsize);
    printf("rk=%d, Xsize=%zu\n", rank, Xsize);
    logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, scattered)",
            cp.datafile.c_str());
    switch (which) {
        case LINGEN_CHECKPOINT_E:
            X.finish_init(ab, m, m+n, runtime_numeric_cast<int>(Xsize));
            break;
        case LINGEN_CHECKPOINT_PI:
            X.finish_init(ab, m+n, m+n, runtime_numeric_cast<int>(Xsize));
            break;
    }
    X.set_size(Xsize);
    ok = cp.load_data_file(X.my_cell());
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fmt::print(stderr, "\nWarning: I/O error while reading {}\n",
                cp.datafile, Xsize);
    }
    bm.t = cp.target_t;
    MPI_Bcast(&bm.t, 1, MPI_UNSIGNED, 0, bm.com[0]);
    return ok;
}/*}}}*/

template<bool is_binary>
static int save_mpi_checkpoint_file_scattered(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 1);
    int ok;
    ok = cp.checkpoint_already_present();
    if (ok) return 1;

    ok = cp.save_aux_file(X.get_size());
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (!ok && !rank) unlink(cp.auxfile.c_str());
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Saving %s (MPI, scattered)",
                cp.datafile.c_str());
        ok = cp.save_data_file(X.my_cell());
        logline_end(&bm.t_cp_io,"");
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok) {
            if (!cp.datafile.empty()) unlink(cp.datafile.c_str());
            if (!rank) unlink(cp.auxfile.c_str());
        }
    }
    if (!ok && !rank) {
        fmt::print(stderr, "Warning: I/O error while saving {}\n", cp.datafile);
    }
    return ok;
}/*}}}*/

template<bool is_binary>
static int load_mpi_checkpoint_file_gathered(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    auto & d = bm.d;
    auto * ab = & d.ab;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    size_t Xsize = 0;
    int ok = 1, error = 0;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (typename lingen_checkpoint<is_binary>::invalid_aux_file const & inv) {
        error = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_INT, MPI_MAX, bm.com[0]);
    ASSERT_ALWAYS(!error);
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (!ok)
        /* This can be a normal situation if the .aux file does not
         * exist.
         */
        return false;
    MPI_Bcast(&Xsize, 1, CADO_MPI_SIZE_T, 0, bm.com[0]);
    int const b = runtime_numeric_cast<int>(m + n);
    MPI_Bcast(bm.delta.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(bm.lucky.data(), b, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, gathered)",
            cp.datafile.c_str());
    do {
        std::unique_ptr<FILE, delete_FILE> data;
        if (!rank) {
            data.reset(fopen(cp.datafile.c_str(), "rb"));
            ok = bool(data);
        }
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok) {
            if (!rank)
                fmt::print(stderr, "Warning: cannot open {}\n", cp.datafile);
            break;
        }

        switch (which) {
            case LINGEN_CHECKPOINT_E:
                X.finish_init(ab, m, m+n, runtime_numeric_cast<int>(Xsize));
                break;
            case LINGEN_CHECKPOINT_PI:
                X.finish_init(ab, m+n, m+n, runtime_numeric_cast<int>(Xsize));
                break;
        }
        X.set_size(Xsize);

        double const avg = ::average_matsize<is_binary>(ab, X.m, X.n, 0);
        size_t const B = iceildiv(lingen_io_matpoly<is_binary>::block_size, avg);

        /* This is only temp storage ! */
        matpoly<is_binary> loc(ab, X.m, X.n, runtime_numeric_cast<int>(B));
        loc.zero_pad(B);

        for(size_t k = 0 ; ok && k < X.get_size() ; k += B) {
            size_t const nc = std::min(B, X.get_size() - k);
            if (!rank) {
                auto rc = lingen_io_matpoly<is_binary>::read(ab, data.get(), loc, 0, nc, 0, 0);
                ok = rc == runtime_numeric_cast<int>(nc);
            }
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            X.scatter_mat_partial(loc, 0, k, nc);
        }
    } while (false);
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    logline_end(&bm.t_cp_io, "");
    if (ok)
        bm.t = cp.target_t;
    else if (!rank)
            fmt::print(stderr, "Warning: I/O error while reading {}\n",
                    cp.datafile);
    MPI_Bcast(&bm.t, 1, MPI_UNSIGNED, 0, bm.com[0]);

    return ok;
}/*}}}*/

template<bool is_binary>
static int save_mpi_checkpoint_file_gathered(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    auto & d = bm.d;
    auto * ab = & d.ab;
    lingen_checkpoint<is_binary> cp(bm, which, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    int ok = cp.checkpoint_already_present();
    if (ok) return 1;
    logline_begin(stdout, SIZE_MAX, "Saving %s (MPI, gathered)",
            cp.datafile.c_str());
    ok = cp.save_aux_file(X.get_size());
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        do {
            std::ofstream data;
            if (!rank) {
                data.open(cp.datafile, std::ios_base::out | std::ios_base::binary);
                ok = (bool) data;
            }
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fmt::print(stderr, "Warning: cannot open {}\n", cp.datafile);
                break;
            }

            double const avg = ::average_matsize<is_binary>(ab, X.m, X.n, 0);
            size_t const B = iceildiv(lingen_io_matpoly<is_binary>::block_size, avg);

            /* This is only temp storage ! */
            matpoly<is_binary> loc(ab, X.m, X.n, runtime_numeric_cast<int>(B));
            loc.zero_pad(B);

            for(size_t k = 0 ; ok && k < X.get_size() ; k += B) {
                size_t const nc = std::min(B, X.get_size() - k);
                X.gather_mat_partial(loc, 0, k, nc);
                if (!rank) {
                    auto rc = lingen_io_matpoly<is_binary>::write(ab, data, loc, 0, nc, 0, 0);
                    ok = rc == runtime_numeric_cast<int>(nc);
                }
                MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            }

            if (!rank) {
                data.close();
                ok = ok && (bool) data;
            }
        } while (false);
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok && !rank) {
            if (!cp.datafile.empty()) unlink(cp.datafile.c_str());
            unlink(cp.auxfile.c_str());
        }
    }
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fmt::print(stderr, "Warning: I/O error while saving {}\n", cp.datafile);
    }
    return ok;
}/*}}}*/

template<bool is_binary>
int load_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* read scattered checkpoint with higher priority if available,
     * because we like distributed I/O. Otherwise, read gathered
     * checkpoint if we could find one.
     */
    if ((t1 - t0) < lingen_checkpoint<is_binary>::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    lingen_checkpoint<is_binary> const cp(bm, which, t0, t1, 1);
    int ok = 0;
    int const aux_ok = rank || access(cp.auxfile.c_str(), R_OK) == 0;
    int const sdata_ok = access(cp.sdatafile.c_str(), R_OK) == 0;
    int scattered_ok = aux_ok && sdata_ok;
    MPI_Allreduce(MPI_IN_PLACE, &scattered_ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (scattered_ok) {
        ok = load_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
        if (ok) return ok;
    }
    int const gdata_ok = rank || access(cp.gdatafile.c_str(), R_OK) == 0;
    int gathered_ok = aux_ok && gdata_ok;
    MPI_Bcast(&gathered_ok, 1, MPI_INT, 0, bm.com[0]);
    if (gathered_ok) {
        ok = load_mpi_checkpoint_file_gathered(bm, which, X, t0, t1);
    }
    return ok;
}/*}}}*/

template<bool is_binary>
int save_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if (lingen_checkpoint<is_binary>::save_gathered) {
        return save_mpi_checkpoint_file_gathered(bm, which, X, t0, t1);
    } else {
        return save_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
    }
}/*}}}*/

#ifdef LINGEN_BINARY
template struct lingen_checkpoint<true>;
// template std::string lingen_checkpoint<true>::default_directory;
template
int load_checkpoint_file<true>(bmstatus<true> & bm, cp_which which, matpoly<true> &, unsigned int t0, unsigned int t1);
template
int save_checkpoint_file<true>(bmstatus<true> & bm, cp_which which, matpoly<true> const &, unsigned int t0, unsigned int t1);
template
int load_checkpoint_file<true>(bmstatus<true> & bm, cp_which which, bigmatpoly<true> &, unsigned int t0, unsigned int t1);
template
int save_checkpoint_file<true>(bmstatus<true> & bm, cp_which which, bigmatpoly<true> const &, unsigned int t0, unsigned int t1);

#else
template struct lingen_checkpoint<false>;
// template std::string lingen_checkpoint<false>::default_directory;
template
int load_checkpoint_file<false>(bmstatus<false> & bm, cp_which which, matpoly<false> &, unsigned int t0, unsigned int t1);
template
int save_checkpoint_file<false>(bmstatus<false> & bm, cp_which which, matpoly<false> const &, unsigned int t0, unsigned int t1);
template
int load_checkpoint_file<false>(bmstatus<false> & bm, cp_which which, bigmatpoly<false> &, unsigned int t0, unsigned int t1);
template
int save_checkpoint_file<false>(bmstatus<false> & bm, cp_which which, bigmatpoly<false> const &, unsigned int t0, unsigned int t1);

#endif
