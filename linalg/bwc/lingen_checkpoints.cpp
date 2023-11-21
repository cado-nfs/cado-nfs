#include "cado.h"
// IWYU pragma: no_include "mpfq_fake.hpp"
// IWYU pragma: no_include "mpfq_layer.h"
// IWYU pragma: no_include <sys/param.h>
#include <cstdio>
#include <istream> // IWYU pragma: keep
#include <fstream> // IWYU pragma: keep
#include <gmp.h>                       // for operator<<, mpz_cmp
#include <stdint.h>                    // for SIZE_MAX
#include <stdlib.h>                    // for abort, free
#include <unistd.h>                    // for unlink, access, R_OK, X_OK
#include <map>                         // for operator!=, map
#include <utility>                     // for move, pair
#include <vector>                      // for vector
#include "fmt/core.h"                  // for check_format_string, char_t
#include "fmt/format.h"                // for basic_buffer::append, basic_pa...
#include "arith-hard.hpp"
#include "lingen_bw_dimensions.hpp"
#include "lingen_call_companion.hpp"   // for lingen_call_companion, lingen_...
#include "lingen_hints.hpp"            // for operator<<, lingen_hints, oper...
#include "lingen_matpoly_select.hpp"   // for matpoly
#include "macros.h"                    // for ASSERT_ALWAYS, MIN, iceildiv
#include "params.h"                    // for cxx_param_list, param_list_loo...
#include "select_mpi.h"                // for MPI_Allreduce, MPI_Bcast, MPI_...
#include "tree_stats.hpp"              // for operator<<, operator>>
#include "lingen_bmstatus.hpp"
#include "lingen_checkpoints.hpp"
#include "lingen_io_matpoly.hpp"
#include "lingen_average_matsize.hpp"
#include "logline.h"
#include "fmt/printf.h" // IWYU pragma: keep
#include "cxx_mpz.hpp"

/* Checkpoints */

const char * lingen_checkpoint::directory = NULL;
unsigned int lingen_checkpoint::threshold = 100;
int lingen_checkpoint::save_gathered = 0;
constexpr unsigned long lingen_checkpoint::format;

/* There's much copy-paste here */


void lingen_checkpoint::decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "checkpoint-directory",
            "where to save checkpoints");
    param_list_decl_usage(pl, "checkpoint-threshold",
            "threshold for saving checkpoints");
    param_list_decl_usage(pl, "lingen_checkpoint_save_gathered",
            "save global checkpoints files, instead of per-job files");
}

void lingen_checkpoint::lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "checkpoint-directory");
    param_list_lookup_string(pl, "checkpoint-threshold");
    param_list_lookup_string(pl, "lingen_checkpoint_save_gathered");
}

void lingen_checkpoint::interpret_parameters(cxx_param_list & pl)
{
    lingen_checkpoint::directory = param_list_lookup_string(pl, "checkpoint-directory");
    param_list_parse_uint(pl, "checkpoint-threshold", &lingen_checkpoint::threshold);
    param_list_parse_int(pl, "lingen_checkpoint_save_gathered", &lingen_checkpoint::save_gathered);

    if (!lingen_checkpoint::directory) return;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!lingen_checkpoint::save_gathered) {
        int ok = access(lingen_checkpoint::directory, X_OK) == 0;
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (!ok) {
            if (!rank)
                printf("# Checkpoint directory %s/ does not exist at all ranks, falling back to gathered checkpoints\n", lingen_checkpoint::directory);
            lingen_checkpoint::save_gathered = 1;
        }
    }
    if (lingen_checkpoint::save_gathered) {
        int ok = access(lingen_checkpoint::directory, X_OK) == 0;
        MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!ok) {
            if (!rank)
                printf("# Checkpoint directory %s/ does not exist at rank zero, checkpoint settings ignored\n", lingen_checkpoint::directory);
            lingen_checkpoint::directory = NULL;
        }
    }
}

lingen_checkpoint::lingen_checkpoint(bmstatus & bm, unsigned int t0, unsigned int t1, int mpi, std::string base)
    : bm(bm), t0(t0), t1(t1), mpi(mpi)
{
    level = bm.depth();
    if (mpi)
        MPI_Comm_rank(bm.com[0], &(rank));
    else
        rank = 0;
    auxfile = base + ".aux";
    gdatafile = base + ".single.data";
    sdatafile = base + fmt::sprintf(".%d.data", rank);
    if (directory) {
        auxfile = fmt::sprintf("%s/%s", directory, auxfile);
        gdatafile = fmt::sprintf("%s/%s", directory, gdatafile);
        sdatafile = fmt::sprintf("%s/%s", directory, sdatafile);
    }

    datafile = mpi ? sdatafile : gdatafile;
}

std::string lingen_checkpoint::get_cp_basename(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1)
{
    int level = bm.depth();
    std::string base;
    switch(which) {
        case LINGEN_CHECKPOINT_E:
            base = fmt::sprintf("E.%d.%u", level, t0);
            break;
        case LINGEN_CHECKPOINT_PI:
            base = fmt::sprintf("pi.%d.%u.%u", level, t0, t1);
            break;
    }
    return base;
}


bool lingen_checkpoint::save_aux_file(size_t Xsize) const /*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    if (rank) return 1;
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
    bool ok = os.good();
    if (!ok)
        unlink(auxfile.c_str());
    return ok;
}/*}}}*/

/* See load_mpi_checkpoint_file_* for the way we reconcile the aux file
 * that is read at rank 0 and the in-memory structure at all ranks.
 */

bool lingen_checkpoint::checkpoint_already_present() const/*{{{*/
{
    lingen_checkpoint test = *this;
    size_t Xsize;
    int ok=0;
    int exc=0;

    /* Do we have an aux file at rank 0 ? */
    try {
        ok = test.load_aux_file(Xsize);     /* rank>0 always returns 1 */
    } catch (lingen_checkpoint::invalid_aux_file const & inv) {
        if (!rank)
        fmt::fprintf(stderr, "Overwriting bogus checkpoint file %s [%s]\n",
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

bool lingen_checkpoint::load_aux_file(size_t & Xsize)/*{{{*/
{
    bmstatus nbm = bm;
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    if (rank) return 1;
    std::ifstream is(auxfile);
    if (!is.good()) return false;
    std::string hfstring;
    unsigned long hformat;
    is >> hfstring >> hformat;
    if (hfstring != "format")
        throw invalid_aux_file("checkpoint file cannot be used (version < 1)");
    if (hformat != format)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (version %lu < %lu)", hformat, format));
    unsigned int xm,xn;
    is >> xm >> xn;
    if (xm != m || xn != n)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (made for (m,n)=(%u,%u)", xm, xn));
    int xlevel;
    unsigned int xt0, xt1;
    if (!(is >> xlevel >> xt0 >> xt1 >> target_t))
        throw invalid_aux_file("checkpoint file cannot be used (parse error)");
    if (xlevel != level || t0 != xt0 || t1 != xt1)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (made for depth=%d t0=%u t1=%u", xlevel, xt0, xt1));
    ASSERT_ALWAYS(target_t <= t1);
    is >> Xsize;
    cxx_mpz xp;
    is >> xp;
    if (mpz_cmp(xp, bm.d.ab.characteristic()) != 0)
        throw invalid_aux_file("checkpoint file cannot be used (made for wrong p)");
    for(unsigned int i = 0 ; i < m + n ; i++) {
        is >> nbm.delta[i];
    }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        is >> nbm.lucky[i];
    }
    is >> nbm.done;
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
                fmt::sprintf(
                    "checkpoint file cannot be used since"
                    " it was made for another set of schedules (stats"
                    " would be incorrect)\n"
                ) + "textual description of the schedule set that we"
                    " expect to find:\n"
                + os.str());
    }

    if (!(is >> nbm.stats))
        throw invalid_aux_file("stats in checkpoint file cannot be merged with existing stats");
   
    bm = std::move(nbm);

    return is.good();
}/*}}}*/

/* TODO: adapt for GF(2) */
int lingen_checkpoint::load_data_file(matpoly & X)/*{{{*/
{
    bw_dimensions & d = bm.d;
    matpoly::arith_hard & ab = d.ab;
    FILE * data = fopen(datafile.c_str(), "rb");
    int rc;
    if (data == NULL) {
        fmt::fprintf(stderr, "Warning: cannot open %s\n", datafile);
        return 0;
    }
    rc = matpoly_read(&ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
/* I think we always have Xsize == X.size, the only questionable
 * situation is when we're saving part of a big matrix */
int lingen_checkpoint::save_data_file(matpoly const & X)/*{{{*/
{
    matpoly::arith_hard & ab = bm.d.ab;
    std::ofstream data(datafile, std::ios_base::out | std::ios_base::binary);
    int rc;
    if (!data) {
        fmt::fprintf(stderr, "Warning: cannot open %s\n", datafile);
        unlink(auxfile.c_str());
        return 0;
    }
    rc = matpoly_write(&ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) goto lingen_checkpoint_save_data_file_bailout;
    if (data.good()) return 1;
lingen_checkpoint_save_data_file_bailout:
    unlink(datafile.c_str());
    unlink(auxfile.c_str());
    return 0;
}/*}}}*/

template<>
int load_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    bw_dimensions & d = bm.d;
    matpoly::arith_hard * ab = & d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;

    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;

    lingen_checkpoint cp(bm, which, t0, t1, 0);

    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = 1;
    try {
        ok = cp.load_aux_file(Xsize);
        if (!ok) return false;
    } catch (lingen_checkpoint::invalid_aux_file const & inv) {
        fmt::fprintf(stderr, "Error with checkpoint file %s:\n%s\n", cp.datafile, inv.what());
        abort();
    }
    logline_begin(stdout, SIZE_MAX, "Reading %s", cp.datafile.c_str());
    switch (which) {
        case LINGEN_CHECKPOINT_E:
            X = matpoly(ab, m, m+n, Xsize);
            break;
        case LINGEN_CHECKPOINT_PI:
            X = matpoly(ab, m+n, m+n, Xsize);
            break;
    }
    X.set_size(Xsize);
    ok = cp.load_data_file(X);
    X.clear_high_word();
    logline_end(&bm.t_cp_io,"");
    if (!ok)
        fmt::fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
    bm.t = cp.target_t;
    return ok;
}/*}}}*/

template<>
int save_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    // coverity[tainted_data_transitive]
    lingen_checkpoint cp(bm, which, t0, t1, 0);
    if (cp.checkpoint_already_present())
        return 1;
    logline_begin(stdout, SIZE_MAX, "Saving %s%s",
            cp.datafile.c_str(),
            cp.mpi ? " (MPI, scattered)" : "");
    int ok = cp.save_aux_file(X.get_size());
    if (ok) ok = cp.save_data_file(X);
    logline_end(&bm.t_cp_io,"");
    if (!ok && !cp.rank)
        fmt::fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
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

int load_mpi_checkpoint_file_scattered(bmstatus & bm, cp_which which, bigmatpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    matpoly::arith_hard * ab = & bm.d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    lingen_checkpoint cp(bm, which, t0, t1, 1);
    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    int ok = 1, error = 0;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (lingen_checkpoint::invalid_aux_file const & inv) {
        if (!rank)
            fmt::fprintf(stderr, "Error with checkpoint file %s:\n%s\n", cp.datafile, inv.what());
        error = 1;
    }
    MPI_Bcast(&error, 1, MPI_INT, 0, bm.com[0]);
    ASSERT_ALWAYS(!error);
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (!ok) {
        return false;
    }
    MPI_Bcast(&Xsize, 1, CADO_MPI_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(bm.delta.data(), m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(bm.lucky.data(), m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    int commsize;
    MPI_Comm_size(bm.com[0], &commsize);
    logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, scattered)",
            cp.datafile.c_str());
    switch (which) {
        case LINGEN_CHECKPOINT_E:
            X.finish_init(ab, m, m+n, Xsize);
            break;
        case LINGEN_CHECKPOINT_PI:
            X.finish_init(ab, m+n, m+n, Xsize);
            break;
    }
    X.set_size(Xsize);
    ok = cp.load_data_file(X.my_cell());
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fmt::fprintf(stderr, "\nWarning: I/O error while reading %s\n",
                cp.datafile, Xsize);
    }
    bm.t = cp.target_t;
    MPI_Bcast(&bm.t, 1, MPI_INT, 0, bm.com[0]);
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_scattered(bmstatus & bm, cp_which which, bigmatpoly const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    lingen_checkpoint cp(bm, which, t0, t1, 1);
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
        fmt::fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
    }
    return ok;
}/*}}}*/

int load_mpi_checkpoint_file_gathered(bmstatus & bm, cp_which which, bigmatpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    matpoly::arith_hard * ab = & d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    lingen_checkpoint cp(bm, which, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    size_t Xsize = 0;
    int ok = 1, error = 0;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (lingen_checkpoint::invalid_aux_file const & inv) {
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
    MPI_Bcast(bm.delta.data(), m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(bm.lucky.data(), m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, gathered)",
            cp.datafile.c_str());
    do {
        FILE * data = NULL;
        if (!rank) ok = (data = fopen(cp.datafile.c_str(), "rb")) != NULL;
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok) {
            if (!rank)
                fprintf(stderr, "Warning: cannot open %s\n", cp.datafile.c_str());
            if (data) free(data);
            break;
        }

        switch (which) {
            case LINGEN_CHECKPOINT_E:
                X.finish_init(ab, m, m+n, Xsize);
                break;
            case LINGEN_CHECKPOINT_PI:
                X.finish_init(ab, m+n, m+n, Xsize);
                break;
        }
        X.set_size(Xsize);

        double avg = average_matsize(ab, X.m, X.n, 0);
        unsigned int B = iceildiv(io_matpoly_block_size, avg);

        /* This is only temp storage ! */
        matpoly loc(ab, X.m, X.n, B);
        loc.zero_pad(B);

        for(unsigned int k = 0 ; ok && k < X.get_size() ; k += B) {
            unsigned int nc = MIN(B, X.get_size() - k);
            if (!rank)
                ok = matpoly_read(ab, data, loc, 0, nc, 0, 0) == (int) nc;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            X.scatter_mat_partial(loc, 0, k, nc);
        }

        if (!rank) {
            int rc = fclose(data);
            ok = ok && (rc == 0);
        }
    } while (0);
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    logline_end(&bm.t_cp_io,"");
    if (ok)
        bm.t = cp.target_t;
    else if (!rank)
            fmt::fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp.datafile);
    MPI_Bcast(&bm.t, 1, MPI_INT, 0, bm.com[0]);

    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_gathered(bmstatus & bm, cp_which which, bigmatpoly const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    matpoly::arith_hard * ab = & d.ab;
    lingen_checkpoint cp(bm, which, t0, t1, 1);
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
                    fmt::fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                break;
            }

            double avg = average_matsize(ab, X.m, X.n, 0);
            unsigned int B = iceildiv(io_matpoly_block_size, avg);

            /* This is only temp storage ! */
            matpoly loc(ab, X.m, X.n, B);
            loc.zero_pad(B);

            for(unsigned int k = 0 ; ok && k < X.get_size() ; k += B) {
                unsigned int nc = MIN(B, X.get_size() - k);
                X.gather_mat_partial(loc, 0, k, nc);
                if (!rank)
                    ok = matpoly_write(ab, data, loc, 0, nc, 0, 0) == (int) nc;
                MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            }

            if (!rank) {
                data.close();
                ok = ok && (bool) data;
            }
        } while (0);
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok && !rank) {
            if (!cp.datafile.empty()) unlink(cp.datafile.c_str());
            unlink(cp.auxfile.c_str());
        }
    }
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fmt::fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
    }
    return ok;
}/*}}}*/

template<>
int load_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* read scattered checkpoint with higher priority if available,
     * because we like distributed I/O. Otherwise, read gathered
     * checkpoint if we could find one.
     */
    if (!lingen_checkpoint::directory) return 0;
    if ((t1 - t0) < lingen_checkpoint::threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    lingen_checkpoint cp(bm, which, t0, t1, 1);
    int ok = 0;
    int aux_ok = rank || access(cp.auxfile.c_str(), R_OK) == 0;
    int sdata_ok = access(cp.sdatafile.c_str(), R_OK) == 0;
    int scattered_ok = aux_ok && sdata_ok;
    MPI_Allreduce(MPI_IN_PLACE, &scattered_ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (scattered_ok) {
        ok = load_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
        if (ok) return ok;
    }
    int gdata_ok = rank || access(cp.gdatafile.c_str(), R_OK) == 0;
    int gathered_ok = aux_ok && gdata_ok;
    MPI_Bcast(&gathered_ok, 1, MPI_INT, 0, bm.com[0]);
    if (gathered_ok) {
        ok = load_mpi_checkpoint_file_gathered(bm, which, X, t0, t1);
    }
    return ok;
}/*}}}*/

template<>
int save_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if (lingen_checkpoint::save_gathered) {
        return save_mpi_checkpoint_file_gathered(bm, which, X, t0, t1);
    } else {
        return save_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
    }
}/*}}}*/
