#include "cado.h"
#include <cstdio>
#include <ostream>
#include <istream>
#include <fstream>
#include <sstream>
#include "lingen_bmstatus.hpp"
#include "lingen_checkpoints.hpp"
#include "lingen_io_matpoly.hpp"
#include "lingen_average_matsize.hpp"
#include "logline.h"
#include "fmt/printf.h"

/* Checkpoints */

static const char * checkpoint_directory;
static unsigned int checkpoint_threshold = 100;
static int save_gathered_checkpoints = 0;

/* There's much copy-paste here */

struct cp_info {
    bmstatus & bm;
    int level;
    unsigned int t0;
    unsigned int t1;
    unsigned int target_t;
    int mpi;
    int rank;
    std::string auxfile;
    std::string sdatafile;
    std::string gdatafile;
    std::string datafile; /* a copy of either sdatafile or gdatafile */
    /* be sure to change when needed */
    static constexpr unsigned long format = 4;
    cp_info(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1, int mpi);
    bool save_aux_file(size_t Xsize) const;
    bool load_aux_file(size_t & Xsize);
    bool checkpoint_already_present() const;
    int load_data_file(matpoly & X);
    int save_data_file(matpoly const & X);
    struct invalid_aux_file : public std::runtime_error {
        invalid_aux_file(std::string const & s) : std::runtime_error(s) {}
    };
};

constexpr unsigned long cp_info::format;

void lingen_checkpoints_decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "checkpoint-directory",
            "where to save checkpoints");
    param_list_decl_usage(pl, "checkpoint-threshold",
            "threshold for saving checkpoints");
    param_list_decl_usage(pl, "save_gathered_checkpoints",
            "save global checkpoints files, instead of per-job files");
}

void lingen_checkpoints_lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "checkpoint-directory");
    param_list_lookup_string(pl, "checkpoint-threshold");
    param_list_lookup_string(pl, "save_gathered_checkpoints");
}

void lingen_checkpoints_interpret_parameters(cxx_param_list & pl)
{
    checkpoint_directory = param_list_lookup_string(pl, "checkpoint-directory");
    param_list_parse_uint(pl, "checkpoint-threshold", &checkpoint_threshold);
    param_list_parse_int(pl, "save_gathered_checkpoints", &save_gathered_checkpoints);
}

cp_info::cp_info(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1, int mpi)
    : bm(bm), t0(t0), t1(t1), mpi(mpi)
{
    if (mpi)
        MPI_Comm_rank(bm.com[0], &(rank));
    else
        rank = 0;
    level = bm.depth();
    switch(which) {
        case LINGEN_CHECKPOINT_E:
            auxfile = fmt::sprintf("%s/E.%d.%u.aux",
                    checkpoint_directory, level, t0);
            gdatafile = fmt::sprintf("%s/E.%d.%u.single.data",
                    checkpoint_directory, level, t0);
            sdatafile = fmt::sprintf("%s/E.%d.%u.%d.data",
                    checkpoint_directory, level, t0, rank);
            break;
        case LINGEN_CHECKPOINT_PI:
            auxfile = fmt::sprintf("%s/pi.%d.%u.%u.aux",
                    checkpoint_directory, level, t0, t1);
            gdatafile = fmt::sprintf("%s/pi.%d.%u.%u.single.data",
                    checkpoint_directory, level, t0, t1);
            sdatafile = fmt::sprintf("%s/pi.%d.%u.%u.%d.data",
                    checkpoint_directory, level, t0, t1, rank);
            break;
    }
    datafile = mpi ? sdatafile : gdatafile;
}

bool cp_info::save_aux_file(size_t Xsize) const /*{{{*/
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
    os << abfield_characteristic_srcptr(bm.d.ab) << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.delta[i];
    os << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.lucky[i];
    os << "\n";
    os << bm.done;
    os << "\n";
    os << logline_serialize();
    os << "\n";
    os << bm.hints;
    os << "\n";
    os << bm.stats;
    os << "\n";
    bool ok = os.good();
    if (!ok)
        unlink(auxfile.c_str());
    return ok;
}/*}}}*/

bool cp_info::checkpoint_already_present() const/*{{{*/
{
    cp_info test = *this;
    size_t Xsize;
    try {
        if (!test.load_aux_file(Xsize))
            return false;
        int sdata_ok = access(sdatafile.c_str(), R_OK) == 0;
        int gdata_ok = rank || access(gdatafile.c_str(), R_OK) == 0;
        return sdata_ok || gdata_ok;
    } catch (cp_info::invalid_aux_file const & inv) {
        fmt::fprintf(stderr, "Overwriting bogus checkpoint file %s [%s]\n",
                datafile, inv.what());
        return false;
    }
}/*}}}*/

bool cp_info::load_aux_file(size_t & Xsize)/*{{{*/
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
        throw invalid_aux_file("checkpoint file cannot be used (version < 1)\n");
    if (hformat != format)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (version %lu < %lu)\n", hformat, format));
    unsigned int xm,xn;
    is >> xm >> xn;
    if (xm != m || xn != n)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (made for (m,n)=(%u,%u)\n", xm, xn));
    int xlevel;
    unsigned int xt0, xt1;
    is >> xlevel >> xt0 >> xt1 >> target_t;
    if (xlevel != level || t0 != xt0 || t1 != xt1)
        throw invalid_aux_file(fmt::sprintf("checkpoint file cannot be used (made for depth=%d t0=%u t1=%u\n", xlevel, xt0, xt1));
    ASSERT_ALWAYS(target_t <= t1);
    is >> Xsize;
    cxx_mpz xp;
    is >> xp;
    if (mpz_cmp(xp, abfield_characteristic_srcptr(bm.d.ab)) != 0)
        throw invalid_aux_file("checkpoint file cannot be used (made for wrong p)\n");
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
    is >> nbm.hints;
    for(auto const & x : nbm.hints) {
        if (!x.second.check()) {
            is.setstate(std::ios::failbit);
            throw invalid_aux_file("checkpoint contains invalid schedule information\n");
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
                + os.str()
                + "\n");
    }

    if (!(is >> nbm.stats))
        throw invalid_aux_file("stats in checkpoint file cannot be merged with existing stats\n");
   
    bm = std::move(nbm);

    return is.good();
}/*}}}*/

/* TODO: adapt for GF(2) */
int cp_info::load_data_file(matpoly & X)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    FILE * data = fopen(datafile.c_str(), "rb");
    int rc;
    if (data == NULL) {
        fmt::fprintf(stderr, "Warning: cannot open %s\n", datafile);
        return 0;
    }
    rc = matpoly_read(ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
/* I think we always have Xsize == X.size, the only questionable
 * situation is when we're saving part of a big matrix */
int cp_info::save_data_file(matpoly const & X)/*{{{*/
{
    abdst_field ab = bm.d.ab;
    std::ofstream data(datafile, std::ios_base::out | std::ios_base::binary);
    int rc;
    if (!data) {
        fmt::fprintf(stderr, "Warning: cannot open %s\n", datafile);
        unlink(auxfile.c_str());
        return 0;
    }
    rc = matpoly_write(ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) goto cp_info_save_data_file_bailout;
    if (data.good()) return 1;
cp_info_save_data_file_bailout:
    unlink(datafile.c_str());
    unlink(auxfile.c_str());
    return 0;
}/*}}}*/

template<>
int load_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;

    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;

    cp_info cp(bm, which, t0, t1, 0);

    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = 1;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (cp_info::invalid_aux_file const & inv) {
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
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    cp_info cp(bm, which, t0, t1, 0);
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

int load_mpi_checkpoint_file_scattered(bmstatus & bm, cp_which which, bigmatpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = bm.d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, which, t0, t1, 1);
    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    int ok = 1, error = 0;
    try {
        ok = cp.load_aux_file(Xsize);
    } catch (cp_info::invalid_aux_file const & inv) {
        if (!rank)
            fmt::fprintf(stderr, "Error with checkpoint file %s:\n%s\n", cp.datafile, inv.what());
        error = 1;
    }
    MPI_Bcast(&error, 1, MPI_INT, 0, bm.com[0]);
    ASSERT_ALWAYS(!error);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (!ok) {
        if (!rank)
            fmt::fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
        return false;
    }
    MPI_Bcast(&Xsize, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&bm.delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
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
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    cp_info cp(bm, which, t0, t1, 1);
    int ok;
    ok = cp.checkpoint_already_present();
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (ok) return 1;

    ok = cp.save_aux_file(X.get_size());
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
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
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, which, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    size_t Xsize;
    int ok = 1, error = 0;
    try {
        cp.load_aux_file(Xsize);
    } catch (cp_info::invalid_aux_file const & inv) {
        error = 1;
    }
    MPI_Bcast(&error, 1, MPI_INT, 0, bm.com[0]);
    ASSERT_ALWAYS(!error);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (!ok)
        /* This can be a normal situation if the .aux file does not
         * exist.
         */
        return false;
    MPI_Bcast(&Xsize, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&bm.delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, gathered)",
            cp.datafile.c_str());
    do {
        FILE * data = NULL;
        if (!rank) ok = (data = fopen(cp.datafile.c_str(), "rb")) != NULL;
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
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

        double avg = average_matsize(ab, m + n, m + n, 0);
        unsigned int B = iceildiv(io_matpoly_block_size, avg);

        /* This is only temp storage ! */
        matpoly loc(ab, m + n, m + n, B);
        loc.zero_pad(B);

        for(unsigned int k = 0 ; ok && k < X.get_size() ; k += B) {
            unsigned int nc = MIN(B, X.get_size() - k);
            if (!rank)
                ok = matpoly_read(ab, data, loc, 0, nc, 0, 0) == (int) nc;
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            X.scatter_mat_partial(loc, 0, k, nc);
        }

        if (!rank) {
            int rc = fclose(data);
            ok = ok && (rc == 0);
        }
    } while (0);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
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
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, which, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    int ok;
    ok = cp.checkpoint_already_present();
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
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
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fmt::fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                break;
            }

            double avg = average_matsize(ab, m + n, m + n, 0);
            unsigned int B = iceildiv(io_matpoly_block_size, avg);

            /* This is only temp storage ! */
            matpoly loc(ab, m + n, m + n, B);
            loc.zero_pad(B);

            for(unsigned int k = 0 ; ok && k < X.get_size() ; k += B) {
                unsigned int nc = MIN(B, X.get_size() - k);
                X.gather_mat_partial(loc, 0, k, nc);
                if (!rank)
                    ok = matpoly_write(ab, data, loc, 0, nc, 0, 0) == (int) nc;
                MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            }

            if (!rank) {
                data.close();
                ok = ok && (bool) data;
            }
        } while (0);
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
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
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    cp_info cp(bm, which, t0, t1, 1);
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
    if (save_gathered_checkpoints) {
        return save_mpi_checkpoint_file_gathered(bm, which, X, t0, t1);
    } else {
        return save_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
    }
}/*}}}*/



