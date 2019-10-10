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
    int mpi;
    int rank;
    char * auxfile;
    char * sdatafile;
    char * gdatafile;
    const char * datafile;
    /* be sure to change when needed */
    static constexpr unsigned long format = 3;
    FILE * aux;
    FILE * data;
    cp_info(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1, int mpi);
    ~cp_info();
    bool save_aux_file(size_t Xsize) const;
    bool load_aux_file(size_t & Xsize);
    int load_data_file(matpoly & X, size_t Xsize);
    int save_data_file(matpoly const & X, size_t Xsize);
};

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
    int rc;
    level = bm.depth();
    switch(which) {
        case LINGEN_CHECKPOINT_E:
            rc = asprintf(&auxfile, "%s/E.%d.%u.aux",
                    checkpoint_directory, level, t0);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&gdatafile, "%s/E.%d.%u.single.data",
                    checkpoint_directory, level, t0);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&sdatafile, "%s/E.%d.%u.%d.data",
                    checkpoint_directory, level, t0, rank);
            ASSERT_ALWAYS(rc >= 0);
            break;
        case LINGEN_CHECKPOINT_PI:
            rc = asprintf(&auxfile, "%s/pi.%d.%u.%u.aux",
                    checkpoint_directory, level, t0, t1);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&gdatafile, "%s/pi.%d.%u.%u.single.data",
                    checkpoint_directory, level, t0, t1);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&sdatafile, "%s/pi.%d.%u.%u.%d.data",
                    checkpoint_directory, level, t0, t1, rank);
            ASSERT_ALWAYS(rc >= 0);
            break;
    }
    datafile = mpi ? sdatafile : gdatafile;
}

cp_info::~cp_info()
{
    free(sdatafile);
    free(gdatafile);
    free(auxfile);
}

bool cp_info::save_aux_file(size_t Xsize) const /*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    if (rank) return 1;
    std::ofstream os(auxfile);
    os << "format " << format << "\n";
    os << Xsize << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.delta[i];
    os << "\n";
    for(unsigned int i = 0 ; i < m + n ; i++) os << " " << bm.lucky[i];
    os << "\n";
    os << bm.done;
    os << "\n";
    os << logline_query_timer();
    os << "\n";
    os << bm.hints;
    os << "\n";
    os << bm.stats;
    os << "\n";
    bool ok = os.good();
    if (!ok)
        unlink(auxfile);
    return ok;
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
    if (hfstring != "format") {
        fprintf(stderr, "Warning: checkpoint file cannot be used (version < 1)\n");
        return false;
    }
    if (hformat != format) {
        fprintf(stderr, "Warning: checkpoint file cannot be used (version %lu < %lu)\n", hformat, format);
        return false;
    }

    is >> Xsize;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        is >> nbm.delta[i];
    }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        is >> nbm.lucky[i];
    }
    is >> nbm.done;
    double tt;
    is >> tt;
    logline_reset_timer(tt);
    is >> nbm.hints;
    for(auto const & x : nbm.hints) {
        if (!x.second.check()) {
            fprintf(stderr, "Warning: checkpoint contains invalid schedule information\n");
            is.setstate(std::ios::failbit);
            return false;
        }
    }

    if (bm.hints != nbm.hints) {
        is.setstate(std::ios::failbit);
        fprintf(stderr, "Warning: checkpoint file cannot be used since it was made for another set of schedules (stats would be incorrect)\n");
        std::stringstream os;
        os << bm.hints;
        fprintf(stderr, "textual description of the schedule set that we expect to find:\n%s\n", os.str().c_str());
        return false;
    }

    if (!(is >> nbm.stats))
        return false;
   
    bm = std::move(nbm);

    return is.good();
}/*}}}*/

/* TODO: adapt for GF(2) */
int cp_info::load_data_file(matpoly & X, size_t Xsize)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    FILE * data = fopen(datafile, "rb");
    int rc;
    if (data == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", datafile);
        return 0;
    }
    X = matpoly(ab, m+n, m+n, Xsize);
    X.set_size(Xsize);
    rc = matpoly_read(ab, data, X, 0, X.get_size(), 0, 0);
    if (rc != (int) X.get_size()) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
/* I think we always have Xsize == X.size, the only questionable
 * situation is when we're saving part of a big matrix */
int cp_info::save_data_file(matpoly const & X, size_t Xsize)/*{{{*/
{
    abdst_field ab = bm.d.ab;
    std::ofstream data(datafile, std::ios_base::out | std::ios_base::binary);
    int rc;
    if (!data) {
        fprintf(stderr, "Warning: cannot open %s\n", datafile);
        unlink(auxfile);
        return 0;
    }
    rc = matpoly_write(ab, data, X, 0, Xsize, 0, 0);
    if (rc != (int) X.get_size()) goto cp_info_save_data_file_bailout;
    if (data.good()) return 1;
cp_info_save_data_file_bailout:
    unlink(datafile);
    unlink(auxfile);
    return 0;
}/*}}}*/

template<>
int load_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;

    cp_info cp(bm, which, t0, t1, 0);

    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = cp.load_aux_file(Xsize);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading %s", cp.datafile);
        ok = cp.load_data_file(X, Xsize);
        logline_end(&bm.t_cp_io,"");
        if (!ok)
            fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
    }
    if (ok) bm.t = t1;
    return ok;
}/*}}}*/

template<>
int save_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly const & X, unsigned int t0, unsigned int t1)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    cp_info cp(bm, which, t0, t1, 0);
    logline_begin(stdout, SIZE_MAX, "Saving %s%s",
            cp.datafile,
            cp.mpi ? " (MPI, scattered)" : "");
    int ok = cp.save_aux_file(X.get_size());
    if (ok) ok = cp.save_data_file(X, X.get_size());
    logline_end(&bm.t_cp_io,"");
    if (!ok && !cp.rank)
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
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
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, which, t0, t1, 1);
    ASSERT_ALWAYS(X.check_pre_init());
    size_t Xsize;
    int ok = cp.load_aux_file(Xsize);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&Xsize, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&bm.delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, scattered)",
                cp.datafile);
        do {
            FILE * data = fopen(cp.datafile, "rb");
            int rc;
            ok = data != NULL;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                if (data) free(data);
                break;
            }
            X.finish_init(ab, m+n, m+n, Xsize);
            X.set_size(Xsize);
            rc = matpoly_read(ab, data, X.my_cell(), 0, X.get_size(), 0, 0);
            ok = ok && rc == (int) X.get_size();
            rc = fclose(data);
            ok = ok && (rc == 0);
        } while (0);
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        logline_end(&bm.t_cp_io,"");
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp.datafile);
        }
    } else if (!rank) {
        fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
    }
    if (ok) bm.t = t1;
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
    int ok = cp.save_aux_file(X.get_size());
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (!ok && !rank) unlink(cp.auxfile);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Saving %s (MPI, scattered)",
                cp.datafile);
        ok = cp.save_data_file(X.my_cell(), X.get_size());
        logline_end(&bm.t_cp_io,"");
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok) {
            if (cp.datafile) unlink(cp.datafile);
            if (!rank) unlink(cp.auxfile);
        }
    }
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
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
    int ok = cp.load_aux_file(Xsize);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&Xsize, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&bm.delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading %s (MPI, gathered)",
                cp.datafile);
        do {
            FILE * data = NULL;
            if (!rank) ok = (data = fopen(cp.datafile, "rb")) != NULL;
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                if (data) free(data);
                break;
            }

            X.finish_init(ab, m+n, m+n, Xsize);
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
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp.datafile);
        }
    } else if (!rank) {
        fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
    }
    if (ok) bm.t = t1;
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
    logline_begin(stdout, SIZE_MAX, "Saving %s (MPI, gathered)",
            cp.datafile);
    int ok = cp.save_aux_file(X.get_size());
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
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
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
            if (cp.datafile) unlink(cp.datafile);
            unlink(cp.auxfile);
        }
    }
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
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
    int aux_ok = rank || access(cp.auxfile, R_OK) == 0;
    int sdata_ok = access(cp.sdatafile, R_OK) == 0;
    int scattered_ok = aux_ok && sdata_ok;
    MPI_Allreduce(MPI_IN_PLACE, &scattered_ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (scattered_ok) {
        ok = load_mpi_checkpoint_file_scattered(bm, which, X, t0, t1);
        if (ok) return ok;
    }
    int gdata_ok = rank || access(cp.gdatafile, R_OK) == 0;
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



