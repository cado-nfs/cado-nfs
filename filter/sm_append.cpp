/* Shirokauer maps 

   This is largely copied from sm_simple.cpp ; however the purpose of
   sm_simple was to be kept simple and stupid, so let's keep it like
   this.

   This program is to be used as an mpi accelerator for reconstructlog-dl

   Given a relation file where each line begins with an a,b pair IN
   HEXADECIMAL FORMAT, output the same line, appending the SM values in
   the end.

   SM computation is offloaded to the (many) MPI jobs.

*/

#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <vector>
#include <string>
#include <cstdarg>     // for va_end, va_list, va_start
#include <cstdint>     // for int64_t, uint64_t
#include <iosfwd>       // for std
#include <memory>       // for allocator_traits<>::value_type
#include <gmp.h>
#include "cado_poly.h"  // for NB_POLYS_MAX, cado_poly_clear, cado_poly_init
#include "gzip.h"       // fopen_maybe_compressed
#include "macros.h"
#include "mpz_poly.h"   // for mpz_poly_clear, mpz_poly_init, mpz_poly, mpz_...
#include "params.h"
#include "select_mpi.h"
#include "sm_utils.h"   // sm_side_info
#include "timing.h"     // seconds
#include "verbose.h"    // verbose_output_print

using namespace std;

struct ab_pair {
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
};

typedef vector<ab_pair> ab_pair_batch;

unsigned int batch_size = 128;
unsigned int debug = 0;

struct task_globals {
    int nsm_total;
    size_t limbs_per_ell;
    FILE * in;
    FILE * out;
    size_t nrels_in;
    size_t nrels_out;
};

struct peer_status {
    MPI_Request req;
    ab_pair_batch batch;
    vector<string> rels;
    void receive(task_globals & tg, int peer, int turn);
    void send_finish(task_globals & tg, int peer, int turn);
    /* returns 1 on eof, 0 normally */
    int create_and_send_batch(task_globals& tg, int peer, int turn);
};

static int debug_fprintf(FILE * out, const char * fmt, ...)
    ATTR_PRINTF(2, 3);
static int debug_fprintf(FILE * out, const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int rc = debug ? vfprintf(out, fmt, ap) : 1;
    va_end(ap);
    return rc;
}

void peer_status::receive(task_globals & tg, int peer, int turn)
{
    unsigned long bsize = batch.size();

    if (!bsize) return;

    MPI_Wait(&req, MPI_STATUS_IGNORE);
    batch.clear();

    mp_limb_t * returns = new mp_limb_t[bsize * tg.nsm_total * tg.limbs_per_ell];
    // [bsize][tg.nsm_total][tg.limbs_per_ell];

    MPI_Recv(returns, bsize * tg.nsm_total * tg.limbs_per_ell * sizeof(mp_limb_t), MPI_BYTE, peer, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    ASSERT_ALWAYS(rels.size() == bsize);

    for(unsigned long i = 0 ; i < bsize ; i++) {
        fputs(rels[i].c_str(), tg.out);
        bool comma = false;
        for (int j = 0 ; j < tg.nsm_total ; j++) {
            mp_limb_t * rij = returns + ((i * tg.nsm_total) + j) * tg.limbs_per_ell;
            gmp_fprintf(tg.out, "%c%Nd", comma ? ',' : ':', rij, tg.limbs_per_ell);
            comma=true;
        }
        fputc('\n', tg.out);
        tg.nrels_out++;
    }
    rels.clear();

    delete[] returns;
}

void peer_status::send_finish(task_globals &, int peer, int turn)
{
    /* tell this slave to finish */
    unsigned long bsize = 0;
    MPI_Send(&bsize, 1, MPI_UNSIGNED_LONG, peer, turn, MPI_COMM_WORLD);
}

int peer_status::create_and_send_batch(task_globals& tg, int peer, int turn)
{
    int eof = 0;
    char buf[1024];

    ASSERT_ALWAYS(batch.empty());
    while (!eof && batch.size() < batch_size && fgets(buf, 1024, tg.in)) {
        int n = strlen(buf);
        if (!n) {
            fprintf(stderr, "Got 0-sized buffer in fgets, shouldn't happen. Assuming EOF.\n");
            eof=true;
            break;
        }
        buf[n-1]='\0';

        if (buf[0] == '#') {
            fputs(buf, tg.out);
            fputc('\n', tg.out);
            continue;
        }

        tg.nrels_in++;
        rels.push_back(string(buf));

        char * p = buf;
        ab_pair ab;
        int64_t sign = 1;
        if (*p=='-') {
            sign=-1;
            p++;
        }
        if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &ab.a, &ab.b) < 2) {
            fprintf(stderr, "Parse error at line: %s\n", buf);
            exit(EXIT_FAILURE);
        }
        ab.a *= sign;

        batch.push_back(ab);
    }
    if (!eof && batch.size() < batch_size) {
        eof=true;
        if (ferror(stdin)) {
            fprintf(stderr, "Error on stdin\n");
        }
    }
    /* 0 bsize will be recognized by slaves as a 
     * reason to stop processing */
    unsigned long bsize = batch.size();
    MPI_Send(&bsize, 1, MPI_UNSIGNED_LONG, peer, turn, MPI_COMM_WORLD);
    if (bsize)
        MPI_Isend((char*) batch.data(), bsize * sizeof(ab_pair), MPI_BYTE, peer, turn, MPI_COMM_WORLD, &req);
    return eof;
}

#define CSI_RED "\033[00;31m"
#define CSI_GREEN "\033[00;32m"
#define CSI_YELLOW "\033[00;33m"
#define CSI_BLUE "\033[00;34m"
#define CSI_PINK "\033[00;35m"
#define CSI_BOLDRED "\033[01;31m"
#define CSI_BOLDGREEN "\033[01;32m"
#define CSI_BOLDYELLOW "\033[01;33m"
#define CSI_BOLDBLUE "\033[01;34m"
#define CSI_BOLDPINK "\033[01;35m"
#define CSI_RESET "\033[m"

static void sm_append_master(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys, int size)
{
    /* need to know how many mp_limb_t's we'll get back from each batch */
    task_globals tg;
    tg.limbs_per_ell = 0;
    tg.nsm_total=0;
    tg.in = in;
    tg.out = out;
    tg.nrels_in = 0;
    tg.nrels_out = 0;
    for(int side = 0; side < nb_polys; side++) {
        tg.nsm_total += sm_info[side]->nsm;
        if (sm_info[side]->nsm) 
            tg.limbs_per_ell = mpz_size(sm_info[side]->ell);
    }

    std::vector<peer_status> peers(size);

    int eof = 0;
    /* eof = 1 on first time. eof = 2 when all receives are done */

    fprintf(stderr, "# running master with %d slaves, batch size %u\n",
            size-1, batch_size);
    fprintf(stderr, "# make sure you use \"--bind-to core\" or equivalent\n");

    double t0 = wct_seconds();
    int turn;
    for(turn = 0 ; eof <= 2 ; turn++, eof += !!eof) {
        double t = wct_seconds();
        debug_fprintf(stderr, "%.3f " CSI_BOLDRED "start turn %d" CSI_RESET "\n", t0, turn);
        for(int peer = 1; peer < size; peer++) {
            if (eof && peers[peer].batch.empty()) {
                /* Our last send was a 0-send, so we have nothing to do */
                continue;
            }

            double dt = wct_seconds();
            debug_fprintf(stderr, "%.3f start turn %d receive from peer %d\n", wct_seconds(), turn - 1, peer);
            peers[peer].receive(tg, peer, turn - 1);
            dt = wct_seconds() - dt;
            debug_fprintf(stderr, "%.3f done turn %d receive from peer %d [taken %.1f]\n", wct_seconds(), turn - 1, peer, dt);

            if (eof) {
                debug_fprintf(stderr, "%.3f start turn %d send finish to peer %d\n", wct_seconds(), turn, peer);
                peers[peer].send_finish(tg, peer, turn);
            } else if (!eof) {
                dt = wct_seconds();
                debug_fprintf(stderr, "%.3f start turn %d send to peer %d\n", wct_seconds(), turn, peer);
                eof = peers[peer].create_and_send_batch(tg, peer, turn);
                dt = wct_seconds() - dt;
                debug_fprintf(stderr, "%.3f done turn %d send to peer %d [taken %.1f]\n", wct_seconds(), turn, peer, dt);
            }
        }
        debug_fprintf(stderr, "%.3f " CSI_BOLDRED "done turn %d " CSI_RESET "[taken %.1f] s\n", wct_seconds(), turn, wct_seconds()-t);
        if (turn && !(turn & (turn+1))) {
            /* print only when turn is a power of two */
            fprintf(stderr, "# printed %zu rels in %.1f s"
                    " (%.1f / batch, %.1f rels/s)\n",
                    tg.nrels_out, wct_seconds()-t0,
                    (wct_seconds()-t0) / turn,
                    tg.nrels_out / (wct_seconds()-t0));
        }
    }
    fprintf(stderr, "# final: printed %zu rels in %.1f s"
            " (%.1f / batch, %.1f rels/s)\n",
            tg.nrels_out, wct_seconds()-t0,
            (wct_seconds()-t0) / turn,
            tg.nrels_out / (wct_seconds()-t0));
}

static void sm_append_slave(sm_side_info *sm_info, int nb_polys)
{
    /* need to know how many mp_limb_t's we'll get back from each batch */
    size_t limbs_per_ell = 0;
    int nsm_total=0;
    int maxdeg = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(int side = 0; side < nb_polys; side++) {
        nsm_total += sm_info[side]->nsm;
        maxdeg = MAX(maxdeg, sm_info[side]->f->deg);
        if (sm_info[side]->nsm) limbs_per_ell = mpz_size(sm_info[side]->ell);
    }


    for(int turn = 0 ; ; turn++) {
        unsigned long bsize;
        MPI_Recv(&bsize, 1, MPI_UNSIGNED_LONG, 0, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (bsize == 0) {
            debug_fprintf(stderr, "%.3f turn %d peer %d receive finish\n", wct_seconds(), turn, rank);
            break;
        }
        ab_pair_batch batch(bsize);
        MPI_Recv((char*) batch.data(), bsize * sizeof(ab_pair), MPI_BYTE, 0, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double t0 = wct_seconds();
        debug_fprintf(stderr, "%.3f turn %d peer %d start on batch of size %lu\n", wct_seconds(), turn, rank, bsize);
        mp_limb_t * returns = new mp_limb_t[bsize * nsm_total * limbs_per_ell];
        memset(returns, 0, bsize*nsm_total*limbs_per_ell*sizeof(mp_limb_t));

#ifdef HAVE_OPENMP
// #pragma omp parallel
#endif
        {
            cxx_mpz_poly smpol(maxdeg), pol(1);
#ifdef HAVE_OPENMP
// #pragma omp for
#endif
            for(unsigned long i = 0 ; i < bsize ; i++) {
                mpz_poly_setcoeff_int64(pol, 0, batch[i].a);
                mpz_poly_setcoeff_int64(pol, 1, -(int64_t) batch[i].b);
                int smidx = 0;
                for (int side = 0; side < nb_polys; ++side) {
                    compute_sm_piecewise(smpol, pol, sm_info[side]);
                    for(int k = 0 ; k < sm_info[side]->nsm ; k++, smidx++) {
                        if (k <= smpol->deg) {
                            mp_limb_t * rix = returns + (i * nsm_total + smidx) * limbs_per_ell;
                            for(size_t j = 0 ; j < limbs_per_ell ; j++) {
                                rix[j] = mpz_getlimbn(smpol->coeff[k], j);
                            }
                        }
                    }
                }
            }
        }
        debug_fprintf(stderr, "%.3f " CSI_BLUE "turn %d peer %d done batch of size %lu" CSI_RESET " [taken %.1f]\n", wct_seconds(), turn, rank, bsize, wct_seconds() - t0);
        if (rank == 1 && turn == 2)
            fprintf(stderr, "# peer processes batch of %lu in %.1f [%.1f SMs/s]\n",
                    bsize,
                    wct_seconds() - t0,
                    bsize / (wct_seconds() - t0));

        t0 = wct_seconds();
        MPI_Send(returns, bsize * nsm_total * limbs_per_ell * sizeof(mp_limb_t), MPI_BYTE, 0, turn, MPI_COMM_WORLD);
        delete[] returns;
        debug_fprintf(stderr, "%.3f turn %d peer %d send return took %.1f\n", wct_seconds(), turn, rank, wct_seconds() - t0);
    }
}

static void sm_append_sync(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys)
{
    char buf[1024];
    mpz_poly pol, smpol;
    int maxdeg = sm_info[0]->f->deg;
    for(int side = 1; side < nb_polys; side++)
        maxdeg = MAX(maxdeg, sm_info[side]->f->deg);
    mpz_poly_init(pol, maxdeg);
    mpz_poly_init(smpol, maxdeg);
    while (fgets(buf, 1024, in)) {
        int n = strlen(buf);
        if (!n) break;
        buf[n-1]='\0';

        if (buf[0] == '#') {
            fputs(buf, out);
            fputc('\n', out);
            continue;
        }

        char * p = buf;
        int64_t a;		/* only a is allowed to be negative */
        uint64_t b;
        int64_t sign = 1;
        if (*p=='-') {
            sign=-1;
            p++;
        }
        if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &a, &b) < 2) {
            fprintf(stderr, "Parse error at line: %s\n", buf);
            exit(EXIT_FAILURE);
        }

        mpz_poly_init_set_ab(pol, a*sign, b);

        fputs(buf, out);
        fputc(':', out);
        for (int side = 0; side < nb_polys; ++side) {
            compute_sm_piecewise(smpol, pol, sm_info[side]);
            print_sm2(out, sm_info[side], smpol, ",");
            if (side == 0 && sm_info[0]->nsm > 0 && sm_info[1]->nsm > 0)
                fputc(',', out);
        }
        fputc('\n', out);
        mpz_poly_clear(pol);
    }
    mpz_poly_clear(smpol);
}


static void sm_append(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys)
{
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size > 1) {
        if (rank == 0) {
            sm_append_master(in, out, sm_info, nb_polys, size);
        } else {
            sm_append_slave(sm_info, nb_polys);
        }
    } else {
        sm_append_sync(in, out, sm_info, nb_polys);
    }
}


static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "nsm", "number of SMs to use per side");
    param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
    param_list_decl_usage(pl, "in", "data input (defaults to stdin)");
    param_list_decl_usage(pl, "out", "data output (defaults to stdout)");
    param_list_decl_usage(pl, "b", "batch size for MPI loop");
    verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit (EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

// coverity[root_function]
int main (int argc, char **argv)
{
    MPI_Init(&argc, & argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char *argv0 = argv[0];

    const char *polyfile = NULL;

    param_list pl;
    cado_poly cpoly;

    mpz_t ell;

    /* read params */
    param_list_init(pl);
    declare_usage(pl);

    if (argc == 1)
        usage (argv[0], NULL, pl);

    argc--,argv++;
    for ( ; argc ; ) {
        if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
        usage (argv0, NULL, pl);
    }

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read ell from command line (assuming radix 10) */
    mpz_init (ell);
    if (!param_list_parse_mpz(pl, "ell", ell)) {
        fprintf(stderr, "Error: parameter -ell is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    cado_poly_init (cpoly);
    cado_poly_read(cpoly, polyfile);

    std::vector<mpz_poly_srcptr> F(cpoly->nb_polys, NULL);

    for(int side = 0; side < cpoly->nb_polys; side++)
        F[side] = cpoly->pols[side];

    std::vector<int> nsm_arg(cpoly->nb_polys, -1);
    param_list_parse_int_args_per_side(pl, "nsm",
            nsm_arg.data(), cpoly->nb_polys,
            ARGS_PER_SIDE_DEFAULT_AS_IS);

    FILE * in = rank ? NULL : stdin;
    FILE * out = rank ? NULL: stdout;
    const char * infilename = param_list_lookup_string(pl, "in");
    const char * outfilename = param_list_lookup_string(pl, "out");

    if (!rank && infilename) {
        in = fopen_maybe_compressed(infilename, "r");
        ASSERT_ALWAYS(in != NULL);
    }
    if (!rank && outfilename) {
        out = fopen_maybe_compressed(outfilename, "w");
        ASSERT_ALWAYS(out != NULL);
    }

    param_list_parse_uint(pl, "b", &batch_size);

    verbose_interpret_parameters(pl);

    const char * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

    if (param_list_warn_unused(pl))
        usage (argv0, NULL, pl);

    if (!rank)
        param_list_print_command_line (stdout, pl);

    sm_side_info * sm_info = new sm_side_info[cpoly->nb_polys];

    for(int side = 0 ; side < cpoly->nb_polys; side++) {
        sm_side_info_init(sm_info[side], F[side], ell);
        sm_side_info_set_mode(sm_info[side], sm_mode_string);
        if (nsm_arg[side] >= 0)
            sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
        if (!rank)
            printf("# Using %d SMs on side %d\n", sm_info[side]->nsm, side);
    }

    /*
       if (!rank) {
       for (int side = 0; side < cpoly->nb_polys; side++) {
       printf("\n# Polynomial on side %d:\nF[%d] = ", side, side);
       mpz_poly_fprintf(stdout, F[side]);

       printf("# SM info on side %d:\n", side);
       sm_side_info_print(stdout, sm_info[side]);

       fflush(stdout);
       }
       }
       */

    sm_append(in, out, sm_info, cpoly->nb_polys);

    /* Make sure we print no footer line, because reconstructlog-dl won't
     * grok it */
    if (!rank) {
        fflush(stdout);
    }

    if (!rank && infilename) fclose_maybe_compressed(in, infilename);
    if (!rank && out != stdout) fclose_maybe_compressed(out, outfilename);

    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
        sm_side_info_clear(sm_info[side]);
    }
    delete[] sm_info;

    mpz_clear(ell);
    cado_poly_clear(cpoly);
    param_list_clear(pl);

    MPI_Finalize();

    return 0;
}
