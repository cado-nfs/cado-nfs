/*
 * Program: crtalgsqrt
 * Authors: E. Thomé.
 * Purpose: computing the squareroots and finishing the factorization
 *
 */

/* TODO list.
 *
 * This program seems to work, but is not complete.
 *
 * For correctness:
 *
 * - finish binding with the rest: check also the rational square root,
 *   and produce the factorization. Not necessarily independent from the
 *   above, since the elementary check is whether we get x^2=a mod N.
 * - update this TODO list.
 *
 * For speed / memory:
 *
 * - The program is probably overzealous with mpz_mod's sometimes. Some
 *   can be s(h)aved.
 * - The peak memory usage is apparently quite high. It is linked to the
 *   fact that we may have several threads doing full-length
 *   multiplications. Notwithstanding this dominant effect for scratch
 *   space usage, the amount of useful data should not exceed 3 times the
 *   ram_gb parameter.
 * - parallelize multiply_all_shares. It's not hard. Presently it
 *   accounts for 40% of the WCT, which is not acceptable.
 *
 * Important in terms of functionality, but not critical:
 *
 * - Use fpLLL for solving the knapsack. Should make it possible to go to
 *   12 primes or so in degree 6 (at least).
 * - Be a lot more flexible with regard to MPI setup. In particular, it
 *   should be possible to fix the s,t parameters not exactly equal to
 *   their computed minimal value.
 *
 * Not important:
 *
 * - understand and, if relevant, fix the memory leak diagnosed by
 *   valgrind. I don't understand.
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <complex.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>
#include <ctype.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "utils.h"
#include "modul_poly.h"
#include "powers_of_p.h"
#include "polyroots.h"
#include "knapsack.h"

#include "select_mpi.h"
#include "gmp-hacks.h"          // TODO: REMOVE !

/* {{{ time */
double program_starttime;
#define WCT     (wct_seconds() - program_starttime)

#define STOPWATCH_DECL	        					\
    double t0, t1;							\
    double w0, w1;							\
    double rate

#define STOPWATCH_GO()	        					\
    t0 = seconds();							\
    w0 = WCT;

#define STOPWATCH_GET()							\
    t1 = seconds();							\
    w1 = WCT;								\
    rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;

#define log_step(step) logprint("%s%s\n", __func__, step)
#define log_begin() log_step(" starts")
#define log_step_time(step)                                             \
    logprint("%s%s. %.2lf, wct %.2lf, %.1f%%\n",                        \
            __func__, step, t1-t0, w1-w0, rate)
#define log_end() log_step_time(" ends")

/* }}} */

static int verbose = 0;

double print_delay = 1;
double ram_gb = 3.0;    // Number of gigabytes. Note that this is the
// maximum size of product that are obtained.
// The actual memory footprint can be quite
// considerably larger due to FFT allocation (by
// a constant factor, though).

#define MPI_MY_SIZE_T   MPI_UNSIGNED_LONG
#define MPI_MY_UINT64_T   MPI_UNSIGNED_LONG
#define MPI_MY_INT64_T   MPI_LONG
#define MPI_MY_MP_SIZE_T        MPI_LONG
#define MPI_MY_MP_LIMB_T        MPI_UNSIGNED_LONG
#define MPI_MY_GMP_INTERNAL_SIZE_FIELD_T MPI_INT

static void usage()
{
    fprintf(stderr, "usage: crtalgsqrt algdepfile ratdepfile polyfile\n");
    exit(1);
}

/* {{{ logging */
int max_loglevel=99;
char prefix[20]={'\0'};

int logprint(const char * fmt, ...)
{
    va_list ap;
    int level=0;
    int s = strlen(fmt);
    const char * pfmt = fmt;
    if (s >= 3 && fmt[0] == '<' && isdigit(fmt[1]) && fmt[2]=='>') {
        level=fmt[1]-'0';
        pfmt += 3;
    }
    if (level > max_loglevel)
        return 0;

    static pthread_mutex_t obuf_lock = PTHREAD_MUTEX_INITIALIZER;
    static size_t st = 0;
    static char * t = NULL;
    pthread_mutex_lock(&obuf_lock);

    size_t wt = strlen(prefix) + strlen(fmt) + 80;
    if (wt > st) {
        t = realloc(t, wt);
        st = wt;
    }
    snprintf(t, st, "# [%2.2lf] %s%s", WCT, prefix, pfmt);
    
    va_start(ap, fmt);
    int rc = vfprintf(stderr, t, ap);
    va_end(ap);

    pthread_mutex_unlock(&obuf_lock);

    return rc;
}
/* }}} */

/* {{{ wrappers for some gmp operations, so as to report timings */
// above this threshold, we report each multiplication we do.
#define MUL_REPORT_THRESHOLD    8000000

#define REPORT_THIS(na, nb)     \
    (((na) > 10) && ((nb) > 10) && ((na) + (nb) > MUL_REPORT_THRESHOLD))

static void WRAP_mpz_mul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_mul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

#if 0
static void WRAP_mpz_invert(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_invert(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_inv %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
#endif

static void WRAP_mpz_addmul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_addmul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_submul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_submul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_barrett_mod(mpz_ptr c, mpz_srcptr a, mpz_srcptr p, mpz_srcptr q)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(p);
    barrett_mod(c,a,p,q);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb) && na > nb + 10) {
        logprint("<9> mpz_mod %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_mod(mpz_ptr c, mpz_srcptr a, mpz_srcptr p)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(p);
    mpz_mod(c,a,p);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb) && na > nb + 10) {
        logprint("<9> mpz_mod %zu %zu (%.1f) %.1f %.1f (%.1f%%)\n", na, nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
/* }}} */

/* {{{ cache files */
int rcache = 1;
int wcache = 1;
#define CACHEDIR        "/tmp"
#define CACHEPREFIX        "CRTALGSQRT."

struct cachefile_s {
    char basename[128];
    FILE * f;
    int writing;
};

typedef struct cachefile_s cachefile[1];
typedef struct cachefile_s * cachefile_ptr;

void cachefile_vinit(cachefile_ptr c, const char * fmt, va_list ap)
{
    memset(c, 0, sizeof(c));
    vsnprintf(c->basename, sizeof(c->basename), fmt, ap);
}

void cachefile_init(cachefile_ptr c, const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    cachefile_vinit(c, fmt, ap);
    va_end(ap);
}

int cachefile_open_w(cachefile_ptr c)
{
    char tname[128];
    snprintf(tname, sizeof(tname), CACHEDIR "/.pre." CACHEPREFIX "%s", c->basename);
    c->f = fopen(tname, "w");
    c->writing = 0;
    if (c->f == NULL) {
        fprintf(stderr, "%s: %s\n", tname, strerror(errno));
    }
    c->writing = 1;
    return c->f != NULL;
}

int cachefile_open_r(cachefile_ptr c)
{
    char tname[128];
    snprintf(tname, sizeof(tname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    c->f = fopen(tname, "r");
    c->writing = 0;
    return c->f != NULL;
}

int cachefile_exists(const char * fmt, ...)
{
    cachefile c;
    char tname[128];
    va_list ap;
    va_start(ap, fmt);
    cachefile_vinit(c, fmt, ap);
    va_end(ap);
    snprintf(tname, sizeof(tname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    return access(tname, R_OK) == 0;
}

void cachefile_close(cachefile_ptr c)
{
    char tname[128];
    char fname[128];
    fclose(c->f);
    if (c->writing == 0)
        return;
    snprintf(tname, sizeof(tname), CACHEDIR "/.pre." CACHEPREFIX "%s", c->basename);
    snprintf(fname, sizeof(tname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    rename(tname, fname);
}
/* }}} */

void
polymodF_mul_monic (poly_t Q, const poly_t P1, const poly_t P2,
        const poly_t F);

/* {{{ mpi-gmp helpers */
static int mpi_data_agrees(void *buffer, int count, MPI_Datatype datatype,/*{{{*/
        MPI_Comm comm)
{
    int s;
    MPI_Type_size(datatype, &s);
    void * b_and = malloc(count * s);
    void * b_or = malloc(count * s);
    MPI_Allreduce(buffer, b_and, count, datatype, MPI_BAND, comm);
    MPI_Allreduce(buffer, b_or, count, datatype, MPI_BOR, comm);
    int ok = memcmp(b_and, b_or, count * s) == 0;
    free(b_and);
    free(b_or);
    return ok;
}/*}}}*/
static void broadcast_mpz(mpz_ptr z, int root, MPI_Comm comm) /*{{{*/
{
    // that's much, much uglier than it should.
    int k;
    MPI_Comm_rank(comm, &k);
    mp_size_t nlimbs = mpz_size(z);
    MPI_Bcast(&nlimbs, 1, MPI_MY_MP_SIZE_T, root, comm);
    if (k != root) _mpz_realloc(z, nlimbs);
    MPI_Bcast(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, root, comm);
    MPI_Bcast(&z->_mp_size, 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, root, comm);
}/*}}}*/

/* {{{ share_1mpz -- two nodes communicate with eachother. integer z is
 * present on input for both nodes. Node s sends it to node r, which of
 * course receives it in integer t. The integer t is unused by node s. It
 * is an error to call this routine if the current node is not one of
 * {r,s}
 */
static void share_1mpz(mpz_ptr z, mpz_ptr t, int r, int s, MPI_Comm comm)
{
    // r <-- s
    int k;
    MPI_Comm_rank(comm, &k);
    mp_size_t nlimbs = mpz_size(z);

    if (k == r) {
        MPI_Recv(&nlimbs, 1, MPI_MY_MP_SIZE_T, s, (r<<4), comm, MPI_STATUS_IGNORE);
        if (nlimbs == 0) {
            return;
        }

        _mpz_realloc(t, nlimbs);
        MPI_Recv(t->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, s, 1+(r<<4), comm, MPI_STATUS_IGNORE);
        MPI_Recv(&(t->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, s, 2+(r<<4), comm, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&nlimbs, 1, MPI_MY_MP_SIZE_T, r, (r<<4), comm);
        if (nlimbs == 0)
            return;

        MPI_Send(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, r, 1+(r<<4), comm);
        MPI_Send(&(z->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, r, 2+(r<<4), comm);
    }
}/*}}}*/
#if 0
/* {{{ share_2mpz is the more complete version, where t is filled with
 * the peer integer on both nodes */
static void share_2mpz(mpz_ptr z, mpz_ptr t, int r, int s, MPI_Comm comm)
{
    int k;
    MPI_Comm_rank(comm, &k);
    int peer = r^s^k;

    mp_size_t nlimbs = mpz_size(z);
    mp_size_t nlimbs_peer = 0;

    MPI_Sendrecv(&nlimbs, 1, MPI_MY_MP_SIZE_T, peer, k<<4,
            &nlimbs_peer, 1, MPI_MY_MP_SIZE_T, peer, peer<<4,
            comm, MPI_STATUS_IGNORE);             
    if (nlimbs == 0 || nlimbs_peer == 0) {
        return;
    }
    _mpz_realloc(t, nlimbs_peer);
    MPI_Sendrecv(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, peer, 1 + (k<<4),
            t->_mp_d, nlimbs_peer, MPI_MY_MP_LIMB_T, peer, 1 + (peer<<4),
            comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&(z->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, peer, 2 + (k<<4),
            &(t->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, peer, 2 + (peer<<4),
            comm, MPI_STATUS_IGNORE);
}
/* }}} */
#endif
void reduce_mulmod_mpz(mpz_ptr z, int recv, MPI_Comm comm, mpz_srcptr px)/*{{{*/
{
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    ASSERT_ALWAYS(recv == 0);
    mpz_t t;
    mpz_init(t);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        if (me % done_size) continue;
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        share_1mpz(z, t, receiver, sender, comm);
        // could be off-loaded to a working queue.
        if (me == receiver) {
            if (!mpz_size(t) || !mpz_size(z)) {
                /* It's not fundamentally wrong, but in our case, it
                 * indicates a bug, for sure. */
                fprintf(stderr, "Warning: received 0 from job %d\n", s);
                mpz_set_ui(z,0);
            }
            WRAP_mpz_mul(z,z,t);
            WRAP_mpz_mod(z,z, px);
        }
    }
    mpz_clear(t);
}/*}}}*/
void allreduce_mulmod_mpz(mpz_ptr z, MPI_Comm comm, mpz_srcptr px)/*{{{*/
{
#if 0
    /* This works only when the communicator size is a power of two.
     * If it isn't, then some parts don't get the result. Examples:
     *
     * 4 : (1 2 3 4) --> (12 12 34 34) --> (1234 1234 1234 1234).
     * 3 : (1 2 3) --> (12 12 3) --> (123 12 123)
     * 
     * So we would need 3 rounds for a communicator of size 3, while the
     * circular algortihm needs two rounds only.
     */
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    mpz_t t;
    mpz_init(t);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        share_2mpz(z, t, receiver, sender, comm);
        if (!mpz_size(t) || !mpz_size(z)) {
            /* It's not fundamentally wrong, but in our case, it
             * indicates a bug, for sure. */
            fprintf(stderr, "Warning: received 0 from job %d\n", s);
            mpz_set_ui(z,0);
        }
        // could be off-loaded to a working queue.
        WRAP_mpz_mul(z,z,t);
        WRAP_mpz_mod(z,z, px);
    }
    mpz_clear(t);
#else
    /* completely stupid, but works */
    reduce_mulmod_mpz(z, 0, comm, px);
    broadcast_mpz(z, 0, comm);
#endif
}/*}}}*/

void reduce_polymodF_mul_monic(poly_t P, int recv, MPI_Comm comm, poly_t F)/*{{{*/
{
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    ASSERT_ALWAYS(recv == 0);
    poly_t Q;
    poly_alloc(Q, F->deg-1);
    ASSERT_ALWAYS(P->deg < F->deg);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        if (me % done_size) continue;
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        for(int j = 0 ; j < F->deg ; j++) {
            mpz_ptr z = P->coeff[j];
            mpz_ptr t = Q->coeff[j];
            if (j > P->deg)
                mpz_set_ui(z,0);
            share_1mpz(z, t, receiver, sender, comm);
        }
        // could be off-loaded to a working queue.
        if (me == receiver) {
            cleandeg(Q, F->deg - 1);
            polymodF_mul_monic(P, P, Q, F);
        }
    }
    poly_free(Q);
}/*}}}*/
static void broadcast_poly(poly_t P, int maxdeg, int root, MPI_Comm comm) /*{{{*/
{
    /* maxdeg must be <= all allocation degrees. */
    ASSERT_ALWAYS(maxdeg + 1 <= P->alloc);
    for(int j = 0 ; j < maxdeg + 1 ; j++) {
        mpz_ptr z = P->coeff[j];
        if (j > P->deg)
            mpz_set_ui(z,0);
        broadcast_mpz(z, root, comm);
    }
    cleandeg(P, maxdeg);
}/*}}}*/
void allreduce_polymodF_mul_monic(poly_t P, MPI_Comm comm, poly_t F)/*{{{*/
{
    reduce_polymodF_mul_monic(P, 0, comm, F);
    broadcast_poly(P, F->deg - 1, 0, comm);
}
/* }}} */

/* }}} */

// {{{ interface for reading the list of (a,b)'s, with sort of a random
// access (for the starting point only).

#define ABFILE_MAX_LINE_LENGTH  256

struct ab_source_s {
    const char * fname0;
    char * sname;
    size_t sname_len;
    size_t nab;
    char prefix[100];
    int depnum;
    // int cado; // use !numfiles instead.

    /* this relates to rough estimations based on the size(s) of the
     * file(s). Note however that apparently the guess logic isn't too
     * good at taking the leading coefficient into account, so it's not
     * meant to be used before it gets fixed. Doing the accurate
     * evaluation via complex embeddings is ultra-fast anyway */
    size_t nab_estim;
    // size_t digitbytes_estim;

    /* This relates to the different files (if several), their sizes, and
     * the current position. */
    size_t * file_bases;
    int nfiles; // 0 for cado format.
    size_t totalsize;

    FILE * f;
    int c;
    size_t cpos;        // position within current file.
    size_t tpos;        // position within totality.
};

typedef struct ab_source_s ab_source[1];
typedef struct ab_source_s * ab_source_ptr;

void ab_source_init(ab_source_ptr ab, const char * fname, int rank, int root, MPI_Comm comm)
{
    memset(ab, 0, sizeof(ab_source));
    ab->fname0 = fname;
    char * magic;
    if ((magic = strstr(fname, ".prep.")) != NULL) {
        // then assume kleinjung format.
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        int fnum;
        if (sscanf(magic, "prep.%d.rel.%d", &ab->depnum, &fnum) == 2) {
            ab->nfiles = -1;  // to be determined later on.
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((magic = strstr(fname, ".dep.alg.")) != NULL) {
        // assume cado format (means only one file, so we don't need to
        // parse, really.
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        if (sscanf(magic, "dep.alg.%d", &ab->depnum) == 1) {
            ab->nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else {
        FATAL_ERROR_CHECK(1, "error in parsing filename");
    }

    // do some size estimations;
    size_t tsize = 0;
    struct stat sbuf[1];
    int rc;
    if (ab->nfiles == 0) {
        rc = stat(fname, sbuf);
        ASSERT_ALWAYS(mpi_data_agrees(&rc, 1, MPI_INT, comm));
        ASSERT_ALWAYS(rc == 0);
        tsize=sbuf->st_size;
        // we have 2.5 non-digit bytes per file line. However we
        // don't know the line count, so we can't subtract. As a
        // guess, we read the first 16kb, and count the number of
        // lines in there.
        if (rank == root) {
            char buf[16384];
            FILE * f = fopen(fname, "r");
            rc = fread(buf, 1, sizeof(buf), f);
            ASSERT_ALWAYS(rc == sizeof(buf));
            fclose(f);
            int nrows_16k = 0;
            for(unsigned int i = 0 ; i < sizeof(buf) ; i++) {
                nrows_16k += buf[i] == '\n';
            }
            ab->nab_estim = (double) tsize * nrows_16k / sizeof(buf);
            // ab->digitbytes_estim = tsize - 2.5 * ab->nab_estim;
        }
        // XXX OK, this requires endianness consistency.
        MPI_Bcast(&ab->nab_estim, 1, MPI_MY_SIZE_T, root, comm);
        ab->file_bases = malloc(2 * sizeof(size_t));
        ab->file_bases[0] = 0;
        ab->file_bases[1] = tsize;
        ab->totalsize = tsize;
    } else {
        size_t dummy;
        size_t hdrbytes;
        char line[ABFILE_MAX_LINE_LENGTH];
        if (rank == root) {
            FILE * f = fopen(fname, "r");
            char * xx = fgets(line, sizeof(line), f);
            DIE_ERRNO_DIAG(xx == NULL, "fgets", fname);
            rc = sscanf(line, "AB %zu %zu", &ab->nab, &dummy);
            DIE_ERRNO_DIAG(rc != 2, "parse", fname);
            hdrbytes = ftell(f);
            fclose(f);
        }
        MPI_Bcast(&hdrbytes, 1, MPI_MY_SIZE_T, root, comm);

        ab->sname_len = strlen(fname) + 10;
        ab->sname = malloc(ab->sname_len);
        for(ab->nfiles = 0 ; ; ab->nfiles++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, ab->nfiles);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0 || errno == ENOENT);
            ASSERT_ALWAYS(mpi_data_agrees(&rc, 1, MPI_INT, comm));
            if (rc < 0) break;
            tsize += sbuf->st_size - hdrbytes;
        }
        ASSERT_ALWAYS(ab->nfiles > 0);
        ab->nab_estim = ab->nab;
        // ab->digitbytes_estim = tsize - 5 * ab->nab_estim;
        ab->file_bases = malloc((ab->nfiles+1) * sizeof(size_t));
        ab->file_bases[0] = 0;
        for(int i = 0 ; i < ab->nfiles ; i++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, i);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0);
            ASSERT_ALWAYS(mpi_data_agrees(&rc, 1, MPI_INT, comm));
            ab->file_bases[i+1]=ab->file_bases[i] + sbuf->st_size;
        }
        ab->totalsize = ab->file_bases[ab->nfiles];
    }
    /*
       fprintf(stderr, "# [%2.2lf] %s: roughly %zu rows,"
       " %zu digits (%zu bits/c)\n",
       WCT,
       ab->nfiles ? "kleinjung" : "cado",
       ab->nab_estim, ab->digitbytes_estim,
       (size_t) (ab->digitbytes_estim * M_LN10 / M_LN2));
       */
}

void ab_source_rewind(ab_source_ptr ab)
{
    if (ab->f) fclose(ab->f);
    ab->f = NULL;
    ab->c = 0;
    ab->nab = 0;
    ab->cpos = 0;
    ab->tpos = 0;
}

void ab_source_init_set(ab_source_ptr ab, ab_source_ptr ab0)
{
    memcpy(ab, ab0, sizeof(struct ab_source_s));
    ab->sname = malloc(ab->sname_len);
    ab->file_bases = malloc((ab->nfiles+1) * sizeof(size_t));
    memcpy(ab->file_bases, ab0->file_bases, (ab->nfiles+1) * sizeof(size_t));
    ab->f = NULL;
    ab_source_rewind(ab);
}

void ab_source_clear(ab_source_ptr ab)
{
    ab_source_rewind(ab);
    free(ab->sname);
    free(ab->file_bases);
    memset(ab, 0, sizeof(ab_source));
}

int ab_openfile_internal(ab_source_ptr ab)
{
    const char * s;
    if (ab->nfiles == 0) {
        ab->f = fopen(s=ab->fname0, "r");
        ab->tpos = 0;
    } else {
        snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                ab->prefix, ab->depnum, ab->c);
        ab->f = fopen(s=ab->sname, "r");
        if (ab->f == NULL && errno == ENOENT)
            return 0;
        ab->tpos = ab->file_bases[ab->c];
        ab->cpos = 0;
        char header[80];
        char * rp = fgets(header, sizeof(header), ab->f);
        ASSERT_ALWAYS(rp);
    }
    ab->cpos = ftell(ab->f);
    ab->tpos += ab->cpos;
    DIE_ERRNO_DIAG(ab->f == NULL, "fopen", s);
    return 1;
}

int ab_source_next(ab_source_ptr ab, int64_t * a, uint64_t * b)
{
    if (ab->f) {
        int rc;
        char line[ABFILE_MAX_LINE_LENGTH];
        char * xx = fgets(line, sizeof(line), ab->f);
        size_t cpos = ftell(ab->f);
        if (xx) {
            if (ab->nfiles == 0) {
                rc = sscanf(line, "%" SCNd64 " %" SCNu64, a, b);
                ASSERT(rc == 2);
            } else {
                int dummy;
                rc = sscanf(line, "%d %" SCNd64 " %" SCNu64, &dummy, a, b);
                ASSERT(rc == 3);
            }
            ab->tpos += cpos - ab->cpos;
            ab->cpos = cpos;
            ab->nab++;
            return 1;
        }
        fclose(ab->f); ab->f = NULL;
        ab->tpos += cpos - ab->tpos;
        // don't update cpos, it is not defined in this situation.
        if (ab->nfiles == 0)
            return 0;
        ab->c++;
    }
    if (ab_openfile_internal(ab) == 0)
        return 0;
    return ab_source_next(ab, a, b);
}

void ab_source_move_afterpos(ab_source_ptr ab, size_t offset)
{
    /* move the file pointer to the earliest non-header line starting at
     * an offset which is greater than or equal to offset.
     *
     * IF the current file position is already >= offset, do nothing.
     *
     * IF the current file position is < offset, then the returned offset
     * is >= offset. This is achieved by seeking to some pre_offset <
     * offset, and
     * advancing to the next line starting at >= offset.
     *
     * This is used to read a collection of files in chunks whose size is
     * governed by the file size.
     */

    /* note that the test below always succeeds for offset == 0 */
    if (ab->tpos >= offset)
        return;
    // otherwise it's really, really a can of worms.
    ASSERT_ALWAYS(ab->f == NULL);

    // which file ?
    FATAL_ERROR_CHECK(offset >= ab->totalsize,
            "attempt to seek beyond end of files");

    size_t pre_offset = MAX(offset, 10) - 10;
    ASSERT_ALWAYS(pre_offset < offset);

    if (ab->nfiles == 0) {
        // well, does not really make a lot of sense here, but anyway.
        ab_openfile_internal(ab);
        fseek(ab->f, pre_offset, SEEK_SET);
    } else {
        for( ; ab->file_bases[ab->c+1] <= pre_offset ; ab->c++) ;
        ab_openfile_internal(ab);
        // so we know that
        // ab->file_bases[ab->c] <= pre_offset < ab->file_bases[ab->c+1]
        fseek(ab->f, pre_offset - ab->file_bases[ab->c], SEEK_SET);
        // note that it is possible that tpos, as obtained after
        // openfile, is >= offset. but we know that we'll be able to seek
        // to pre_offset, and that one is appropriately < offset, so no
        // special case.
    }
    char line[ABFILE_MAX_LINE_LENGTH];
    char * xx = fgets(line, sizeof(line), ab->f);
    DIE_ERRNO_DIAG(xx == NULL, "fgets", ab->nfiles ? ab->sname : ab->fname0);
    size_t cpos = ftell(ab->f);
    ab->tpos += cpos - ab->cpos;
    ab->cpos = cpos;
    for(int n_adjust = 0 ; ab->tpos < offset ; n_adjust++) {
        FATAL_ERROR_CHECK(n_adjust > 10, "adjustment on the runaway");
        int64_t a;
        uint64_t b;
        int r = ab_source_next(ab, &a, &b);
        FATAL_ERROR_CHECK(r == 0, "adjustment failed");
    }
    logprint("(a,b) rewind to %s, pos %zu\n",
            ab->nfiles ? ab->sname : ab->fname0, ab->cpos);
}
// }}}

/* {{{ POSIX threads stuff : work queues */

typedef void * (*wq_func_t)(void *);

struct wq_task {
    wq_func_t f;
    void * arg;

    // the remaining fields are reserved. In particular, access to done
    // must be mutex protected.
    int done;

    pthread_mutex_t m_[1];
    pthread_cond_t c_[1];

    struct wq_task * next;
    // struct wq_task * prev;
};

struct work_queue {
    pthread_mutex_t m[1];
    pthread_cond_t c[1];   // used by waiters.
    pthread_t * clients;
    // TODO: Use a single tasklist item for a doubly linked list.
    struct wq_task * head;
};

// puts a task on the pending list.
// returns a handle which might be waited for by a join.
struct wq_task * wq_push(struct work_queue * wq, wq_func_t f, void * arg)
{

    struct wq_task * t = malloc(sizeof(struct wq_task));
    t->f = f;
    t->arg = arg;
    t->done = 0;
    if (t->f) {
        // t->f == NULL is used as a special marker. In this case the
        // mutexes are not initialized
        pthread_mutex_init(t->m_, NULL);
        pthread_cond_init(t->c_, NULL);
    }

    pthread_mutex_lock(wq->m);
    t->next = wq->head;
    wq->head = t;
    pthread_cond_signal(wq->c);
    pthread_mutex_unlock(wq->m);
    return t;
}

// grabs a task as soon as one is available.
struct wq_task * wq_pop_wait(struct work_queue * wq)
{
    struct wq_task * t;
    pthread_mutex_lock(wq->m);
    for( ; wq->head == NULL ; ) {
        pthread_cond_wait(wq->c, wq->m);
    }
    t = wq->head;
    wq->head = t->next;
    pthread_mutex_unlock(wq->m);
    t->next = NULL;
    return t;
}

void * wq_waiter(void * wq)
{
    for( ; ; ) {
        struct wq_task * t = wq_pop_wait((struct work_queue *) wq);
        if (t->f == NULL) {
            /* t == NULL is an indication that we're reaching the
             * end-of-work signal. For this special case, t->m_ and t->c_
             * have not been initialized.  */
            free(t);
            break;
        }
        void * res = (*t->f)(t->arg);
        pthread_mutex_lock(t->m_);
        t->done = 1;
        t->arg = res;   /* We reuse the arg field */
        pthread_cond_signal(t->c_);     // signal waiters.
        pthread_mutex_unlock(t->m_);
    }
    /* we have to obey the pthread_create proto. */
    return NULL;
}

void * wq_join(struct wq_task * t)
{
    pthread_mutex_lock(t->m_);
    for( ; t->done == 0 ; ) {
        pthread_cond_wait(t->c_, t->m_);
    }
    void * res = t->arg;
    pthread_mutex_unlock(t->m_);
    pthread_cond_destroy(t->c_);
    pthread_mutex_destroy(t->m_);
    free(t);
    return res;
}

void wq_init(struct work_queue * wq, unsigned int n)
{
    pthread_mutex_init(wq->m, NULL);
    pthread_cond_init(wq->c, NULL);
    wq->head = NULL;
    wq->clients = malloc(n * sizeof(pthread_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_create(wq->clients + i, NULL, &wq_waiter, wq);
    }
}

void wq_clear(struct work_queue * wq, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        // schedule end-of-work for everybody.
        wq_push(wq, NULL, NULL);
    }
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_join(wq->clients[i], NULL);
    }

    ASSERT_ALWAYS(wq->head == NULL);
    pthread_mutex_destroy(wq->m);
    pthread_cond_destroy(wq->c);
}

/* }}} */

// some global variables (sheeh !)

struct sqrt_globals {
    int m;      // nprimes
    int n;      // degree
    int s;      // number of shares of A
    int t;      // number of prime groups
    int r;      // number of primes in each prime group
    int prec;
    mpz_t P;    // prime product (not to the power prec)
    size_t nbits_sqrt;
    // size_t nbits_a;
    mpz_t * f_hat;
    mpz_t * f_hat_diff;
    double f_hat_coeffs;
    cado_poly pol;
    poly_t F;
    ab_source ab;
    int lll_maxdim;
    int rank;
    int nprocs;
    int ncores;
    struct work_queue wq[1];
    MPI_Comm acomm;     // same share of A
    MPI_Comm pcomm;     // same sub-product tree
    int arank, asize;
    int prank, psize;
    barrier_t barrier[1];
};

struct sqrt_globals glob = { .lll_maxdim=50, .ncores = 2 };

// {{{ trivial utility
static const char * size_disp(size_t s, char buf[16])
{
    char * prefixes = "bkMGT";
    double ds = s;
    const char * px = prefixes;
    for( ; px[1] && ds > 500.0 ; ) {
        ds /= 1024.0;
        px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}
// }}}



// {{{ TODO: Now that the v field is gone, replace the polymodF layer.
// Here's the only fragments which need to remain.
    static int
poly_normalized_p (const poly_t f)
{
    return (f->deg == -1) || mpz_cmp_ui (f->coeff[f->deg], 0) != 0;
}

static void
poly_from_ab_monic(poly_t tmp, long a, unsigned long b) {
    tmp->deg = b != 0;
    mpz_set_ui (tmp->coeff[1], b);
    mpz_neg (tmp->coeff[1], tmp->coeff[1]);
    mpz_set_si (tmp->coeff[0], a);
    mpz_mul(tmp->coeff[0], tmp->coeff[0], glob.pol->f[glob.n]);
}

    static void
poly_reducemodF_monic(poly_t P, poly_t p, const poly_t F)
{
    if (p->deg < F->deg) {
        poly_copy(P, p);
        return;
    }
    const int d = F->deg;
    while (p->deg >= d) {
        const int k = p->deg;
        for (int i = 0; i < d; ++i) 
            mpz_submul (p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

        cleandeg (p, k-1);
    }

    poly_copy(P, p);
}

    void
polymodF_mul_monic (poly_ptr Q, poly_srcptr P1, poly_srcptr P2,
        poly_srcptr F)
{
    poly_t prd;
    poly_alloc(prd, P1->deg+P2->deg);
    ASSERT_ALWAYS(poly_normalized_p (P1));
    ASSERT_ALWAYS(poly_normalized_p (P2));
    poly_mul(prd, P1, P2);
    poly_reducemodF_monic(Q, prd, F);
    poly_free(prd);
}

void poly_swap(poly_ptr a, poly_ptr b)
{
    ASSERT_ALWAYS(a->deg + 1 <= b->alloc);
    ASSERT_ALWAYS(b->deg + 1 <= a->alloc);
    for(int i = 0 ; i <= a->deg || i<=b->deg ; i++) {
        mpz_swap(a->coeff[i], b->coeff[i]);
    }
    int d = a->deg;
    a->deg = b->deg;
    b->deg = d;
}

// }}}

// {{{ floating point stuff
// {{{ getting the coefficients of the lagrange interpolation matrix.
long double complex lagrange_polynomial(long double complex * res, double * f, int deg, long double complex r)
{
    long double complex z = f[deg];
    long double complex y = deg * f[deg];
    res[deg-1] = z;
    for(int i = deg-1 ; i > 0 ; i--) {
        z *= r; z += f[i];
        y *= r; y += i * f[i];
        res[i-1] = z;
    }
    for(int i = deg-1 ; i >= 0 ; i--) {
        res[i] /= y;
    }
    return z;
}

double lagrange_polynomial_abs(double * res, double * f, int deg, long double complex r)
{
    long double complex * cres = malloc(deg * sizeof(long double complex));
    long double complex z = lagrange_polynomial(cres, f, deg, r);
    for(int i = deg-1 ; i >= 0 ; i--) {
        res[i] = cabsl(cres[i]);
    }
    free(cres);
    return cabsl(z);
}
// }}}

void estimate_nbits_sqrt(size_t * sbits, ab_source_ptr ab) // , int guess)
{
    size_t abits[1];
    /*
       if (guess) {
     *abits = ab->digitbytes_estim * M_LN10 / M_LN2;
     *sbits = *abits / 2;
    // when doing gross estimates like this, we can hardly avoid
    // taking a safety margin.
     *abits += *abits / 10;
     *sbits += *sbits / 10;
     logprint("coefficients of A"
     " have at most %zu bits (ESTIMATED)\n", WCT, *abits);
     logprint("square root coefficients"
     " have at most %zu bits (ESTIMATED)\n", WCT, *sbits);
     return;
     }
     */

    int n = glob.pol->degree;

    long double complex * eval_points = malloc(n * sizeof(long double complex));
    double * double_coeffs = malloc((n+1) * sizeof(double));
    double * evaluations = malloc(n * sizeof(double));

    // take the roots of f, and multiply later on to obtain the roots of
    // f_hat. Otherwise we encounter precision issues.
    for(int i = 0 ; i <= n ; i++) {
        double_coeffs[i] = mpz_get_d(glob.pol->f[i]);
    }

    int rc = poly_roots_longdouble(double_coeffs, n, eval_points);
    if (rc) {
        logprint("Warning: rootfinder had accuracy problem with %d roots\n", rc);
    }


    // {{{ compress the list of roots.
    int nreal = 0, ncomplex = 0, rs = 0;
    for(int i = 0 ; i < n ; i++) {
        if (cimagl(eval_points[i]) > 0) {
            eval_points[rs] = eval_points[i];
            rs++;
            ncomplex++;
        } else if (cimagl(eval_points[i]) < 0) {
            continue;
        } else {
            eval_points[rs] = creall(eval_points[i]);
            // eval_points[rs] = eval_points[i];
            rs++;
            nreal++;
        }
    }
    // }}}

    // {{{ post-scale to roots of f_hat, and store f_hat instead of f
    for(int i = 0 ; i < rs ; i++) {
        eval_points[i] *= mpz_get_d(glob.pol->f[n]);
    }
    for(int i = 0 ; i <= n ; i++) {
        double_coeffs[i] = mpz_get_d(glob.f_hat[i]);
    }
    // }}}

    // {{{ print the roots.
    if (nreal) {
        fprintf(stderr, "# [%2.2lf]", WCT);
        fprintf(stderr, " real");
        for(int i = 0 ; i < rs ; i++) {
            double r = cimagl(eval_points[i]);
            if (r == 0) {
                fprintf(stderr, " %.4Lg", creall(eval_points[i]));
            }
        }
    }
    if (ncomplex) {
        fprintf(stderr, " complex");
        for(int i = 0 ; i < rs ; i++) {
            if (cimagl(eval_points[i]) > 0) {
                fprintf(stderr, " %.4Lg+i*%.4Lg", creall(eval_points[i]), cimagl(eval_points[i]));
            }
        }
    }
    fprintf(stderr, "\n");
    // }}}

    // {{{ now evaluate the product.
    for(int i = 0 ; i < n ; i++) {
        evaluations[i] = 0;
    }

    // Consider the product A of all a-b\alpha's, which is a
    // polynomial in \alpha. This polynomial has _rational_ coefficients,
    // since each reduction involves dividing out by f_d.
    //
    // Easier to handle is f_d^nab*A (product of all f_d*a-b*f_d*\alpha),
    // which is an element of the order Z[f_d alpha] (f_d
    // alpha is an algebraic integer). It can be expressed with integer
    // coefficients in the powers of f_d\alpha. Its square root, however,
    // is not necessarily in this sub-order of the ring of integers.
    // Therefore we multiply by f_hat'(f_d\alpha), where f_hat is the
    // minimal polynomial of f_d\alpha.
    //
    // We ensure that we've rounded up nab to the next even multiple.

    int64_t a;
    uint64_t b;
    double w1,wt;
    w1 = WCT;
    ab_source_rewind(ab);
    for( ; ab_source_next(ab, &a, &b) ; ) {
        for(int i = 0 ; i < rs ; i++) {
            long double complex y = a * mpz_get_d(glob.pol->f[n]);   
            long double complex w = eval_points[i] * b;
            y = y - w;
            evaluations[i] += log(cabsl(y));
        }
        wt = WCT;
        if (wt > w1 + print_delay || !(ab->nab % 10000000)) {
            w1 = wt;
            fprintf(stderr,
                    "# [%2.2lf] floating point evaluation: %zu (%.1f%%)\n",
                    WCT, ab->nab, 100.0*(double)ab->nab/ab->nab_estim);
        }
    }
    // }}}
    // note that now that we've read everything, we know the precise
    // number of (a,b)'s. Thus we can replace the estimation.
    ab->nab_estim = ab->nab;

    // {{{ post-process evaluation: f'(alpha), and even nab. print.

    if (ab->nab & 1) {
        fprintf(stderr, "# [%2.2lf] odd number of pairs !\n", WCT);
        for(int i = 0 ; i < rs ; i++) {
            evaluations[i] += log(fabs(mpz_get_d(glob.pol->f[n])));
        }
    }

    // multiply by the square of f_hat'(f_d\alpha).
    for(int i = 0 ; i < rs ; i++) {
        long complex double s = n;
        for(int j = n - 1 ; j >= 0 ; j--) {
            s *= eval_points[i];
            s += double_coeffs[j] * j;
        }
        evaluations[i] += 2 * clog(s);
    }
    fprintf(stderr, "# [%2.2lf] Log_2(A)", WCT);
    for(int i = 0 ; i < rs ; i++) {
        fprintf(stderr, " %.4Lg", creall(evaluations[i]) / M_LN2);
        if (cimagl(eval_points[i]) > 0) {
            fprintf(stderr, "*2");
        }
    }
    fprintf(stderr, "\n");
    // }}}

    // {{{ deduce the lognorm. print.
    double lognorm = 0;
    for(int i = 0 ; i < rs ; i++) {
        if (cimagl(eval_points[i]) > 0) { 
            lognorm += 2*creall(evaluations[i]);
        } else {
            lognorm += creall(evaluations[i]);
        }
    }
    fprintf(stderr, "# [%2.2lf] log_2(norm(A)) %.4g\n",
            WCT, lognorm / M_LN2);
    // }}}

    // {{{ now multiply this by the Lagrange matrix.
    double * a_bounds = malloc(n * sizeof(double));
    double * sqrt_bounds = malloc(n * sizeof(double));

    for(int j = 0 ; j < n ; j++) {
        a_bounds[j] = 0;
        sqrt_bounds[j] = 0;
    }
    for(int i = 0 ; i < rs ; i++) {
        double * lmat = malloc(n * sizeof(double ));
        lagrange_polynomial_abs(lmat, double_coeffs, n, eval_points[i]);
        for(int j = 0 ; j < n ; j++) {
            double za, zs;
            za = zs = log(lmat[j]);
            za += evaluations[i];
            zs += evaluations[i] / 2;
            if (cimagl(eval_points[i]) > 0) {
                za += log(2);
                zs += log(2);
            }
            if (za > a_bounds[j]) a_bounds[j] = za;
            if (zs > sqrt_bounds[j]) sqrt_bounds[j] = zs;
        }
        free(lmat);
    }
    // }}}

    // {{{ get global logbounds
    double logbound_sqrt = 0;
    double logbound_a = 0;
    for(int j = 0 ; j < n ; j++) {
        // note that we might have added up to n times the same thing.
        // (inequality a+b < 2max(a,b) )
        sqrt_bounds[j] += log(n);
        a_bounds[j] += log(n);

        // safety margin for inaccuracies ?
        sqrt_bounds[j] += 100 * M_LN2;
        a_bounds[j] += 100 * M_LN2;
        if (sqrt_bounds[j] > logbound_sqrt) logbound_sqrt = sqrt_bounds[j];
        if (a_bounds[j] > logbound_a) logbound_a = a_bounds[j];
    }
    // }}}

#if 0
    // {{{ print logbounds
    fprintf(stderr, "# [%2.2lf] logbounds per coeff of A", WCT);
    for(int j = 0 ; j < n ; j++) {
        fprintf(stderr, " %.4Lg", a_bounds[j] / M_LN2);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "# [%2.2lf] logbounds per coeff of sqrt(A)", WCT);
    for(int j = 0 ; j < n ; j++) {
        fprintf(stderr, " %.4Lg", sqrt_bounds[j] / M_LN2);
    }
    fprintf(stderr, "\n");
    // }}}
#endif

    // Since we're being very gross, we use a stupid bound on the last
    // coefficient.

    *sbits = ceil(logbound_sqrt / M_LN2);
    *abits = ceil(logbound_a / M_LN2);
    fprintf(stderr, "# [%2.2lf] coefficients of A"
            " have at most %zu bits\n", WCT, *abits);
    fprintf(stderr, "# [%2.2lf] square root coefficients"
            " have at most %zu bits\n", WCT, *sbits);

    free(a_bounds);
    free(sqrt_bounds);
    free(eval_points);
    free(double_coeffs);
    free(evaluations);
}
/* }}} */

/*{{{ have to pre-declare prime_data */
struct alg_ptree_s;

#if 0
struct individual_contribution {
    uint64_t ratio;
    mpz_t modN;
};
#endif

struct prime_data {
    unsigned long p;
    unsigned long * r;
    // unsigned long rj;
    void * powers;      // see .cpp file.

    // computed somewhat late.
    mpz_t iHx;

    struct alg_ptree_s * T;     // always NULL.

    // after the rational ptree reduction, this contains the share of A
    // with coefficients reduced. Short-lived.
    poly_t A;

    poly_t evals;       // set of evaluations.
    poly_t lroots;      // (lifted) roots.
    poly_t sqrts;       // square roots of A(x)     (only in the end !)

};/* }}} */

/* {{{ product tree (algebraic) */
struct alg_ptree_s {
    struct alg_ptree_s * t0;
    struct alg_ptree_s * t1;
    poly_t s;
};
typedef struct alg_ptree_s alg_ptree_t;

#if 0
alg_ptree_t * alg_ptree_build(struct prime_data * p, int i0, int i1)
{
    /* Everything being done the naive way, the algebraic ptree wins
     * nothing: the count of multiplications is exactly the same (n^2-n),
     * and the algebraic ptree incurs an extra storage for n
     * coefficients. It's thus disabled.
     */
    return NULL;

    mpz_srcptr px = power_lookup(p->powers, glob.prec);
    ASSERT_ALWAYS(i0 < i1);
    alg_ptree_t * res = malloc(sizeof(alg_ptree_t));
    memset(res, 0, sizeof(alg_ptree_t));
    poly_alloc(res->s, i1-i0);
    if (i1-i0 == 1) {
        res->s->deg = 1;
        mpz_set_ui(res->s->coeff[1], 1);
        mpz_neg(res->s->coeff[0], p->rx);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = alg_ptree_build(p, i0, i0+d);
    res->t1 = alg_ptree_build(p, i0+d, i1);
    poly_mul(res->s, res->t0->s, res->t1->s);
    poly_reducemodF_monic(res->s, res->s, glob.F);
    poly_reduce_mod_mpz (res->s, res->s, px);
    return res;
}
#endif

void alg_ptree_clear(alg_ptree_t * T)
{
    if (T == NULL) return;
    alg_ptree_clear(T->t0);
    alg_ptree_clear(T->t1);
    free(T);
}

/* }}} */

/* {{{ finding the CRT primes */

typedef int (*sortfunc_t) (const void *, const void *);
int modul_cmp_sortfunc(const unsigned long * a, const unsigned long * b)
{
    return (*a>*b) - (*b>*a);
}

struct prime_data * suitable_crt_primes()
{
    unsigned int m = glob.m;
    struct prime_data * res = malloc(m * sizeof(struct prime_data));
    unsigned int i = 0;

    // Note that p0 has to be well above the factor base bound. Thus on
    // 32-bit machines, it's possibly problematic.
    unsigned long p0 = ((~0UL)>>1)+1;    // 2^31 or 2^63
    unsigned long p = p0;
    unsigned long * roots = malloc(glob.n * sizeof(unsigned long));

    if (glob.rank == 0)
        fprintf(stderr, "# [%2.2lf] Searching for CRT primes\n", WCT);
    // fprintf(stderr, "# [%2.2lf] p0=%lu\n", WCT, p0);

    for( ; i < m ; ) {
        p = ulong_nextprime(p);
        modulusul_t q;
        modul_initmod_ul(q, p);
        memset(roots, 0, glob.n * sizeof(unsigned long));
        int nr = modul_poly_roots_ulong(roots, glob.pol->f, glob.n, q);
        if (nr != glob.n) continue;
        memset(&(res[i]), 0, sizeof(struct prime_data));
        res[i].r = malloc(glob.n * sizeof(unsigned long));
        memcpy(res[i].r, roots, glob.n * sizeof(unsigned long));
        res[i].p = p;
        // res[i].log2_p = log(p)/M_LN2;
        // fprintf(stderr, "# [%2.2lf] p0+%lu", WCT, p-p0);
        // if (glob.rank == 0) fprintf(stderr, "#[%2.2lf] %lu\n", WCT, p);
        residueul_t fd;
        modul_init(fd, q);
        modul_set_ul_reduced(fd, mpz_fdiv_ui(glob.pol->f[glob.n], modul_getmod_ul (q)), q);
        for(int j = 0 ; j < glob.n ; j++) {
            residueul_t r;
            modul_init(r, q);
            modul_set_ul_reduced(r, res[i].r[j], q);
            modul_mul(r,r,fd,q);
            res[i].r[j] = modul_get_ul(r, q);
            // fprintf(stderr, " %lu", res[i].r[j]);
            modul_clear(r, q);
        }
        // sort with respect to the values of the scaled roots.
        qsort(res[i].r, glob.n, sizeof(unsigned long), (sortfunc_t) &modul_cmp_sortfunc);
        modul_clear(fd, q);
        modul_clearmod(q);
        // fprintf(stderr, "\n");
        i++;
    }
    if (glob.rank == 0)
        fprintf(stderr, "# [%2.2lf] Found all CRT primes\n", WCT);
    free(roots);

    return res;
}
/* }}} */

/* {{{ everything that happens only modulo one prime */

/* {{{ tonelli-shanks */
void modul_find_ts_gen(residueul_t z, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = ctzl(pp);
    pp >>= e;
    unsigned long s = 1UL << (e-1);
    residueul_t r;
    modul_init(r, p);
    do {
        modul_set_ul(z, random(), p);
        modul_pow_ul(z, z, pp, p);
        modul_pow_ul(r, z, s, p);
        modul_add_ul(r, r, 1, p);
    } while (!modul_is0(r, p));
    modul_clear(r, p);
}

int modul_field_sqrt(residueul_t z, residueul_t a, residueul_t g, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = ctzl(pp);
    pp >>= e+1;
    if (modul_is0(a, p)) {
        modul_set0(z, p);
        return 1;
    }
    if (modul_is0(g, p)) {
        modul_find_ts_gen(g, p);
    }
    residueul_t b, x, y, t;
    modul_init(b, p);
    modul_init(x, p);
    modul_init(y, p);
    modul_init(t, p);
    int r = e;
    unsigned long s = 1UL << (e-1);

    // modul_set(x, a, p);
    modul_set(y, g, p);

    modul_pow_ul(x, a, pp, p);
    modul_sqr(b, x, p);
    modul_mul(x, x, a, p);
    modul_mul(b, b, a, p);

    int m;
    for(;;) {
        modul_set(t, b, p);
        for(m=0; !modul_is1(t, p); m++)
            modul_sqr(t, t, p);
        assert(m<=r);

        if (m==0 || m==r)
            break;

        s = 1UL << (r-m-1);
        r = m;

        modul_pow_ul(t, y, s, p);
        modul_sqr(y, t, p);
        modul_mul(x, x, t, p);
        modul_mul(b, b, y, p);
    }

    modul_set(z, x, p);
    modul_clear(t, p);
    modul_clear(x, p);
    modul_clear(y, p);
    modul_clear(b, p);
    return (m==0);
}
/* }}} */
/* {{{ sqrt / invsqrt */
void invsqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1) {
        residueul_t z, a;
        modulusul_t q;
        modul_initmod_ul(q, p->p);
        modul_init(z, q);
        modul_init(a, q);
        modul_set_ul_reduced(a, mpz_get_ui(A), q);
        int issquare = modul_field_sqrt(a, a, z, q);
        ASSERT_ALWAYS(issquare);
        modul_inv(a, a, q);
        if (modul_get_ul(a, q) & 1) modul_neg(a, a, q);
        mpz_set_ui(sx, modul_get_ul(a, q));
        modul_clear(a, q);
        modul_clear(z, q);
        modul_clearmod(q);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup_const(p->powers, lower);

    // we're going to recurse ; arrange so that we recurse on something
    // acceptably small.  The problem with this approach is that we store
    // many reductions in memory.
    mpz_t A_save;
    mpz_init(A_save);
    WRAP_mpz_mod(A_save, A, pl);
    // recurse.
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_t ta;
    mpz_init(ta);

    WRAP_mpz_mul(ta, sx, sx);
    WRAP_mpz_mod(ta, ta, pk);
    WRAP_mpz_mul(ta, ta, A);
    WRAP_mpz_mod(ta, ta, pk);
    mpz_sub_ui(ta, ta, 1);
    if (mpz_odd_p(ta)) mpz_add(ta, ta, pk);
    mpz_div_2exp(ta, ta, 1);
    mpz_submul(sx, ta, sx);
    WRAP_mpz_mod(sx, sx, pk);

    mpz_clear(ta);
}

void sqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
{
    double w0 = WCT;
    int lower = precision - precision / 2;
    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    mpz_srcptr pl = power_lookup_const(p->powers, lower);
    mpz_t A_save;
    mpz_init(A_save);
    WRAP_mpz_mod(A_save, A, pl);
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);
    // inverse square root now in sx.

    mpz_t tmp;
    mpz_init(tmp);

    WRAP_mpz_mul(tmp, A, sx);
    WRAP_mpz_mod(tmp, tmp, pl);
    // XXX This destroys A !!!
    mpz_submul(A, tmp, tmp);
    if (mpz_odd_p(sx)) mpz_add(sx, sx, pk);
    mpz_div_2exp(sx, sx, 1);
    mpz_addmul(tmp, sx, A);
    WRAP_mpz_mod(sx, tmp, pk);

    mpz_clear(tmp);
}
/* }}} */

/* {{{ evaluation of polynomials (one or two polynomials) */
static void mp_poly_eval_mod(mpz_ptr r, mpz_t * poly, int deg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
        WRAP_mpz_mul(r, r, a);
        // falls back on mpz_mod if qx == NULL
        WRAP_barrett_mod(r, r, q, qx);
        // mpz_mod(r, r, q);
        mpz_add(r, r, poly[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
}

static void mp_2poly_eval_mod(mpz_ptr r, mpz_ptr s, mpz_t * f, mpz_t * g, int degf, int degg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    if (r == NULL) {
        mp_poly_eval_mod(s, g, degg, a, q, qx);
        return;
    }

    if (s == NULL) {
        mp_poly_eval_mod(r, f, degf, a, q, qx);
        return;
    }

    mpz_t w;
    mpz_init(w);
    mpz_set(r,f[0]);
    mpz_set(s,g[0]);
    mpz_set(w, a);
    WRAP_mpz_addmul(r, w, f[1]);
    WRAP_mpz_addmul(s, w, g[1]);
    for(i = 2 ; i <= degf && i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    for( ; i <= degf ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
    }
    for( ; i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
    WRAP_barrett_mod(s, s, q, qx);
    mpz_clear(w);
}
/* }}} */

/* {{{ This wraps around barrett reduction, depending on whether we
 * choose to use it or not */
mpz_ptr my_barrett_init(mpz_srcptr px MAYBE_UNUSED)
{
#ifndef WITH_BARRETT
    return NULL;
#else
    if (mpz_size(px) < 10000)
        return NULL;
    mpz_ptr qx = malloc(sizeof(mpz_t));
    mpz_init(qx);
    barrett_init(qx, px);
    return qx;
#endif
}

void my_barrett_clear(mpz_ptr qx)
{
    if (!qx) return;
    mpz_clear(qx);
    free(qx);
}
/* }}} */

void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/* {{{ */
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1) {
        return;
    }
    int lower = precision - precision / 2;

    // recurse.
    root_lift(p, rx, irx, lower);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);
    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    mpz_srcptr pl = power_lookup_const(p->powers, lower);
    mpz_ptr qk = my_barrett_init(pk);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    mpz_t fprime;
    mpz_init(fprime);

    // we know r to half-precision, 1/f'(r) to quarter-precision. Compute
    // f(r) (to full precision), and obtain 1/f'(r) to half-precision. We
    // compute f'(r) to half-precision as a side-effect of the
    // computation of f(r).
    mp_2poly_eval_mod(tb, fprime, glob.f_hat, glob.f_hat_diff, glob.n, glob.n-1, rx, pk, qk);
    /* use irx. only one iteration of newton.  */
    mpz_srcptr ql = NULL;       /* FIXME if we use barrett reduction */
    WRAP_barrett_mod(fprime, fprime, pl, ql);
    WRAP_mpz_mul(ta, irx, fprime);
    WRAP_barrett_mod(ta, ta, pl, ql);
    mpz_sub_ui(ta,ta,1);
    WRAP_mpz_submul(irx, irx, ta);
    WRAP_barrett_mod(irx, irx, pl, ql);

    mpz_clear(fprime);

    WRAP_mpz_mul(tb, irx, tb);
    mpz_sub(rx, rx, tb);
    WRAP_barrett_mod(rx, rx, pk, qk);

    mpz_clear(ta);
    mpz_clear(tb);

    my_barrett_clear(qk);
}/* }}} */

void prime_initialization(struct prime_data * p)/* {{{ */
{
    // trigger computation of p^glob.prec
    p->powers = power_lookup_table_init(p->p);
    mpz_init(p->iHx);

    poly_alloc(p->A, glob.n-1);
    poly_alloc(p->evals, glob.n-1);
    poly_alloc(p->lroots, glob.n-1);
    poly_alloc(p->sqrts, glob.n-1);
}
/* }}} */

void prime_cleanup(struct prime_data * p)/* {{{ */
{
    power_lookup_table_clear(p->powers);
    mpz_clear(p->iHx);
    free(p->r);
    alg_ptree_clear(p->T);

    poly_free(p->A);
    poly_free(p->evals);
    poly_free(p->lroots);
    poly_free(p->sqrts);

    memset(p, 0, sizeof(struct prime_data));
}
/* }}} */

/* }}} */

/* {{{ knapsack stuff */
struct crtalgsqrt_knapsack {
    knapsack_object ks;
    mpz_t Px;
    mpz_t lcx;
    mpz_t fhdiff_modN;
    mpz_t sqrt_modN;
    const mp_limb_t * tabN;
};

void crtalgsqrt_knapsack_init(struct crtalgsqrt_knapsack * cks)
{
    knapsack_object_init(cks->ks);
    mpz_init(cks->Px);
    mpz_init(cks->lcx);
    mpz_ptr z = cks->fhdiff_modN;
    mpz_init(z);
    mpz_init(cks->sqrt_modN);
}

int crtalgsqrt_knapsack_callback(struct crtalgsqrt_knapsack * cks,
        unsigned long v, int64_t x)
{
    knapsack_object_srcptr ks = cks->ks;
    const int64_t * c64 = cks->ks->tab;
    const mp_limb_t * cN = cks->tabN;
    unsigned int s;
    char * signs = malloc(ks->nelems+1);

    memset(signs, 0, ks->nelems+1);
    for(s = 0 ; s < ks->nelems ; s++) {
        signs[s]=(v & (1UL << s)) ? '+' : '-';
    }

    // now try to look at the solution more closely */
    mpz_t e;
    mpz_init_set_ui(e, 0);
    int spurious = 0;
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        mpz_t z, t, zN, qz, rz;
        mpz_init_set_ui(z, cks->ks->bound);
        mpz_init_set_ui(zN, 0);
        mpz_init_set_ui(t, 0);
        mpz_init(qz);
        mpz_init(rz);
        int64_t sk = ks->bound;
        for(int j = 0 ; j < glob.m * glob.n ; j++) {
            mpz_set_uint64(t, c64[j * glob.n + k]);
            mpz_t w;
            mp_size_t sN = mpz_size(glob.pol->n);
            MPZ_INIT_SET_MPN(w, cN + sN * (j * glob.n + k), sN);
            if (v & (1UL << j)) {
                sk+=c64[j * glob.n + k];
                mpz_add(z, z, t);
                mpz_add(zN, zN, w);
            } else {
                sk-=c64[j * glob.n + k];
                mpz_sub(z, z, t);
                mpz_sub(zN, zN, w);
            }
            mpz_clear(w);
        }
        if (sk >= 2 * ks->bound || sk < 0) {
            fprintf(stderr, "[%s] recombination of coeff in X^%d yields noise (%"PRIu64")\n",signs, k, sk);
            spurious++;
            // break;
        }
        // so we have (sum of r's) mod p^k = something * p^k + small
        // the ``something'' is in the quotient of the division.
        // FIXME
        // The problem is that the ``small'' thing might be small
        // enough that we won't be able to tell apart s and p-s.
        // This implies that we will have to try 2^degree
        // combinations: those with the quotients as given, and
        // the other combinations with one or several quotients
        // lowered by one unit.
        mpz_set(qz,z);
        if (mpz_cmp_ui(qz,0) >= 0) {
            mpz_fdiv_q_2exp(qz, qz, 63);
            mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
            mpz_fdiv_q_2exp(qz, qz, 1);
        } else {
            mpz_neg(qz,qz);
            mpz_fdiv_q_2exp(qz, qz, 63);
            mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
            mpz_fdiv_q_2exp(qz, qz, 1);
            mpz_neg(qz,qz);
        }
        mpz_mul_2exp(rz, qz, 64);
        mpz_sub(rz, z, rz);
        // gmp_printf("[%d] %Zd = %Zd * 2^64 + %Zd\n", k, z, qz, rz);
        mpz_submul(zN, qz, cks->Px);
        mpz_mod(zN, zN, glob.pol->n);
        // gmp_printf("[X^%d] %Zd\n", k, zN);
        // good. we have the coefficient !
        mpz_mul(e, e, glob.pol->m);
        mpz_mul(e, e, glob.pol->f[glob.n]);
        mpz_add(e, e, zN);
        mpz_mod(e, e, glob.pol->n);
        mpz_clear(z);
        mpz_clear(t);
        mpz_clear(zN);
        mpz_clear(qz);
        mpz_clear(rz);

    }

    mpz_mul(e,e,cks->lcx);
    mpz_mod(e,e,glob.pol->n);

    mpz_mul(e,e,cks->fhdiff_modN);
    mpz_mod(e,e,glob.pol->n);

    mpz_set(cks->sqrt_modN, e);

    mpz_clear(e);

    if (spurious) {
        fprintf(stderr, "# [%2.2lf] %lx (%s) %"PRId64" [SPURIOUS]\n",
                WCT, v, signs, (int64_t) x);
    } else {
        gmp_fprintf(stderr, "# [%2.2lf] %lx (%s) %"PRId64" [%Zd]\n",
                WCT, v, signs, (int64_t) x, cks->sqrt_modN);
    }
    free(signs);

    return spurious == 0;
}

/* Compute all the stuff which will be needed for checking the solutions
*/
void crtalgsqrt_knapsack_prepare(struct crtalgsqrt_knapsack * cks, size_t lc_exp)
{
#if 0 /* {{{ display the knapsack contents */
    {
        for(int i = 0 ; i < glob.m ; i++) {
            for(int j =  0 ; j < glob.n ; j++) {
                for(int k = 0 ; k < glob.n ; k++) {
                    int64_t c = contribs64[ ((i * glob.n) + j) * glob.n + k];
                    printf(" %"PRId64, c);
                }
                for(int s = 0 ; s < glob.m * glob.n ; s++) {
                    printf(" %d", s == (i * glob.n + j));
                }
                printf("\n");
            }
        }
        mpz_t z;
        mpz_init(z);
        mpz_ui_pow_ui(z, 2, 64);
        for(int k = 0 ; k < glob.n ; k++) {
            for(int s = 0 ; s < glob.n + glob.m * glob.n ; s++) {
                if (s == k) {
                    gmp_printf(" %Zd", z);
                } else {
                    printf(" 0");
                }
            }
            printf("\n");
        }
        mpz_clear(z);
    }
#endif/*}}}*/

    mpz_powm_ui(cks->Px, glob.P, glob.prec, glob.pol->n);

    mpz_set(cks->lcx, glob.pol->f[glob.n]);
    mpz_powm_ui(cks->lcx, glob.pol->f[glob.n], lc_exp, glob.pol->n);
    mpz_invert(cks->lcx, cks->lcx, glob.pol->n);

    mpz_ptr z = cks->fhdiff_modN;
    // evaluate the derivative of f_hat in alpha_hat mod N, that is lc*m.
    mpz_set_ui(z, 0);
    for(int k = glob.n ; k > 0 ; k--) {
        mpz_mul(z, z, glob.pol->m);
        mpz_mul(z, z, glob.pol->f[glob.n]);
        mpz_addmul_ui(z, glob.f_hat[k], k);
        mpz_mod(z, z, glob.pol->n);
    }
    mpz_invert(z, z, glob.pol->n);

    cks->ks->cb_arg = cks;
    cks->ks->cb = (knapsack_object_callback_t) crtalgsqrt_knapsack_callback;

    cks->ks->bound = 1+ceil(log(glob.m*glob.n)/log(2));

    cks->ks->stride = glob.n;

    unsigned int nelems = cks->ks->nelems = glob.m * glob.n;
    unsigned int k1 = nelems / 2;
    unsigned int k2 = nelems - k1;
    uint64_t n1 = UINT64_C(1) << k1;
    uint64_t n2 = UINT64_C(1) << k2;
    char buf[16];

    fprintf(stderr,
            "# [%2.2lf] Recombination: dimension %u = %u + %u, %s needed\n",
            WCT,
            nelems, k1, k2,
            size_disp((n1+n2) * sizeof(uint64_t), buf)
           );

}

void crtalgsqrt_knapsack_clear(struct crtalgsqrt_knapsack * cks)
{
    mpz_clear(cks->sqrt_modN);
    mpz_clear(cks->Px);
    mpz_clear(cks->lcx);
    mpz_clear(cks->fhdiff_modN);
    knapsack_object_clear(cks->ks);
}
/* }}} */

void mpi_custom_types_check()/*{{{*/
{
    int ret;
    /* If this ever fails, we have to do something to the defines above
     * (and maybe touch the cmake files) */
    MPI_Type_size(MPI_MY_SIZE_T,&ret);   ASSERT_ALWAYS(ret==sizeof(size_t));
    MPI_Type_size(MPI_MY_MP_SIZE_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_size_t));
    MPI_Type_size(MPI_MY_MP_LIMB_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_limb_t));
    MPI_Type_size(MPI_MY_UINT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(uint64_t));
    MPI_Type_size(MPI_MY_INT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(int64_t));
    {
        mpz_t z;
        MPI_Type_size(MPI_MY_GMP_INTERNAL_SIZE_FIELD_T,&ret);
        ASSERT_ALWAYS(ret==sizeof(z->_mp_size));
    }
}/*}}}*/

// {{{ parallelism parameters selection
// how many parts of A shall we consider ?
// now how do we choose the appropriate number of primes ? Let:
//
// R: the number of bits of the ram guide amount given.
// n: the degree
// T: the number of bits of one coefficient of the square root.
//
// Note thus that A occupies a memory space of 2Tn bits.
//
// and
//
// s: the number of parts of A considered
// B: the number of bits of lifted primes p^lambda
// r: the number of primes stored together in a product tree.
// t: the number of product trees considered.
//
// We have:
//
// [1] sR >= 2Tn (parts of A must form the totality of A)
// [2] nrB <= R  (sets of ``reduced'' A's must fit in RAM)
// [3] rtB >= T  (lifted primes must be large enough for the reconstruction)
//
// A further condition is given by the maximum knapsack lattice
// dimension which can be efficiently handled. This is somewhat
// arbitrary, let's say:
//
// [4] rtn <= 80
//
// Note that [1] readily gives s, and that [2+3] give t. Then per
// [2,3], the product rB is placed within an interval [T/t, R/n].
// Obeying [3], either by choosing the smallest/largest admissible
// power of two for r, or simply the smallest/largest admissible
// value, we deduce r and B. (large r will yield more
// parallelism-multithreading, but OTOH small r will yield less
// time).

void get_parameters(int * pr, int * ps, int * pt, int asked_r)
{
    int n=0,r=0,s=0,t=0;
    size_t R = ram_gb * (1UL << 30) * 8UL;
    n = glob.n;
    size_t T = glob.nbits_sqrt + 80;    // some margin for knapsack.
    double ramfraction = 2*T*n / (double) R;
    s = ceil(ramfraction);
    t = ceil((double)n*T / R);
    int dmax = glob.lll_maxdim - n;
    double rmax = (double) dmax / n / t;
    r = rmax;
    if (r == 0) {
        r = 1;
        fprintf(stderr, "Warning: resorting to choosing r=1, which will give a reconstruction effort of dimension %d\n",
                n*t);
    }
    // for(int mask = ~0 ; r & mask ; mask <<= 1) r &= mask;
    // try r=1 always.
    if (asked_r > r) {
        fprintf(stderr, "value %d requested for r is too large, using %d instead\n", asked_r, r);
    } else if (asked_r) {
        r = asked_r;
    } else {
        /* Do we provide a default ? */
        r = 1;
    }
    size_t B = ceil((double)T/t/r);

    char buf[16];
    char buf2[16];
    fprintf(stderr, "# [%2.2lf] A is %s, %.1f%% of %s."
            " Splitting in %d parts\n",
            WCT,
            size_disp(2*T*n/8.0, buf),
            100.0 * ramfraction, size_disp(ram_gb*1024.0*1048576.0, buf2), s);

    fprintf(stderr, "# [%2.2lf] %d groups of %d prime powers,"
            " each of %zu bits\n",
            WCT, t, r, B);

    /* If these asserts fail, it means that the constraints are
     * impossible to satisfy */
    ASSERT_ALWAYS((double)T/t < (double) R/n);
    ASSERT_ALWAYS(B*n <= R);

    fprintf(stderr, "# [%2.2lf] number of pairs is %zu\n", WCT,
            glob.ab->nab);
    *pr = r;
    *ps = s;
    *pt = t;
}
/* }}} */

int rat_red_caches_ok(struct prime_data * primes, int i0, int i1, size_t off0, size_t off1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        nok += cachefile_exists("a_%zu_%zu_mod_%lu_%lu", off0, off1, primes[i].p, glob.prec);
    }
    MPI_Allreduce(MPI_IN_PLACE, &nok, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.m * glob.s;
}/*}}}*/

int alg_red_caches_ok(struct prime_data * primes, int i0, int i1, size_t off0, size_t off1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        struct prime_data * p = &primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("a_%zu_%zu_mod_%lu_%lu_%lu",
                    off0, off1, p->p, p->r[j], glob.prec);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nok, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

int a_shared_caches_ok(struct prime_data * primes, int i0, int i1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        struct prime_data * p = &primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("a_mod_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nok, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

int sqrt_caches_ok(struct prime_data * primes, int i0, int i1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        struct prime_data * p = &primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("sqrt_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nok, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

struct subtask_info_t {
    struct prime_data * p;
    int i0, i1;
    int j;
    size_t off0, off1;
    poly_ptr P;
    poly_ptr P0;
    poly_ptr P1;
    size_t nab_loc;
    struct wq_task * handle;
};

/*{{{ a_poly_read_share and companion */ 
#define ABPOLY_OFFSET_THRESHOLD        16
// NOTE: This does not depend on p (nor r of course).
// the accumulation is done for all data between:
// the first data line starting at offset >= off0 (inclusive)
// the first data line starting at offset >= off1 (exclusive)
size_t accumulate_ab_poly(poly_ptr P, ab_source_ptr ab, size_t off0, size_t off1, poly_t tmp)
{
    size_t res = 0;
    mpz_set_ui(P->coeff[0], 1);
    P->deg = 0;
    if (off1 - off0 < ABPOLY_OFFSET_THRESHOLD) {
        ab_source_move_afterpos(ab, off0);
        for( ; ab->tpos < off1 ; res++) {
            int64_t a;
            uint64_t b;
            int r = ab_source_next(ab, &a, &b);
            FATAL_ERROR_CHECK(!r, "dep file ended prematurely\n");
            poly_from_ab_monic(tmp, a, b);
            polymodF_mul_monic(P, P, tmp, glob.F);
        }
        return res;
    }
    size_t d = (off1 - off0) / 2;
    poly_t Pl, Pr;
    poly_alloc(Pl, glob.n);
    poly_alloc(Pr, glob.n);
    res += accumulate_ab_poly(Pl, ab, off0, off0 + d, tmp);
    res += accumulate_ab_poly(Pr, ab, off0 + d, off1, tmp);
    polymodF_mul_monic(P, Pl, Pr, glob.F);
    poly_free(Pl);
    poly_free(Pr);
    return res;
}


void * a_poly_read_share_child(struct subtask_info_t * info)
{
    size_t off0 = info->off0;
    size_t off1 = info->off1;
        int rc;
    poly_ptr P = info->P;
    cachefile c;

    cachefile_init(c, "a_%zu_%zu", off0, off1);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        for(int i = 0 ; i < glob.n ; i++) {
            rc = gmp_fscanf(c->f, "%Zx", info->P->coeff[i]);
            ASSERT_ALWAYS(rc == 1);
        }
        cleandeg(info->P, glob.n - 1);
        cachefile_close(c);
        return NULL;
    }

    ab_source ab;
    ab_source_init_set(ab, glob.ab);
    poly_t tmp;
    poly_alloc(tmp, 1);
    info->nab_loc = accumulate_ab_poly(P, ab, off0, off1, tmp);
    poly_free(tmp);
    ab_source_clear(ab);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        for(int i = 0 ; i < glob.n ; i++) {
            gmp_fprintf(c->f, "%Zx\n", info->P->coeff[i]);
        }
        cachefile_close(c);
    }

    return NULL;
}

void * a_poly_read_share_child2(struct subtask_info_t * info)
{
    polymodF_mul_monic(info->P0, info->P0, info->P1, glob.F);
    return NULL;
}

size_t a_poly_read_share(poly_t P, size_t off0, size_t off1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    size_t nab_loc = 0;
    int rc;
    log_begin();

    cachefile c;

    cachefile_init(c, "a_%zu_%zu", off0, off1);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &nab_loc);
        ASSERT_ALWAYS(rc == 1);
        for(int i = 0 ; i < glob.n ; i++) {
            rc = gmp_fscanf(c->f, "%Zx", P->coeff[i]);
            ASSERT_ALWAYS(rc == 1);
        }
        cleandeg(P, glob.n - 1);
        cachefile_close(c);
        return nab_loc;
    }

    size_t nparts = glob.ncores * glob.t;

    struct subtask_info_t * a_tasks;
    a_tasks = malloc(glob.ncores *sizeof(struct subtask_info_t));

    int j0 = glob.ncores * glob.arank;
    int j1 = j0 + glob.ncores;

    poly_t * pols = malloc((j1-j0) * sizeof(poly_t));
    for(int j = j0 ; j < j1 ; j++) poly_alloc(pols[j-j0], glob.n);

    for(int j = j0 ; j < j1 ; j++) {
        struct subtask_info_t * task = a_tasks + (j-j0);
        task->P = pols[j-j0];
        task->off0 = off0 + (off1 - off0) * j / nparts;
        task->off1 = off0 + (off1 - off0) * (j+1) / nparts;
        task->nab_loc = 0;
        wq_func_t f = (wq_func_t) &a_poly_read_share_child;
        task->handle = wq_push(glob.wq, f, task);
    }
    /* we're doing nothing, only waiting. */
    for(int j = j0 ; j < j1 ; j++) {
        wq_join(a_tasks[j-j0].handle);
        nab_loc += a_tasks[j-j0].nab_loc;
    }
    /* XXX freeing a_tasks is deferred because of course we still need A ! */

    STOPWATCH_GET();
    log_step_time(" done on leaf threads");
    log_step(": sharing among threads");

    struct subtask_info_t * a_tasks2;
    a_tasks2 = malloc(glob.ncores *sizeof(struct subtask_info_t));
    for(int done = 1 ; done < glob.ncores ; done<<=1) {
        for(int j = 0 ; j < glob.ncores ; j += done << 1) {
            if (j + done >= glob.ncores)
                break;
            struct subtask_info_t * task = a_tasks2 + j;
            task->P0 = a_tasks[j].P;
            task->P1 = a_tasks[j+done].P;
            wq_func_t f = (wq_func_t) &a_poly_read_share_child2;
            task->handle = wq_push(glob.wq, f, task);
        }
        /* we're doing nothing, only waiting. */
        for(int j = 0 ; j < glob.ncores ; j += done << 1) {
            if (j + done >= glob.ncores)
                break;
            wq_join(a_tasks2[j].handle);
        }
    }
    free(a_tasks2);

    poly_swap(P, a_tasks[0].P);
    for(int j = j0 ; j < j1 ; j++) poly_free(pols[j-j0]);

    free(a_tasks);

    STOPWATCH_GET();
    log_step_time(" done locally");
    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    allreduce_polymodF_mul_monic(P, glob.acomm, glob.F);
    MPI_Allreduce(MPI_IN_PLACE, &nab_loc, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

    STOPWATCH_GET();
    log_end();

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", nab_loc);
        for(int i = 0 ; i < glob.n ; i++) {
            gmp_fprintf(c->f, "%Zx\n", P->coeff[i]);
        }
        cachefile_close(c);
    }

    return nab_loc;
}/*}}}*/

#if 0
struct tree_like_subtask_info_t {
    struct prime_data * p;
    int j;
    int i0, i1;
    struct wq_task * handle;

    /* In tree-like mode, tasks populate their argument with some info on
     * the child tasks they've spawned. It is of course mandatory to join
     * these tasks as well.
     */
    struct wq_task * t0;
    struct wq_task * t1;
};
#endif

/* {{{ precompute_powers */

void * precompute_powers_child(struct subtask_info_t * info)/* {{{ */
{
    struct prime_data * p = info->p;

    logprint("Precomputing p^%lu, p=%lu\n", glob.prec, p->p);
    // this triggers the whole precomputation.
    power_lookup(p->powers, glob.prec);

    return NULL;
}

/* }}} */

void precompute_powers(struct prime_data * primes, int i0, int i1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    logprint("precompute_powers starts\n");

    {
        struct subtask_info_t * tasks;
        tasks = malloc((i1-i0) * sizeof(struct subtask_info_t));
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            // if (k % glob.psize != glob.prank) continue;
            struct subtask_info_t * task = tasks + k;
            task->p = primes + i;
            task->j = 0;
            wq_func_t f = (wq_func_t) &precompute_powers_child;
            task->handle = wq_push(glob.wq, f, task);
        }
        /* we're doing nothing, only waiting. */
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            // if (k % glob.psize != glob.prank) continue;
            struct wq_task * handle = tasks[k].handle;
            wq_join(handle);
        }
        free(tasks);
    }

    STOPWATCH_GET();

    logprint("precompute_powers ends. wct %.2lf.\n", w1-w0);
}
/*}}}*/

/* {{{ lifting roots */

void * lifting_roots_child(struct subtask_info_t * info)/* {{{ */
{
    struct prime_data * p = info->p;
    int j = info->j;

    mpz_srcptr p1 = power_lookup_const(p->powers, 1);
    mpz_ptr rx = p->lroots->coeff[j];

    cachefile c;

    cachefile_init(c, "lroot_%lu_%lu_%lu", p->p, p->r[j], glob.prec);


    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        gmp_fscanf(c->f, "%Zx", rx);
        cachefile_close(c);
        return NULL;
    }

    mpz_set_ui(rx, p->r[j]);

    mpz_t irx;
    mpz_init(irx);

    mp_poly_eval_mod(irx, glob.f_hat_diff, glob.n-1, rx, p1, NULL);
    mpz_invert(irx, irx, p1);

    logprint("lifting p=%lu, r=%lu\n", p->p, p->r[j]);
    root_lift(p, rx, irx, glob.prec);

    mpz_clear(irx);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        gmp_fprintf(c->f, "%Zx\n", rx);
        cachefile_close(c);
    }

    return NULL;
}

/* }}} */

void lifting_roots(struct prime_data * primes, int i0, int i1)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    //  lifting all roots mod all primes (per pgroup)
    log_begin();

    {
        struct subtask_info_t * lift_tasks;
        lift_tasks = malloc((i1-i0)*n*sizeof(struct subtask_info_t));
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                struct subtask_info_t * task = lift_tasks + k;
                task->p = primes + i;
                task->j = j;
                wq_func_t f = (wq_func_t) &lifting_roots_child;
                task->handle = wq_push(glob.wq, f, task);
            }
        }
        /* we're doing nothing, only waiting. */
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                struct wq_task * handle = lift_tasks[k].handle;
                wq_join(handle);
            }
        }
        free(lift_tasks);
    }

    STOPWATCH_GET();

    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    {
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                struct prime_data * p = primes + i;
                mpz_ptr z = p->lroots->coeff[j];
                broadcast_mpz(z, k % glob.psize, glob.pcomm);
            }
        }
    }
    STOPWATCH_GET();
    log_end();

    for(int i = i0 ; i < i1 ; i++) {
        for(int j = 0 ; j < n ; j++) {
            ASSERT_ALWAYS(mpz_size(primes[i].lroots->coeff[j]));
        }
    }
}/*}}}*/

/* {{{ rational reduction */

/* {{{ building rational product tree */
struct rat_ptree_s {
    struct rat_ptree_s * t0;
    struct rat_ptree_s * t1;
    mpz_t z;
    mpz_srcptr zx;
    struct prime_data * p;
};
typedef struct rat_ptree_s rat_ptree_t;

//  rational product tree: all prime powers.
// it's quite easy to set up.
//
// As such the code does not use the worker threads at all.
// It is possible to do it in a distributed manner, however it is not
// clear that a lot is to be gained.

rat_ptree_t * rat_ptree_build_inner(struct prime_data * p, int i0, int i1)
{
    ASSERT_ALWAYS(i0 < i1);
    rat_ptree_t * res = malloc(sizeof(rat_ptree_t));
    memset(res, 0, sizeof(rat_ptree_t));
    if (i1-i0 == 1) {
        res->p = p + i0;
        res->zx = power_lookup_const(res->p->powers, glob.prec);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = rat_ptree_build_inner(p, i0, i0+d);
    res->t1 = rat_ptree_build_inner(p, i0+d, i1);
    mpz_init(res->z);
    res->zx = res->z;
    WRAP_mpz_mul(res->z, res->t0->zx, res->t1->zx);
    return res;
}

rat_ptree_t * rat_ptree_build(struct prime_data * p, int i0, int i1)
{
    if (i1 == i0)
        return NULL;

    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    rat_ptree_t * res = rat_ptree_build_inner(p, i0, i1);
    STOPWATCH_GET();
    log_end();
    return res;
}

void rat_ptree_clear(rat_ptree_t * t)
{
    if (t == NULL) return;
    rat_ptree_clear(t->t0);
    rat_ptree_clear(t->t1);
    if (t->p == NULL) mpz_clear(t->z);
    free(t);
}

#if 0
unsigned int rat_ptree_nleaves(rat_ptree_t * t)
{
    if (t->p)
        return 1;
    unsigned int n0 = rat_ptree_nleaves(t->t0);
    unsigned int n1 = rat_ptree_nleaves(t->t1);
    return n0 + n1;
}
#endif
/* }}} */

void reduce_poly_mod_rat_ptree(poly_ptr P, rat_ptree_t * T)/*{{{*/
{
    if (T == NULL)
        return;

    // of course this destroys P.
    // we assume that P is reduced mod the top-level T.
    if (T->p) {
        poly_swap(T->p->A, P);
        return;
    }
    poly_t temp;
    poly_alloc(temp, glob.n);
    temp->deg = glob.n - 1;
    for(int i = 0 ; i < glob.n ; i++) {
        WRAP_mpz_mod(temp->coeff[i], P->coeff[i], T->t0->zx);
        WRAP_mpz_mod(P->coeff[i], P->coeff[i], T->t1->zx);
    }
    cleandeg(temp, glob.n - 1);
    cleandeg(P, glob.n - 1);
    reduce_poly_mod_rat_ptree(temp, T->t0);
    reduce_poly_mod_rat_ptree(P, T->t1);
    poly_free(temp);
}/*}}}*/

void * rational_reduction_child(struct subtask_info_t * info)
{
    struct prime_data * primes = info->p;
    int i0 = info->i0;
    int i1 = info->i1;
    size_t off0 = info->off0;
    size_t off1 = info->off1;

    ASSERT_ALWAYS(info->P->deg >= 0);

    rat_ptree_t * ptree = rat_ptree_build(primes, i0, i1);

    if (ptree == NULL)
        ASSERT_ALWAYS(glob.r < glob.ncores);

    if (glob.r <= glob.ncores) {
        ASSERT_ALWAYS(ptree == NULL || ptree->p != NULL);
        if (ptree) {
            for(int i = 0 ; i < glob.n ; i++) {
                WRAP_mpz_mod(ptree->p->A->coeff[i], info->P->coeff[i], ptree->zx);
            }
            cleandeg(ptree->p->A, glob.n-1);
        }
        if (barrier_wait(glob.barrier, NULL, NULL, NULL) == BARRIER_SERIAL_THREAD) {
            poly_t foo;
            poly_alloc(foo, glob.n);
            poly_swap(info->P, foo);
            poly_free(foo);
        }
    } else {
        ASSERT_ALWAYS(ptree != NULL);

        poly_t temp;
        poly_alloc(temp, glob.n);
        for(int i = 0 ; i < glob.n ; i++) {
            WRAP_mpz_mod(temp->coeff[i], info->P->coeff[i], ptree->zx);
        }
        cleandeg(temp, glob.n-1);
        if (barrier_wait(glob.barrier, NULL, NULL, NULL) == BARRIER_SERIAL_THREAD) {
            poly_t foo;
            poly_alloc(foo, glob.n);
            poly_swap(info->P, foo);
            poly_free(foo);
        }
        reduce_poly_mod_rat_ptree(temp, ptree);
        poly_free(temp);
    }

    rat_ptree_clear(ptree);

    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(primes[i].A->deg >= 0);
    }
    if (wcache) {
        for(int i = i0 ; i < i1 ; i++) {
            cachefile c;
            cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu",
                    off0, off1, primes[i].p, glob.prec);
            if (cachefile_open_w(c)) {
                logprint("writing cache %s\n", c->basename);
                /* XXX nab_loc is a misnomer here. In reality we have nab_total
                 * there in this context */
                fprintf(c->f, "%zu\n", info->nab_loc);
                for(int j = 0 ; j < glob.n ; j++)
                    gmp_fprintf(c->f, "%Zx\n", primes[i].A->coeff[j]);
                cachefile_close(c);
            }
        }
    }

    return NULL;
}

void rational_reduction(struct prime_data * primes, int i0, int i1, poly_ptr P, size_t off0, size_t off1, size_t * p_nab_total)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    int rc;

    ASSERT_ALWAYS(P->deg >= 0);

    if (rat_red_caches_ok(primes, i0, i1, off0, off1)) {
        ASSERT_ALWAYS(*p_nab_total == 0);
        for(int i = i0 ; i < i1 ; i++) {
            cachefile c;
            cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu",
                    off0, off1, primes[i].p, glob.prec);
            if (cachefile_open_r(c)) {
                logprint("reading cache %s\n", c->basename);
                rc = fscanf(c->f, "%zu", p_nab_total);
                ASSERT_ALWAYS(rc == 1);
                for(int j = 0 ; j < glob.n ; j++) {
                    rc = gmp_fscanf(c->f, "%Zx", primes[i].A->coeff[j]);
                    ASSERT_ALWAYS(rc == 1);
                }
                cachefile_close(c);
            }
        }
        return;
    }

    {
        struct subtask_info_t * tasks;
        tasks = malloc(glob.ncores*sizeof(struct subtask_info_t));
        for(int k = 0 ; k < glob.ncores ; k++) {
            struct subtask_info_t * task = tasks + k;
            task->p = primes;
            task->i0 = i0 + (i1-i0) * k / glob.ncores;
            task->i1 = i0 + (i1-i0) * (k+1) / glob.ncores;
            task->off0 = off0;
            task->off1 = off1;
            task->P = P;
            task->nab_loc = *p_nab_total;
            wq_func_t f = (wq_func_t) &rational_reduction_child;
            task->handle = wq_push(glob.wq, f, task);
        }
        /* we're doing nothing, only waiting. */
        for(int k = 0 ; k < glob.ncores ; k++) wq_join(tasks[k].handle);
        free(tasks);
    }

    STOPWATCH_GET();
    log_end();

    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(primes[i].A->deg >= 0);
    }
}
/* }}} */

/* {{{ algebraic reduction */

void * algebraic_reduction_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int rc;
    logprint("alg_red (%lu, x-%lu) starts\n", p->p, p->r[j]);

    cachefile c;

    cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu_%lu",
            info->off0, info->off1, info->p->p, info->p->r[j], glob.prec);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        rc = gmp_fscanf(c->f, "%Zx", p->evals->coeff[j]);
        ASSERT_ALWAYS(rc == 1);
        cachefile_close(c);
        return NULL;
    }

    mpz_t ta;
    mpz_init(ta);
    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);
    mpz_set(ta, p->A->coeff[glob.n - 1]);
    for(int k = glob.n - 2 ; k >= 0 ; k--) {
        WRAP_mpz_mul(ta, ta, p->lroots->coeff[j]);
        mpz_add(ta, ta, p->A->coeff[k]);
        WRAP_mpz_mod(ta, ta, px);
    }
    WRAP_mpz_mul(p->evals->coeff[j], p->evals->coeff[j], ta);
    if (mpz_size(p->evals->coeff[j]) >= mpz_size(px) * 3/2)
        WRAP_mpz_mod(p->evals->coeff[j], p->evals->coeff[j], px);
    mpz_clear(ta);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        gmp_fprintf(c->f, "%Zx\n", p->evals->coeff[j]);
        cachefile_close(c);
    }
    return NULL;
}

void algebraic_reduction(struct prime_data * primes, int i0, int i1, size_t off0, size_t off1, size_t * p_nab_total)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();

    if (!alg_red_caches_ok(primes, i0, i1, off0, off1)) {
        for(int i = i0 ; i < i1 ; i++) {
            ASSERT_ALWAYS(primes[i].A->deg >= 0);
        }
    }

    struct subtask_info_t * tasks;
    tasks = malloc((i1-i0) * glob.n * sizeof(struct subtask_info_t));
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            struct subtask_info_t * task = tasks + k;
            task->off0 = off0;
            task->off1 = off1;
            task->p = primes + i;
            task->j = j;
            task->nab_loc = *p_nab_total;
            wq_func_t f = (wq_func_t) &algebraic_reduction_child;
            task->handle = wq_push(glob.wq, f, task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
            * p_nab_total = tasks[k].nab_loc;
        }
    }
    free(tasks);
    for(int i = i0 ; i < i1 ; i++) {
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            mpz_realloc(primes[i].A->coeff[k], 0);
        }
    }

    STOPWATCH_GET();
    log_end();
}/*}}}*/

void multiply_all_shares(struct prime_data * primes, int i0, int i1, size_t * p_nab_total)/*{{{*/
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    int rc;

    // ok. for all eval points, collect the reduced stuff. we'll do the
    // collection in a tree-like manner.
    // this all happens within pcomm, which has size s.

    if (a_shared_caches_ok(primes, i0, i1)) {
        ASSERT_ALWAYS(*p_nab_total == 0);
        for(int i = i0 ; i < i1 ; i++) {
            struct prime_data * p = &primes[i];
            for(int j = 0 ; j < glob.n ; j++) {
                cachefile c;
                cachefile_init(c, "a_mod_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
                int ok = cachefile_open_r(c);
                ASSERT_ALWAYS(ok);
                logprint("reading cache %s\n", c->basename);
                rc = fscanf(c->f, "%zu", p_nab_total);
                ASSERT_ALWAYS(rc == 1);
                rc = gmp_fscanf(c->f, "%Zx", primes[i].evals->coeff[j]);
                ASSERT_ALWAYS(rc == 1);
                cachefile_close(c);
            }
        }
        return ;
    }

    for(int i = i0 ; i < i1 ; i++) {
        mpz_srcptr px = power_lookup_const(primes[i].powers, glob.prec);
        for(int j = 0 ; j < glob.n ; j++) {
            mpz_ptr z = primes[i].evals->coeff[j];
            allreduce_mulmod_mpz(z, glob.pcomm, px);
        }
    }

    STOPWATCH_GET();
    log_end();

    if (wcache) {
        for(int i = i0 ; i < i1 ; i++) {
            struct prime_data * p = &primes[i];
            for(int j = 0 ; j < glob.n ; j++) {
                cachefile c;
                cachefile_init(c, "a_mod_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
                int ok = cachefile_open_w(c);
                ASSERT_ALWAYS(ok);
                logprint("writing cache %s\n", c->basename);
                fprintf(c->f, "%zu\n", * p_nab_total);
                gmp_fprintf(c->f, "%Zx\n", primes[i].evals->coeff[j]);
                cachefile_close(c);
            }
        }
    }
}/*}}}*/

/* {{{ local square roots*/

void * local_square_roots_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int rc;
    cachefile c;
    cachefile_init(c, "sqrt_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
    if (rcache && cachefile_open_r(c)) {
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        logprint("reading cache %s\n", c->basename);
        rc = gmp_fscanf(c->f, "%Zx", p->sqrts->coeff[j]);
        ASSERT_ALWAYS(rc == 1);
        cachefile_close(c);
        return NULL;
    }
    fprintf(stderr, "# [%2.2lf] [P%dA%d] lifting sqrt (%lu, x-%lu)\n", WCT, glob.arank, glob.prank, p->p, p->r[j]);
    sqrt_lift(p, p->evals->coeff[j], p->sqrts->coeff[j], glob.prec);
    // fprintf(stderr, "# [%2.2lf] done\n", WCT);
    mpz_realloc(p->evals->coeff[j], 0);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        gmp_fprintf(c->f, "%Zx\n", p->sqrts->coeff[j]);
        cachefile_close(c);
    }

    return NULL;
}

void local_square_roots(struct prime_data * primes, int i0, int i1, size_t * p_nab_total)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    struct subtask_info_t * tasks;
    tasks = malloc((i1-i0)*glob.n*sizeof(struct subtask_info_t));
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            if (k % glob.psize != glob.prank) continue;
            struct subtask_info_t * task = tasks + k;
            task->p = primes + i;
            task->j = j;
            task->nab_loc = * p_nab_total;
            wq_func_t f = (wq_func_t) &local_square_roots_child;
            task->handle = wq_push(glob.wq, f, task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
            * p_nab_total = tasks[k].nab_loc;
        }
    }
    free(tasks);

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);
    log_step(": sharing");

    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            struct prime_data * p = primes + i;
            mpz_ptr z = p->sqrts->coeff[j];
            broadcast_mpz(z, k % glob.psize, glob.pcomm);
        }
    }

    STOPWATCH_GET();
    log_end();
}

/* }}} */

/* {{{ prime inversion lifts */

void inversion_lift(struct prime_data * p, mpz_ptr iHx, mpz_srcptr Hx, int precision)/* {{{ */
{
    double w0 = WCT;

    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    assert(precision > 0);

    if (precision == 1) {
        mpz_invert(iHx, Hx, pk);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup_const(p->powers, lower);

    // we're going to recurse. change the long Hx_mod value stored by a
    // temporary small one. The problem with this approach is that we
    // store many reductions in memory.
    mpz_t Hx_save;
    mpz_init(Hx_save);
    WRAP_mpz_mod(Hx_save, Hx, pl);
    // recurse.
    inversion_lift(p, iHx, Hx_save, lower);
    mpz_clear(Hx_save);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);

    mpz_t ta;
    mpz_init(ta);

    WRAP_mpz_mul(ta, iHx, Hx);
    WRAP_mpz_mod(ta, ta, pk);

    mpz_sub_ui(ta, ta, 1);
    WRAP_mpz_mul(ta, ta, iHx);
    mpz_sub(iHx, iHx, ta);
    WRAP_mpz_mod(iHx, iHx, pk);

    mpz_clear(ta);
    // gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", WCT, p->iHx_mod);
}/* }}} */

void prime_inversion_lifts_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);

    mpz_t Hx;
    mpz_init(Hx);
    mpz_divexact_ui(Hx, glob.P, p->p);
    mpz_powm_ui(Hx, Hx, glob.prec, px);
    /*
       fprintf(stderr, "# [%2.2lf] size(H^l) %.1f MB\n",
       WCT, mpz_sizeinbase(p->Hx, 256) * 1.0e-6);
       */
    // compute the inverse of H^prec modulo p^prec.
    // need a recursive function for computing the inverse.
    fprintf(stderr, "# [%2.2lf] [P%dA%d] lifting H^-l\n", WCT, glob.arank, glob.prank);
    inversion_lift(p, p->iHx, Hx, glob.prec);
    mpz_clear(Hx);
}

void prime_inversion_lifts(struct prime_data * primes, int i0, int i1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    {
        struct subtask_info_t * tasks;
        tasks = malloc((i1-i0)*sizeof(struct subtask_info_t));
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            if (k % glob.psize != glob.prank)
                continue;
            struct subtask_info_t * task = tasks + k;
            task->p = primes + i;
            task->j = INT_MAX;
            wq_func_t f = (wq_func_t) &prime_inversion_lifts_child;
            task->handle = wq_push(glob.wq, f, task);
        }
        /* we're doing nothing, only waiting. */
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            if (k % glob.psize != glob.prank)
                continue;
            wq_join(tasks[k].handle);
        }
        free(tasks);
    }

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    for(int i = i0 ; i < i1 ; i++) {
        int k = i-i0;
        struct prime_data * p = primes + i;
        broadcast_mpz(p->iHx, k % glob.psize, glob.pcomm);
    }

    STOPWATCH_GET();
    log_end();
}
/* }}} */

/* {{{ prime postcomputations (lagrange reconstruction) */
struct postcomp_subtask_info_t {
    struct prime_data * p;
    int j;
    int64_t * c64;
    mp_limb_t * cN;

    struct wq_task * handle;
};

void prime_postcomputations_child(struct postcomp_subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int64_t * c64 = info->c64;
    mp_limb_t * cN = info->cN;

    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);
    // fprintf(stderr, "# [%2.2lf] done\n", WCT);

    // Lagrange reconstruction.
    //
    // Normally each coefficient has to be divided by the evaluation of
    // the derivative. However we skip this division, effectively
    // reconstructing the polynomial multiplied by the square of the
    // derivative -- which is exactly what we're looking for, in fact.

    // recall that we're working with the number field sieve in mind. So
    // we don't really care about the whole reconstruction in the number
    // field, and from here on we are going to take wild shortcuts.
    // Indeed, even though the ``magical sign combination'' is not known
    // at this point, we do know that the eventual reconstruction will be
    // linear. Thus instead of storing n^2 full length modular integers
    // (n for each root),  and do this for each prime, we store only the
    // pair (quotient mod p^x, residue mod N).

    mpf_t pxf, ratio;
    mpz_t z;

    mpf_init2(pxf, 256);
    mpf_init2(ratio, 256);
    mpz_init(z);

    mpf_set_z(pxf, px);

    mpz_t Hxm;
    mpz_init(Hxm);
    mpz_set_ui(Hxm, p->p);
    mpz_invert(Hxm, Hxm, glob.pol->n);
    mpz_mul(Hxm, Hxm, glob.P);
    mpz_mod(Hxm, Hxm, glob.pol->n);
    mpz_powm_ui(Hxm, Hxm, glob.prec, glob.pol->n);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    /* work mod one root */
    mpz_srcptr rx = p->lroots->coeff[j];
    mpz_ptr sx = p->sqrts->coeff[j];

    ASSERT_ALWAYS(mpz_size(rx));
    ASSERT_ALWAYS(mpz_size(sx));

    // so we have this nice square root. The first thing we do on our
    // list is to scramble it by multiplying it with the inverse of
    // H^x...
    WRAP_mpz_mul(sx, sx, p->iHx);
    WRAP_mpz_mod(sx, sx, px);

    // Now use the evaluation of f_hat mod rx to obtain the lagrange
    // coefficients.
    mpz_set_ui(ta, 1);
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        if (k < glob.n - 1) {
            WRAP_mpz_mul(ta, ta, rx);
            mpz_add(ta, ta, glob.f_hat[k+1]);
            WRAP_mpz_mod(ta, ta, px);
        }
        // multiply directly with H^-x * sqrt
        WRAP_mpz_mul(tb, ta, sx);
        if (k < glob.n - 1) {
            WRAP_mpz_mod(tb, tb, px);
        }
        ASSERT_ALWAYS(mpz_cmp_ui(tb, 0) >= 0);
        ASSERT_ALWAYS(mpz_cmp(tb, px) < 0);

        // now the shortcuts.
        mpf_set_z(ratio, tb);
        mpf_div(ratio, ratio, pxf);
        mpf_mul_2exp(ratio, ratio, 64);
        mpz_set_f(z, ratio);

        uint64_t u;
#if GMP_LIMB_BITS == 64
        u = mpz_get_ui(z);
#else
        u = (uint64_t) mpz_getlimbn(z,1);
        u <<= 32;
        u |= (uint64_t) mpz_getlimbn(z,0);
#endif
        c64[k] = (int64_t) u;

        mpz_mul(tb, tb, Hxm);
        mpz_mod(tb, tb, glob.pol->n);
        mp_size_t sN = mpz_size(glob.pol->n);
        ASSERT_ALWAYS(SIZ(tb) > 0);
        MPN_SET_MPZ(cN + k * sN, sN, tb);
    }

    mpz_clear(ta);
    mpz_clear(tb);
    mpz_clear(Hxm);
    mpz_clear(z);
    mpf_clear(pxf);
    mpf_clear(ratio);
}

void prime_postcomputations(struct prime_data * primes, int i0, int i1, int64_t * contribs64, mp_limb_t * contribsN)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    struct postcomp_subtask_info_t * tasks;
    tasks = malloc((i1-i0)*n*sizeof(struct postcomp_subtask_info_t));
    mp_size_t sN = mpz_size(glob.pol->n);
    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            if (k % glob.psize != glob.prank) continue;
            struct postcomp_subtask_info_t * task = tasks + k;
            task->p = primes + i;
            task->j = j;
            task->c64 = contribs64 + (i * n + j) * n;
            task->cN = contribsN + ((i * n + j) * n) * sN;
            wq_func_t f = (wq_func_t) &prime_postcomputations_child;
            task->handle = wq_push(glob.wq, f, task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
        }
    }
    free(tasks);

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            int64_t * c64 = contribs64 + (i * n + j) * n;
            mp_limb_t * cN = contribsN + ((i * n + j) * n) * sN;
            MPI_Bcast(c64, n, MPI_MY_INT64_T,  k % glob.psize, glob.pcomm);
            MPI_Bcast(cN, n*sN, MPI_MY_MP_LIMB_T, k % glob.psize, glob.pcomm);
        }
    }

    STOPWATCH_GET();
    log_end();
}
void old_prime_postcomputations(int64_t * c64, mp_limb_t * cN, struct prime_data * p)/* {{{ */
{
    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);

    // Lagrange reconstruction.
    //
    // Normally each coefficient has to be divided by the evaluation of
    // the derivative. However we skip this division, effectively
    // reconstructing the polynomial multiplied by the square of the
    // derivative -- which is exactly what we're looking for, in fact.

    // recall that we're working with the number field sieve in mind. So
    // we don't really care about the whole reconstruction in the number
    // field, and from here on we are going to take wild shortcuts.
    // Indeed, even though the ``magical sign combination'' is not known
    // at this point, we do know that the eventual reconstruction will be
    // linear. Thus instead of storing n^2 full length modular integers
    // (n for each root),  and do this for each prime, we store only the
    // pair (quotient mod p^x, residue mod N).

    mpf_t pxf, ratio;
    mpz_t z;

    mpf_init2(pxf, 256);
    mpf_init2(ratio, 256);
    mpz_init(z);

    mpf_set_z(pxf, px);

    mpz_t Hxm;
    mpz_init(Hxm);
    mpz_set_ui(Hxm, p->p);
    mpz_invert(Hxm, Hxm, glob.pol->n);
    mpz_mul(Hxm, Hxm, glob.P);
    mpz_mod(Hxm, Hxm, glob.pol->n);
    mpz_powm_ui(Hxm, Hxm, glob.prec, glob.pol->n);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    // XXX Eh ! mpi-me !
    for(int j = 0 ; j < glob.n ; j++) {
        mpz_srcptr rx = p->lroots->coeff[j];
        mpz_ptr sx = p->sqrts->coeff[j];

        // so we have this nice square root. The first thing we do on our
        // list is to scramble it by multiplying it with the inverse of
        // H^x...
        mpz_mul(sx, sx, p->iHx);
        mpz_mod(sx, sx, px);

        // Now use the evaluation of f_hat mod rx to obtain the lagrange
        // coefficients.
        mpz_set_ui(ta, 1);
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            if (k < glob.n - 1) {
                WRAP_mpz_mul(ta, ta, rx);
                mpz_add(ta, ta, glob.f_hat[k+1]);
                WRAP_mpz_mod(ta, ta, px);
            }
            // multiply directly with H^-x * sqrt
            WRAP_mpz_mul(tb, ta, sx);
            if (k < glob.n - 1) {
                WRAP_mpz_mod(tb, tb, px);
            }
            ASSERT_ALWAYS(mpz_cmp_ui(tb, 0) >= 0);
            ASSERT_ALWAYS(mpz_cmp(tb, px) < 0);

            // now the shortcuts.
            mpf_set_z(ratio, tb);
            mpf_div(ratio, ratio, pxf);
            mpf_mul_2exp(ratio, ratio, 64);
            mpz_set_f(z, ratio);

            uint64_t u;
#if GMP_LIMB_BITS == 64
            u = mpz_get_ui(z);
#else
            u = (uint64_t) mpz_getlimbn(z,1);
            u <<= 32;
            u |= (uint64_t) mpz_getlimbn(z,0);
#endif
            c64[j*glob.n+k] = (int64_t) u;

            mpz_mul(tb, tb, Hxm);
            mpz_mod(tb, tb, glob.pol->n);
            mp_size_t sN = mpz_size(glob.pol->n);
            ASSERT_ALWAYS(SIZ(tb) > 0);
            MPN_SET_MPZ(cN + (j*glob.n+k) * sN, sN, tb);
        }
    }
    mpz_clear(ta);
    mpz_clear(tb);
    mpz_clear(Hxm);
    mpz_clear(z);
    mpf_clear(pxf);
    mpf_clear(ratio);
}
/* }}} */

/* }}} */

void mpi_set_communicators()/*{{{*/
{
    int s = glob.s;
    int t = glob.t;
    if (s * t > glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Error: need at least %d*%d=%d jobs (or more RAM)\n",
                    s,t,s*t);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else if (s * t < glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Warning: only %d*%d=%d jobs needed, have %d\n",
                    s,t,s*t,glob.nprocs);
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank / t, glob.rank, &glob.acomm);
    MPI_Comm_rank(glob.acomm, &glob.arank);
    MPI_Comm_size(glob.acomm, &glob.asize);
    ASSERT_ALWAYS(glob.asize == t);

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank % t, glob.rank, &glob.pcomm);
    MPI_Comm_rank(glob.pcomm, &glob.prank);
    MPI_Comm_size(glob.pcomm, &glob.psize);
    ASSERT_ALWAYS(glob.psize == s);
}
/*}}}*/

void banner()
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (glob.rank == 0) {
        fprintf(stderr, "###########################################################################\n");
    }
}

int main(int argc, char **argv)
{
    int ret, i;
    int asked_r = 0;
    // int size_guess = 0;

    MPI_Init(&argc, &argv);
    program_starttime = wct_seconds();

    MPI_Comm_rank(MPI_COMM_WORLD, &glob.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &glob.nprocs);

    /* {{{ parameter parsing */
    /* print the command line */
    fprintf(stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");

    param_list pl;
    param_list_init(pl);
    int cache=0;
    param_list_configure_knob(pl, "-v", &verbose);
    param_list_configure_knob(pl, "--cache", &cache);
    param_list_configure_knob(pl, "--rcache", &rcache);
    param_list_configure_knob(pl, "--wcache", &wcache);
    // param_list_configure_knob(pl, "--size-guess", &size_guess);
    int wild = 0;
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        if (argv[0][0] != '-' && wild == 0) {
            param_list_add_key(pl, "depfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 1) {
            param_list_add_key(pl, "ratdepfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 2) {
            param_list_add_key(pl, "polyfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    if (cache)  rcache = wcache = 1;

    param_list_parse_double(pl, "print_delay", &print_delay);
    param_list_parse_double(pl, "ram", &ram_gb);

    param_list_parse_int(pl, "ncores", &glob.ncores);
    param_list_parse_int(pl, "r", &asked_r);
    param_list_parse_int(pl, "lll_maxdim", &glob.lll_maxdim);

    if (param_list_lookup_string(pl, "depfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "ratdepfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "polyfile") == NULL)
        usage();
    /* }}} */

    cado_poly_init(glob.pol);
    ret = cado_poly_read(glob.pol, param_list_lookup_string(pl, "polyfile"));
    glob.n = glob.pol->degree;
    ASSERT_ALWAYS(ret == 1);

    /* {{{ create f_hat, the minimal polynomial of alpha_hat = lc(f) *
     * alpha */
    {
        mpz_t tmp;
        mpz_init_set_ui(tmp, 1);

        glob.f_hat = malloc((glob.n + 1) * sizeof(mpz_t));
        glob.f_hat_diff = malloc(glob.n * sizeof(mpz_t));
        mpz_init_set_ui(glob.f_hat[glob.n], 1);
        for(int i = glob.n - 1 ; i >= 0 ; i--) {
            mpz_init(glob.f_hat[i]);
            mpz_init(glob.f_hat_diff[i]);
            mpz_mul(glob.f_hat[i], tmp, glob.pol->f[i]);
            mpz_mul(tmp, tmp, glob.pol->f[glob.n]);
            mpz_mul_ui(glob.f_hat_diff[i], glob.f_hat[i+1], i+1);
        }
        mpz_clear(tmp);
        if (glob.rank == 0)
            fprintf(stderr, "# [%2.2lf] Note: all computations done with polynomial f_hat\n", WCT);
    }
    /* }}} */

    // FIXME: merge these two instances of the algebraic polynomial in
    // one type only.
    poly_alloc(glob.F, glob.n);
    for (int i = glob.pol->degree; i >= 0; --i)
        poly_setcoeff(glob.F, i, glob.f_hat[i]);

    // fprintf(stderr, "# [%2.2lf] A is f_d^%zu*f_hat'(alpha_hat)*prod(f_d a - b alpha_hat)\n", WCT, nab + (nab &1));

    ab_source_init(glob.ab, param_list_lookup_string(pl, "depfile"),
            glob.rank, 0, MPI_COMM_WORLD);

    // note that for rsa768, this estimation takes only 10 minutes, so
    // it's not a big trouble.
    unsigned long toto;
    if (!param_list_parse_ulong(pl, "sqrt_coeffs_bits", &toto)) {
        if (glob.rank == 0) {
            estimate_nbits_sqrt(&glob.nbits_sqrt, glob.ab); //, size_guess);
        }
        MPI_Bcast(&glob.nbits_sqrt, 1, MPI_MY_SIZE_T, 0, MPI_COMM_WORLD);
    } else {
        glob.nbits_sqrt = toto;
    }
    // MPI_Bcast(&glob.nbits_a, 1, MPI_MY_SIZE_T, 0, MPI_COMM_WORLD);
    // we no longer need to know nab, so let's drop it as a proof !
    glob.ab->nab = 0;

    glob.r=0;
    glob.s=0;
    glob.t=0;
    if (glob.rank == 0) {
        get_parameters(&glob.r,&glob.s,&glob.t,asked_r);
    }
    MPI_Bcast(&glob.s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&glob.t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&glob.r, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int r = glob.r;
    int s = glob.s;
    // int n = glob.n;
    int t = glob.t;
    glob.m = r * t;


    // this is sufficiently trivial to do.
    struct prime_data * primes = suitable_crt_primes();

    mpz_init_set_ui(glob.P, 1);
    double log2_P = 0;
    for(int i = 0 ; i < glob.m ; i++) {
        log2_P += log(primes[i].p)/M_LN2;
        mpz_mul_ui(glob.P, glob.P, primes[i].p);
    }
    glob.prec = ceil((glob.nbits_sqrt + 128)/ log2_P);

    ASSERT_ALWAYS(mpi_data_agrees(&glob.prec, 1, MPI_INT, MPI_COMM_WORLD));

    if (glob.rank == 0) {
        char sbuf[32];
        fprintf(stderr, "# [%2.2lf] Lifting to precision l=%d (p^l is approx %s)\n", WCT, glob.prec, size_disp(glob.prec * log(primes[0].p)/M_LN2 / 8, sbuf));
    }

    mpi_set_communicators();

    if (glob.rank == 0) {
        fprintf(stderr, "# [%2.2lf] starting %d worker threads on each node\n", WCT, glob.ncores);
    }
    wq_init(glob.wq, glob.ncores);
    barrier_init(glob.barrier, glob.ncores);

    int pgnum = glob.rank % t;
    int apnum = glob.rank / t;

    snprintf(prefix, sizeof(prefix), "[P%dA%d] ", glob.arank, glob.prank);

    ASSERT_ALWAYS(glob.arank == pgnum);
    ASSERT_ALWAYS(glob.prank == apnum);

    int i0 = pgnum*r;
    int i1 = i0 + r;

    size_t off0 = apnum * glob.ab->totalsize / s;
    size_t off1 = (apnum+1) * glob.ab->totalsize / s;

    for(int i = i0 ; i < i1 ; i++) {
        prime_initialization(&(primes[i]));
        for(int j = 0 ; j < glob.n ; j++) {
            mpz_set_ui(primes[i].evals->coeff[j], 1);
        }
    }

    // unless we have really nothing to do, we will need these two steps:

    banner(); /*********************************************/
    precompute_powers(primes, i0, i1);

    banner(); /*********************************************/
    lifting_roots(primes, i0, i1);

    // find some info on the caches.
    int sqrt_caches = sqrt_caches_ok(primes, i0, i1);
    int shared_caches = a_shared_caches_ok(primes, i0, i1);
    int alg_caches = alg_red_caches_ok(primes, i0, i1, off0, off1);
    int rat_caches = rat_red_caches_ok(primes, i0, i1, off0, off1);

    if (rcache) {
        logprint("Caches: sqrt %s, shared %s, alg %s, rat %s\n",
                sqrt_caches ? "ok" : "no",
                shared_caches ? "ok" : "no",
                alg_caches ? "ok" : "no",
                rat_caches ? "ok" : "no");
    }

    size_t nab = 0;

    poly_t P;
    poly_alloc(P, glob.n);

    int next = 0;
    if (sqrt_caches) {
        next = 4;
    } else if (shared_caches) {
        next = 3;
    } else if (alg_caches) {
        next = 2;
    } else if (rat_caches) {
        next = 1;
    } else {
        next = 0;
    }

    logprint("Next step is: %d\n", next);

    if (next < 1) {
        banner(); /*********************************************/
        nab = a_poly_read_share(P, off0, off1);
    }

    if (next < 2) {
        banner(); /*********************************************/
        rational_reduction(primes, i0, i1, P, off0, off1, &nab);
    }

    if (next < 3) {
        banner(); /*********************************************/
        algebraic_reduction(primes, i0, i1, off0, off1, &nab);
    }

    if (next < 4) {
        banner(); /*********************************************/
        multiply_all_shares(primes, i0, i1, &nab);
        // XXX This is a hack. we are evaluating the products
        // f_d^\epsilon*A*f_hat'(\hat\alpha)
        //
        // where A is the value denoted by the variable named A. It
        // is defined as:
        //
        // A(\hat\alpha) = (\prod_{(a,b)}(f_da-b\hat\alpha)
        //
        // Instead of fixing what's missing in A, we use shortcuts.
        // Extra f_d coefficients are added to the evaluation array,
        // and the derivative is eliminated later on by avoiding the
        // normalization in the Lagrange step.
        if (nab & 1) {
            logprint("Correcting by one leading term\n");
            for(int i = i0 ; i < i1 ; i++) {
                for(int j =  0 ; j < glob.n ; j++) {
                    mpz_ptr x = primes[i].evals->coeff[j];
                    mpz_mul(x, x, glob.pol->f[glob.n]);
                }
            }
        }
    }

    banner(); /*********************************************/
    local_square_roots(primes, i0, i1, &nab);

    logprint("Number of pairs is %lu\n", nab);

    poly_free(P);

    banner(); /*********************************************/
    prime_inversion_lifts(primes, i0, i1);

    banner(); /*********************************************/
    size_t nc = glob.m * glob.n * glob.n;
    int64_t * contribs64 = malloc(nc * sizeof(int64_t));
    mp_size_t sN = mpz_size(glob.pol->n);
    memset(contribs64,0,nc * sizeof(int64_t));
    mp_limb_t * contribsN = malloc(nc * sN * sizeof(mp_limb_t));
    memset(contribsN, 0,  nc * sN * sizeof(mp_limb_t)); 
#if 1
    prime_postcomputations(primes, i0, i1, contribs64, contribsN);
    // now share the contribs !
    for(int i = 0 ; i < glob.m ; i++) {
        int root = i / r;
        int d64 = glob.n * glob.n;
        int dN = d64 * mpz_size(glob.pol->n);
        MPI_Bcast(contribs64 + i * d64, d64, MPI_MY_INT64_T, root, glob.acomm);
        MPI_Bcast(contribsN + i * dN, dN, MPI_MY_MP_LIMB_T, root, glob.acomm);
    }
#else
    if (glob.prank == 0) {

        for(int i = i0 ; i < i1 ; i++) {
            int disp = i * glob.n * glob.n;
            old_prime_postcomputations(contribs64 + disp, contribsN + disp * sN, &primes[i]);
        }
    }

    // now share the contribs !
    if (glob.prank == 0) {
        for(int i = 0 ; i < glob.m ; i++) {
            int root = i / r;
            int d64 = glob.n * glob.n;
            int dN = d64 * mpz_size(glob.pol->n);
            MPI_Bcast(contribs64 + i * d64, d64, MPI_MY_INT64_T, root, glob.acomm);
            MPI_Bcast(contribsN + i * dN, dN, MPI_MY_MP_LIMB_T, root, glob.acomm);
        }
    }
#endif

    if (glob.rank == 0) {
        fprintf(stderr, "# [%2.2lf] clearing work queues\n", WCT);
    }
    wq_clear(glob.wq, glob.ncores);
    barrier_destroy(glob.barrier);

    for(int i = i0 ; i < i1 ; i++) prime_cleanup(&(primes[i]));
    free(primes);

    banner(); /*********************************************/
    // and only the leader does the knapsack.
    if (glob.rank == 0) {
        struct crtalgsqrt_knapsack cks[1];
        crtalgsqrt_knapsack_init(cks);
        crtalgsqrt_knapsack_prepare(cks, (nab+(nab&1)) / 2);
        cks->ks->tab = contribs64;
        cks->tabN = contribsN;
        knapsack_solve(cks->ks);
        crtalgsqrt_knapsack_clear(cks);
    }

    /****************************************************************/

    free(contribs64);
    free(contribsN);

    {
        for(int i = glob.n ; i >= 0 ; i--) {
            mpz_clear(glob.f_hat[i]);
        }
        for(int i = glob.n - 1 ; i >= 0 ; i--) {
            mpz_clear(glob.f_hat_diff[i]);
        }
        free(glob.f_hat);
        free(glob.f_hat_diff);
    }

    mpz_clear(glob.P);
    cado_poly_clear(glob.pol);
    poly_free(glob.F);

    ab_source_clear(glob.ab);
    param_list_clear(pl);

    MPI_Finalize();

    return 0;
}

#if 0/*{{{*/

// Init F to be the algebraic polynomial
poly_alloc(F, degree);
for (i = pol->degree; i >= 0; --i)
poly_setcoeff(F, i, pol->f[i]);

// Init prd to 1.
poly_alloc(prd->p, pol->degree);
mpz_set_ui(prd->p->coeff[0], 1);
prd->p->deg = 0;
prd->v = 0;

// Allocate tmp
poly_alloc(tmp->p, 1);

// Accumulate product
int nab = 0, nfree = 0;
#if 0
// Naive version, without subproduct tree
    while (fscanf(depfile, "%ld %lu", &a, &b) != EOF) {
        if (!(nab % 100000))
            fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab,
                    WCT);
        if ((a == 0) && (b == 0))
            break;
        polymodF_from_ab(tmp, a, b);
        polymodF_mul(prd, prd, tmp, F);
        nab++;
    }
#else
// With a subproduct tree
{
    polymodF_t *prd_tab;
    unsigned long lprd = 1;	/* number of elements in prd_tab[] */
    unsigned long nprd = 0;	/* number of accumulated products in prd_tab[] */
    prd_tab = (polymodF_t *) malloc(lprd * sizeof(polymodF_t));
    poly_alloc(prd_tab[0]->p, F->deg);
    mpz_set_ui(prd_tab[0]->p->coeff[0], 1);
    prd_tab[0]->p->deg = 0;
    prd_tab[0]->v = 0;
    while (fscanf(depfile, "%ld %lu", &a, &b) != EOF) {
        if (!(nab % 100000))
            fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab,
                    WCT);
        if ((a == 0) && (b == 0))
            break;
        polymodF_from_ab(tmp, a, b);
        prd_tab = accumulate_fast(prd_tab, tmp, F, &lprd, nprd++);
        nab++;
        if (b == 0)
            nfree++;
    }
    fprintf(stderr, "# Read %d including %d free relations\n", nab,
            nfree);
    ASSERT_ALWAYS((nab & 1) == 0);
    ASSERT_ALWAYS((nfree & 1) == 0);
    accumulate_fast_end(prd_tab, F, lprd);
    fclose(depfile);

    poly_copy(prd->p, prd_tab[0]->p);
    prd->v = prd_tab[0]->v;
    for (i = 0; i < (long) lprd; ++i)
        poly_free(prd_tab[i]->p);
    free(prd_tab);
}
#endif

fprintf(stderr, "Finished accumulating the product at %2.2lf\n",
        WCT);
fprintf(stderr, "nab = %d, nfree = %d, v = %d\n", nab, nfree, prd->v);
fprintf(stderr, "maximal polynomial bit-size = %lu\n",
        (unsigned long) poly_sizeinbase(prd->p, pol->degree - 1, 2));

p = FindSuitableModP(F);
fprintf(stderr, "Using p=%lu for lifting\n", p);

double tm = seconds();
polymodF_sqrt(prd, prd, F, p);


fprintf(stderr, "Square root lifted in %2.2lf\n", seconds() - tm);
mpz_t algsqrt, aux;
mpz_init(algsqrt);
mpz_init(aux);
poly_eval_mod_mpz(algsqrt, prd->p, pol->m, pol->n);
mpz_invert(aux, F->coeff[F->deg], pol->n);	// 1/fd mod n
mpz_powm_ui(aux, aux, prd->v, pol->n);	// 1/fd^v mod n
mpz_mul(algsqrt, algsqrt, aux);
mpz_mod(algsqrt, algsqrt, pol->n);

gmp_fprintf(stderr, "Algebraic square root is: %Zd\n", algsqrt);

depfile = fopen(param_list_lookup_string(pl, "ratdepfile"), "r");
ASSERT_ALWAYS(depfile != NULL);
gmp_fscanf(depfile, "%Zd", aux);
fclose(depfile);

gmp_fprintf(stderr, "Rational square root is: %Zd\n", aux);
fprintf(stderr, "Total square root time is %2.2lf\n", seconds());

int found = 0;
{
    mpz_t g1, g2;
    mpz_init(g1);
    mpz_init(g2);

    // First check that the squares agree
    mpz_mul(g1, aux, aux);
    mpz_mod(g1, g1, pol->n);
    if (mpz_cmp_ui(pol->g[1], 1) != 0) {
        // case g(X)=m1*X+m2 with m1 != 1
        // we should have prod (a+b*m2/m1) = A^2 = R^2/m1^(nab-nfree)
        // and therefore nab should be even
        if (nab & 1) {
            fprintf(stderr, "Sorry, but #(a, b) is odd\n");
            fprintf(stderr,
                    "Bug: this should be patched! Please report your buggy input\n");
            printf("Failed\n");
            return 0;
        }
        mpz_powm_ui(g2, pol->g[1], (nab - nfree) >> 1, pol->n);
        mpz_mul(algsqrt, algsqrt, g2);
        mpz_mod(algsqrt, algsqrt, pol->n);
    }
    mpz_mul(g2, algsqrt, algsqrt);
    mpz_mod(g2, g2, pol->n);
    if (mpz_cmp(g1, g2) != 0) {
        fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
        //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
    }
    mpz_sub(g1, aux, algsqrt);
    mpz_gcd(g1, g1, pol->n);

    if (mpz_cmp(g1, pol->n)) {
        if (mpz_cmp_ui(g1, 1)) {
            found = 1;
            gmp_printf("%Zd\n", g1);
        }
    }
    mpz_add(g2, aux, algsqrt);
    mpz_gcd(g2, g2, pol->n);
    if (mpz_cmp(g2, pol->n)) {
        if (mpz_cmp_ui(g2, 1)) {
            found = 1;
            gmp_printf("%Zd\n", g2);
        }
    }
    mpz_clear(g1);
    mpz_clear(g2);
}
if (!found)
    printf("Failed\n");

    cado_poly_clear(pol);
    mpz_clear(aux);
    mpz_clear(algsqrt);
    poly_free(F);
    poly_free(prd->p);
    poly_free(tmp->p);

    return 0;
#endif/*}}}*/
