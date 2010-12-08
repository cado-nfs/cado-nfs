#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "parallelizing_info.h"
#include "select_mpi.h"

#include "utils.h"
// #include "mf.h"

#include "balancing.h"
#include "balancing_data_source.h"
#include "balancing_file_source.h"
#include "balancing_mpi_source.h"

#ifdef  HAVE_CURL
#include "balancing_curl_source.h"
#endif

#include "balancing_workhorse.h"

/* TODO: implement file-backed rewind on thread pipes. There's a
 * balancing_use_auxfile parameter for this. */
/* TODO: integrate fully with matmul_top.c, and also call mf_scan and
 * mf_bal when needed (we also want a standalone tool, though).
 */

/**********************************************************************/
/* {{{ doc */
/* Balancing a matrix can be decomposed as a 5-step process. 
 *
 *     1) read the original matrix. collect row and column frequency
 *     data.  This has to be done at some point, but clearly, since it is
 *     balancing independent, for large experiments it is possibly useful
 *     to save it.  This is done within mf_scan (but could optionally be
 *     called from here as well).
 *
 *     2) from this frequency data, and the information on which kind of
 *     balancing to perform, compute the row and column permutation which
 *     must be applied to the original matrix. Note that this step does
 *     _NOT_ require the original matrix to be read ! It is performed by
 *     mf_dobal (but could optionally be called from here as well).
 * 
 * The code here handles the two remaining steps. Precisions are given
 * later on.
 *
 *     3) ``master loop''. From the row and column permutation data, and
 *     from the matrix, read input coefficient by coefficient, and deduce
 *     who goes where.
 *
 *     4) ``slave loop 1''. This loop is executed a priori for each
 *     submatrix corresponding to the balancing. Scan the input quickly
 *     first, just so as to count row and column frequencies in each
 *     submatrix.
 *
 *     5) ``slave loop 2''. Process again the preceding input, but given
 *     the knowledge acquired from the previous step, it is possible to
 *     place pointers for accessing the matrix in core memory, in order
 *     to write to it. Transposing can be performed here at no extra
 *     cost. The output is produced in a format which appears ready for
 *     the matmul_build routines.
 *
 * Some notes on steps 3--5
 *
 *     ``going'', in step 3, can have different meanings depending on
 *     whether we're talking big business or not. For small to moderate
 *     size examples where pthreads are enough, it just means passing on
 *     to the next processing step. Otherwise it may incur sending via
 *     MPI.
 * 
 *     The data flow from step 3 to step 4 may be either a message
 *     passing scenario (thus the input appears from an API call, MPI in
 *     this case), or a shortcut: step 3 output may be hooked directly
 *     to step 4 processing.
 *
 *     Rewinding the input for doing step 5 may be tricky. If the input
 *     is to be kept in RAM, then the memory footprint is not necessarily
 *     acceptable (note that the total output set is (almost) guaranteed
 *     to fit, though, at least provided that the machine (or cluster) on
 *     which the balancing is done corresponds to the one which actually
 *     runs the linear algebra eventually.
 *     
 *     One thing to note regarding rewinding is that although it is by no
 *     means accurate, it is possible to _estimate_ the size of each
 *     submatrix.
 */

/* There are several modes of operations considered. Note that the
 * program has to provide some basic functionality even in the absence of
 * an available MPI environment (threads are required though).
 *
 * 1) The matrix has to be split according to an m*n grid, on a single
 * node. The linear algebra will be performed using m*n threads on that
 * same node. The input used is the matrix file, and the balancing file.
 *    1a) Assume that the binary matrix fits _twice_ in memory.
 *        In this case, we have.
 *          - a master loop, with the file as a source, and a dispatching
 *          destination. The (several) destinations are all (as in all
 *          other cases) handled by the data_dest object. The
 *          implementation in this case flushes its input buffer by
 *          feeding it to a companion thread which works as a consumer.
 *          - the consumer threads do exactly what is done in the step
 *          called ``slave loop 1 above''. Henceforth, the must have the
 *          capability of rewinding their input.
 *          XXX => the data_dest buffers are all parts of big
 *          per-thread buffers. they are tightly bound to the source
 *          which is used for the slave loop 2, hence in reality it is
 *          the same object.
 *    1b) Same, but now we assume that the binary matrix does _not_ fit
 *        twice in memory.
 *        This case is quite similar to the above, and differ only by the
 *        rewinding mechanism. data_dest buffers are tiny scratch space, and
 *        everything that passes through these buffers also goes to a temp
 *        file. (XXX to be done -- likely to be an mmap eventually)
 * 2) The matrix has to be split according to an m*n grid, on m*n nodes.
 *    Then an MPI communication channel is necessary. Thus the slave loops
 *    ``see'' their input as coming from MPI_Recv (one may want to consider
 *    the case where the job which reads the original matrix plays both
 *    roles. Then the above approach may be used instead).
 *    2a) everything fits in RAM twice.
 *        in this case, rewinding is performed by duplicating the source
 *        via a buffer. We thus have a rewindable input stream, built on
 *        top of one which is not. In order to avoid extraneous copies,
 *        we arrange for the mpi_source data to put data *in the buffer
 *        provided* instead of the buffer from the data_source struct in
 *        case the latter has not been allocated.
 *    2b) A temp file is needed.
 *        Then the temp file is created, from the same ``wrapped''
 *        mpi_source buffer.
 * 3) both mpi and pthreads. the situation is complicated by the fact
 * the MPI_PTHREAD_MULTIPLE only barely works most of the time, thus we
 * do not expect that several threads may successfully communicate
 * simultaneously.
 *    (not yet implemented)
 *    per node, only one mpi_source exists. it receives with MPI_TAG_ANY.
 *    mpi_outchannels are merged, and the tags used for sending encode
 *    the receiving sub-thread.
 *    XXX MPI messages are guaranteed not to overtake eachother. thus
 *    it's fair to use the tag for this, apparently.
 *    The difficulty here is that the receiving side can't predict for
 *    whom it's going to receive stuff. At the moment I can't see how to
 *    avoid a copy, which is something one can live with anyway (since
 *    it's downstream after an MPI communication).
 *    * duplicating:
 *    making a rewindable buffer in this situation is not obvious. The
 *    goal is to make slave sloop 2 thread agnostic. Thus buffer
 *    duplicating (or saving to tempfile) must be done in a way which is
 *    similar to case 1 above, presumably.
 *
 * Note that situation 3 above can be simulated by an mpi-only approach,
 * and that's decent enough for a first plan. However, if one wants to
 * use this program as a routine called in a simd-fashion directly from
 * within the matmul-top routines, it's not possible.
 */
/* }}} */

// int transposing = 0;
// const char * bfile;

/* The number of send queues allocated on the master is 2 *
 * n_peer_threads */
// both sizes below are in uint32_t's
const size_t default_queue_size = 1 << 20;

// {{{ trivial utility
static const char *size_disp(size_t s, char buf[16])
{
    char *prefixes = "bkMGT";
    double ds = s;
    const char *px = prefixes;
    for (; px[1] && ds > 500.0;) {
	ds /= 1024.0;
	px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}

// }}}

/* {{{ sources  */

/* }}} */

/* {{{ dests. Impls are task-specific, thus defined later.  */

struct data_dest_s {/*{{{*/
    int (*put)(void*,uint32_t *, size_t);
    // void (*progress)(void*, time_t);
    // number of rows which passed through. Pretty ad-hoc to have it
    // here, indeed.
    uint32_t r;
    size_t pos;
};
typedef struct data_dest_s data_dest[1];
typedef struct data_dest_s *data_dest_ptr;/*}}}*/

/*}}}*/

/* {{{ dual-head thread pipes */
struct thread_pipe_s;
struct thread_source_s;
struct thread_dest_s;

struct thread_pipe_s {
    pthread_cond_t hello;
    pthread_mutex_t mu;
    uint32_t * buf;
    int mmaped;
    size_t batchsize;
    size_t size;
    size_t avail;       // available for reading.
    size_t total;       // set to dst->b->pos on close
    // these are the two facets of the same object.
    struct thread_source_s * src;
    struct thread_dest_s * dst;
};
typedef struct thread_pipe_s *thread_pipe_ptr;
typedef struct thread_pipe_s thread_pipe[1];

struct thread_dest_s {
    data_dest b;
    thread_pipe_ptr t;
};
typedef struct thread_dest_s *thread_dest_ptr;
typedef struct thread_dest_s thread_dest[1];

struct thread_source_s {
    data_source b;
    thread_pipe_ptr t;
};
typedef struct thread_source_s *thread_source_ptr;
typedef struct thread_source_s thread_source[1];

int thread_dest_put(thread_dest_ptr D, uint32_t * p, size_t n)
{
    thread_pipe_ptr t = D->t;
    FATAL_ERROR_CHECK(D->b->pos + n > t->size,
        "Argh, error in estimation of the per-core memory\n");
    /* It is presently *NOT* possible to realloc() here, since mf_pipe
     * keeps the pointer it has received from the other side of the
     * thread_pipe, and passes it to slave_dest_put. Fixing this would
     * require changing the structures a bit (probably doable without too
     * much harm).
     */

    ASSERT_ALWAYS(D->b->pos + n <= t->size);
    if (n) {
        /* This memcpy should be the only data duplication in the
         * critical path. It's avoidable when there is only one thread,
         * but it's not done so yet. When there is only one mpi job, this
         * memcpy is always tiny (within a cache line).
         */
        memcpy(t->buf + D->b->pos, p, n * sizeof(uint32_t));
        D->b->pos += n;
    }
    ASSERT_ALWAYS(t->avail % t->batchsize == 0);

    if (p == NULL && n == 0) {
        pthread_mutex_lock(&t->mu);
        t->avail = D->b->pos;
        t->total = D->b->pos;
        pthread_cond_signal(&t->hello);
        pthread_mutex_unlock(&t->mu);
        return -1;
    } else if (D->b->pos - t->avail > t->batchsize) {
        size_t more = (D->b->pos - D->b->pos % t->batchsize) - t->avail;
        ASSERT_ALWAYS(more > 0);
        pthread_mutex_lock(&t->mu);
        t->avail += more;
        pthread_cond_signal(&t->hello);
        pthread_mutex_unlock(&t->mu);
    }
    return 0;
}

thread_dest_ptr thread_dest_alloc(thread_pipe_ptr t)
{
    thread_dest_ptr res = malloc(sizeof(thread_dest));
    memset(res, 0, sizeof(thread_dest));
    res->t = t;
    res->b->put = (int(*)(void*,uint32_t*,size_t))thread_dest_put;
    return res;
}

size_t thread_source_get(thread_source_ptr s, uint32_t ** p, size_t avail)
{
    ASSERT_ALWAYS(avail == 0);
    ASSERT_ALWAYS(s->b->pos <= s->t->avail);
    size_t n = s->t->avail - s->b->pos;
    *p = s->t->buf+s->b->pos;
    if (n) {
        s->b->pos += n;
        return n;
    }
    pthread_mutex_lock(&s->t->mu);
    for( ; (n = s->t->avail - s->b->pos) == 0 && s->t->avail < s->t->total; ) {
        pthread_cond_wait(&s->t->hello, &s->t->mu);
        ASSERT_ALWAYS(s->b->pos <= s->t->avail);
    }
    s->b->pos += n;
    pthread_mutex_unlock(&s->t->mu);
    return n;
}

void thread_source_rewind(thread_source_ptr s)
{
    s->b->pos = 0;      /* easy enough ;-) */
}

thread_source_ptr thread_source_alloc(thread_pipe_ptr t)
{
    thread_source_ptr res = malloc(sizeof(thread_source));
    memset(res, 0, sizeof(thread_source));
    res->t = t;
    res->b->get = (size_t(*)(void*,uint32_t**,size_t))thread_source_get;
    return res;
}

void thread_source_free(thread_source_ptr x)
{
    free(x);
}

void thread_dest_free(thread_dest_ptr x)
{
    free(x);
}


// this is a RAM-backed pipe. the corresponding source is rewindable.
thread_pipe_ptr thread_pipe_alloc(size_t size, int use_auxfile)
{
    ASSERT_ALWAYS(use_auxfile == 0);    // TODO !
    /* must be called by only _one_ thread ! */
    thread_pipe_ptr res = malloc(sizeof(thread_pipe));
    memset(res, 0, sizeof(thread_pipe));
    pthread_cond_init(&res->hello, NULL);
    pthread_mutex_init(&res->mu, NULL);
    res->batchsize = 65536;
    res->total = SIZE_MAX;
    res->size = size;
    res->mmaped = use_auxfile;
    res->buf = malloc(size * sizeof(uint32_t));
    res->dst = thread_dest_alloc(res);
    res->src = thread_source_alloc(res);
    return res;
}

void thread_pipe_free(thread_pipe_ptr tp)
{
    free(tp->buf);
    pthread_cond_destroy(&tp->hello);
    pthread_mutex_destroy(&tp->mu);
    thread_dest_free(tp->dst);
    thread_source_free(tp->src);
    memset(tp, 0, sizeof(thread_pipe));
    free(tp);
}

/* }}} */

/*{{{ pipes */
void mf_progress(data_source_ptr s, data_dest_ptr d, time_t dt, const char * name)
{
    char buf[16];
    char buf2[16];
    printf("%s: %s, %" PRIu32 " rows in %d s ; %s/s  \r",
            name,
            size_disp(s->pos * sizeof(uint32_t), buf), d->r, (int) dt,
            size_disp(dt > 0 ? (size_t) (s->pos * sizeof(uint32_t) / dt) : 0, buf2));
    fflush(stdout);
}

void mf_pipe(data_source_ptr input, data_dest_ptr output, const char * name)/*{{{*/
{
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;
    time_t t = t0;

    for (;;) {
        uint32_t * ptr = NULL;  // we are not providing a buffer.
        size_t n = input->get(input, &ptr, 0);
        if (n == 0) break;
        int r = output->put(output, ptr, n);
	t = time(NULL);
        if (r < 0)
            break;
	if (t >= t1) {
	    t1 = t + 1;
            mf_progress(input, output, t-t0, name);
	}
    }
    mf_progress(input, output, t-t0, name);
    printf("\n");
}/*}}}*/

/* }}} */




// _outbound_ communication channel./*{{{*/

struct mpi_dest_s;

// all outbound channels use the same type of send queue. There is a
// double-size buffer, and flushes are performed each time one of the two
// gets filled (while the other one continues filling up normally).
// only the flush operations are virtual.
struct mpi_dest_s {
    data_dest b;

    uint32_t *buf[2];
    int cur;

    size_t size;
    size_t pos;
    size_t total;
    size_t tailsize;

    time_t pending;
    time_t waited;
    time_t total_io_wct;	// including offloaded time.

    int tag;
    MPI_Request req;
    int peer;

    void (*complete_flush)(void *);
    void (*try_flush)(void *);
};
typedef struct mpi_dest_s mpi_dest[1];
typedef struct mpi_dest_s *mpi_dest_ptr;

/* helper functions */

void mpi_dest_alloc_one_queue(mpi_dest_ptr C, size_t queue_size)
{
    uint32_t * sq_pool = malloc(2 * queue_size * sizeof(uint32_t));
    C->buf[0] = sq_pool;
    C->buf[1] = sq_pool + queue_size;
    C->size = queue_size;
}

void mpi_dest_free_one_queue(mpi_dest_ptr C)
{
    free(C->buf[0]);
    memset(C, 0, sizeof(mpi_dest));
}

void mpi_dest_alloc_several_queues(mpi_dest_ptr * C, int npeers, size_t queue_size)
{
    if (npeers == 0) return;
    uint32_t * sq_pool = malloc(npeers * 2 * queue_size * sizeof(uint32_t));
    for (int i = 0; i < npeers; i++) {
        // no, don't clear !!!
	// memset(C[i], 0, sizeof(mpi_dest));
	C[i]->buf[0] = sq_pool + queue_size * (2 * i);
	C[i]->buf[1] = sq_pool + queue_size * (2 * i + 1);
        C[i]->size = queue_size;
	// C[i]->pending = 0;
	// C[i]->cur = 0;
    }
}

void mpi_dest_free_several_queues(mpi_dest_ptr * C, int npeers)
{
    if (npeers == 0) return;
    free(C[0]->buf[0]);
    for (int i = 0; i < npeers; i++) {
	memset(C[i], 0, sizeof(mpi_dest));
    }
}

/* Note that output channels are used after dispatching. So data is
 * constructed word by word, so there's no need for something smarter
 * (area writes, or splicing).
 */
// message passing interface for communicating bytes.


// this never reschedules a flush, and returns immediately if none is
// pending.
void mpi_dest_complete_flush(mpi_dest M)
{
    time_t t0 = M->pending;
    if (t0 == 0)
        return;

    time_t t1 = time(NULL);
    MPI_Wait(&M->req, MPI_STATUS_IGNORE);
    time_t t = time(NULL);
    M->waited += t - t1;
    M->total_io_wct += t - t0;
    M->pending = 0;
}

int mpi_dest_try_flush(mpi_dest M)
{
    time_t t0 = M->pending;
    if (t0 != 0) {
	/* If there's a request pending, then we'll do the send later. */
	return 0;
    }
    ASSERT_ALWAYS(M->pos <= M->size - M->tailsize);
    uint32_t * tailptr = M->buf[M->cur] + M->size - M->tailsize;
    memcpy(tailptr, &M->pos, sizeof(size_t));
    struct {
        size_t sz;
        size_t over;
    } sz_info[1];
    memcpy(sz_info, tailptr, sizeof(sz_info));
    M->pending = time(NULL);
    // fprintf(stderr, "send [%d] (%zu, over: %d) ->%d\n", M->tag, sz_info->sz, sz_info->over != 0, M->peer);
    /* note that the closing message contains no significant data.
     * However the receiving side does not know it. */
    MPI_Isend(M->buf[M->cur], M->size * sizeof(uint32_t),
            MPI_BYTE, M->peer, M->tag, MPI_COMM_WORLD, &M->req);
    M->cur ^= 1;
    M->pos = 0;
    return 1;
}

void mpi_dest_push(mpi_dest C, uint32_t v)
{
    if (C->pos == C->size - C->tailsize) {
	/* At this point, if a request was pending (thus on the other
	 * buffer), it is time to wait for its completion, and once it
	 * has completed, initiate a second one.
	 */
	mpi_dest_complete_flush(C);
        uint32_t * tailptr = C->buf[C->cur] + C->size - C->tailsize;
        memset(tailptr, 0, C->tailsize * sizeof(uint32_t));
	mpi_dest_try_flush(C);
    }
    ASSERT(C->pos < C->size - C->tailsize);
    C->buf[C->cur][C->pos] = v;
    ++C->b->pos;
    if (++C->pos == C->size - C->tailsize) {
	/* if there's a request pending, we don't need to block yet
	 * before being able to send our payload. Although admittedly it
         * makes very little difference whether we block now or later...*/
        uint32_t * tailptr = C->buf[C->cur] + C->size - C->tailsize;
        memset(tailptr, 0, C->tailsize * sizeof(uint32_t));
	mpi_dest_try_flush(C);
    }
}

int mpi_dest_put(mpi_dest_ptr C, uint32_t * p, size_t n)
{
    if (p == NULL && n == 0) {
        // close
        mpi_dest_complete_flush(C);
        uint32_t * tailptr = C->buf[C->cur] + C->size - C->tailsize;
        memset(tailptr, -1, C->tailsize * sizeof(uint32_t));
        mpi_dest_try_flush(C);
        mpi_dest_complete_flush(C);
        return -1;
    }
    for(size_t i = 0 ; i < n ; i++) {
        mpi_dest_push(C, p[i]);
    }
    return 0;
}


// this does not allocate the buffer for the send queues.
data_dest_ptr mpi_dest_alloc_partial()
{
    mpi_dest_ptr p = malloc(sizeof(mpi_dest));
    memset(p, 0, sizeof(mpi_dest));
    p->b->put = (int(*)(void*,uint32_t*,size_t))&mpi_dest_put;
    // p->b->progress = NULL;
    ASSERT_ALWAYS(sizeof(size_t) % sizeof(uint32_t) == 0);
    p->tailsize = 2*iceildiv(sizeof(size_t),sizeof(uint32_t));

    return (data_dest_ptr) p;
}

void mpi_dest_free(data_dest_ptr M)
{
    // Normally, we expect the buffer to have been allocated
    // elsewhereby e.g. mpi_dest_alloc_queues, so that freeing must have
    // been performed _before_ entering this function. For this reason,
    // we consider it's an error to have anything non-zero in the
    // mpi_dest struct.
    ASSERT_ALWAYS(((mpi_dest_ptr)M)->buf[0] == NULL);
    ASSERT_ALWAYS(((mpi_dest_ptr)M)->buf[1] == NULL);
    // free(((mpi_dest_ptr)M)->buf);
    free(M);
}
/*}}}*/


struct master_data_s {/*{{{*/
    balancing bal;
    const char * mfile;
    uint32_t row_block0;
    uint32_t col_block0;
    uint32_t row_cellbase;
    uint32_t col_cellbase;

    uint32_t *fw_rowperm;
    uint32_t *fw_colperm;

    uint32_t *sent_rows;
    uint32_t *exp_rows;

    parallelizing_info_ptr pi;
};
typedef struct master_data_s master_data[1];
typedef struct master_data_s *master_data_ptr;/*}}}*/

int who_has_row_bare(master_data m, uint32_t rnum)/*{{{*/
{
    int b;
    ASSERT(rnum < m->bal->trows);
    if (rnum < m->row_block0) {
	b = rnum / (m->row_cellbase + 1);
    } else {
	int q = m->bal->trows % m->bal->h->nh;
	b = (rnum - q) / m->row_cellbase;
    }
    return b;
}

int who_has_row(master_data m, uint32_t rnum)
{
    int b = who_has_row_bare(m, rnum);
    return b * m->bal->h->nv;
}
/*}}}*/
int who_has_col_bare(master_data m, uint32_t cnum) /*{{{*/
{
    int b;
    ASSERT(cnum < m->bal->tcols);
    if (cnum < m->col_block0) {
	b = cnum / (m->col_cellbase + 1);
    } else {
	int q = m->bal->tcols % m->bal->h->nv;
	b = (cnum - q) / m->col_cellbase;
    }
    return b;
}

int who_has_col(master_data m, uint32_t cnum)
{
    int b = who_has_col_bare(m, cnum);
    // int j = b / m->pi->wr[0]->ncores;
    // int t = b % m->pi->wr[0]->ncores;
    // return j * m->pi->m->ncores + t;
    return b;
}
/*}}}*/

void read_bfile(master_data m, const char * bfile)
{
    balancing_init(m->bal);
    balancing_read(m->bal, bfile);

    /* FIXME: this does not really belong here. OTOH the slaves do work
     * with something half-baked, as a matter of fact. */
    m->bal->trows = m->bal->h->nrows;
    m->bal->tcols = m->bal->h->ncols;
    if (m->bal->h->flags & FLAG_PADDING) {
	m->bal->tcols = m->bal->trows = MAX(m->bal->h->nrows, m->bal->h->ncols);
    }
}

void share_bfile_header_data(parallelizing_info_ptr pi, balancing_ptr bal)
{
    global_broadcast(pi->m, bal->h, sizeof(balancing_header), 0, 0);
    global_broadcast(pi->m, &bal->trows, sizeof(uint32_t), 0, 0);
    global_broadcast(pi->m, &bal->tcols, sizeof(uint32_t), 0, 0);
}

/* {{{ slave loops 1 and 2 */

struct slave_data_s {/*{{{*/
    balancing bal;		// only partially filled !!!
    parallelizing_info_ptr pi;
    int my_i, my_j;
    uint32_t my_nrows;
    uint32_t my_ncols;
    uint32_t my_row0;
    uint32_t my_col0;

    // char *tmatfile;		/* This stores the transposed matrix file */
    // char *tmatfile_aux;		/* This stores the transposed matrix file */
    // FILE *auxfile;

    size_t expected_size;       // labelled in data words, not counting
                                // the extra transfer window size.

    // will be dropped at some point.
    // size_t rq_size;
    // uint32_t *rq;
    // uint32_t *rq_ptr;

    /* This is the final result. It's allocated only when we're halfway through.
     */
    matrix_u32 mat;

    thread_pipe_ptr tp;
};
typedef struct slave_data_s slave_data[1];
typedef struct slave_data_s *slave_data_ptr;/*}}}*/

void set_slave_variables(slave_data s, param_list pl, parallelizing_info_ptr pi)/*{{{*/
{
    s->pi = pi;

    share_bfile_header_data(pi, s->bal);

    /* Of course this does not make sense to come here if we have a
     * separate master job.  */
    s->my_j = s->pi->wr[0]->jrank * s->pi->wr[0]->ncores + s->pi->wr[0]->trank;
    s->my_i = s->pi->wr[1]->jrank * s->pi->wr[1]->ncores + s->pi->wr[1]->trank;
    // int gridpos = pi->m->jrank * pi->m->ncores + pi->m->trank;

    uint32_t quo_rows = s->bal->trows / s->bal->h->nh;
    uint32_t quo_cols = s->bal->tcols / s->bal->h->nv;
    int rem_rows = s->bal->trows % s->bal->h->nh;
    int rem_cols = s->bal->tcols % s->bal->h->nv;

    s->my_ncols = quo_cols + (s->my_j < rem_cols);
    s->my_nrows = quo_rows + (s->my_i < rem_rows);
    s->my_row0 = s->my_i * quo_rows + MIN(rem_rows, s->my_i);
    s->my_col0 = s->my_j * quo_cols + MIN(rem_cols, s->my_j);

    /*
    char *suffix;
    int rc = asprintf(&suffix, "tr.h%dv%d", s->my_i, s->my_j);
    ASSERT_ALWAYS(rc >= 0);
    s->tmatfile = build_mat_auxfile(mfile, suffix, ".bin");
    // s->tmatfile_aux = build_mat_auxfile(mfile, suffix, ".bin.aux");
    free(suffix);
    */

    /* Get a rough estimate on the number of coefficients we will see.
     * Since this has to be rewinded, it remains in core memory, thus the
     * necessity of guessing this.
     */
    s->expected_size = s->bal->h->ncoeffs / pi->m->totalsize;
    s->expected_size += 2 * sqrt(s->expected_size);
    s->expected_size += 2 * s->bal->trows;
    s->expected_size += s->expected_size / 10;

    /* For small matrices, the deviation is sometimes quite high. Don't
     * bother for such tiny amounts.
     */
    s->expected_size += 1 << 20;

    // in an mpi environment, the threads will see more data passing
    // through than just the expected size.
    int use_auxfile = 0;
    param_list_parse_int(pl, "balancing_use_auxfile", &use_auxfile);

    size_t queue_size = default_queue_size;
    if (param_list_parse_size_t(pl, "balancing_queue_size", &queue_size)) {
        queue_size /= sizeof(uint32_t);
    }
    s->tp = thread_pipe_alloc(s->expected_size + (pi->m->njobs > 1 ? queue_size : 0), use_auxfile);

    /*
    if (!use_auxfile) {
	s->rq_size = s->my_nrows + s->my_ncols + s->bal->h->ncoeffs / gsize;
	s->rq_size += s->rq_size / 10;
	s->rq_size += default_queue_size;
    } else {
	// we'll use an auxiliary file. In this case, we only care about
	// having the write window reasonably large.
	s->rq_size = disk_io_window_size;
	s->auxfile = fopen(s->tmatfile_aux, "w");
	FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }
    */
    printf("J%uT%u (%2d,%2d) expects"
	   " rows %" PRIu32 "+%" PRIu32 " cols %" PRIu32 "+%" PRIu32 ".\n",
           pi->m->jrank, pi->m->trank,
           s->my_i, s->my_j,
           s->my_row0, s->my_nrows,
           s->my_col0, s->my_ncols);
    /*
    char buf[16];
	   // " Datafile %s."
	   " Allocating %s MB\n",
	   // s->tmatfile,
	   size_disp(s->rq_size * sizeof(uint32_t), buf));
    s->rq = malloc(s->rq_size * sizeof(uint32_t));
    s->rq_ptr = s->rq;
    */
}/*}}}*/

void clear_slave_variables(slave_data s)/*{{{*/
{
    // free(s->rq);
    // s->rq = NULL;
    thread_pipe_free(s->tp);
    // free(s->tmatfile);
    // s->tmatfile = NULL;
    // free(s->tmatfile_aux);
    // s->tmatfile_aux = NULL;
}/*}}}*/

struct slave_dest_s {/*{{{*/
    data_dest b;
    uint32_t *row_weights;
    uint32_t *col_weights;
    uint32_t current_row;
    int incoming_rowindex;
    uint64_t tw;
    slave_data_ptr s;
    /* these two are defined only for the second pass. they serve as a
     * marker for the pass number. */
    uint32_t ** fill;
};
typedef struct slave_dest_s slave_dest[1];
typedef struct slave_dest_s *slave_dest_ptr;

uint32_t slave_dest_put(slave_dest_ptr R, uint32_t * p, size_t n)
{
    for(size_t i = 0 ; i < n ; i++) {
        uint32_t x = p[i];
        if (R->incoming_rowindex) {
            R->current_row = x;
            R->incoming_rowindex = 0;
            // don't increase the row counter upon exit.
            R->b->r++;
            ASSERT(R->current_row <= R->s->bal->trows);
            ASSERT(R->current_row >= R->s->my_row0);
            ASSERT(R->current_row < R->s->my_row0 + R->s->my_nrows);
            continue;
        }
        if (x == UINT32_MAX) {
            R->incoming_rowindex = 1;
            continue;
        }
        ASSERT(x <= R->s->bal->tcols);
        ASSERT(x >= R->s->my_col0);
        ASSERT(x < R->s->my_col0 + R->s->my_ncols);

        uint32_t r = R->current_row - R->s->my_row0;
        uint32_t c = x - R->s->my_col0;
        if (R->fill == NULL) {
            R->row_weights[r]++;
            R->col_weights[c]++;
        } else if (R->s->mat->transpose) {
            *(R->fill[c]++) = r;
        } else {
            *(R->fill[r]++) = c;
        }
        R->tw++;
    }
    R->b->pos += n;
    int rc = 0; // R->current_row == UINT32_MAX ? -1 : 0;
    // fprintf(stderr, "slave_dest_put(%zu) returns %d (last row seen %u)\n", n, rc, R->current_row);
    return rc;
}

void slave_dest_stats(slave_dest_ptr s)
{
    uint32_t * row_weights = s->row_weights;
    uint32_t my_nrows = s->s->my_nrows;
    uint32_t * col_weights = s->col_weights;
    uint32_t my_ncols = s->s->my_ncols;
    // some stats
    double s1 = 0;
    double s2 = 0;
    uint64_t tw = 0;
    uint32_t row_min = UINT32_MAX, row_max = 0;
    double row_avg, row_dev;
    s1 = s2 = 0;
    for (uint32_t i = 0; i < my_nrows; i++) {
	double d = row_weights[i];
	if (row_weights[i] > row_max)
	    row_max = row_weights[i];
	if (row_weights[i] < row_min)
	    row_min = row_weights[i];
	s1 += d;
	s2 += d * d;
	tw += row_weights[i];
    }
    row_avg = s1 / my_nrows;
    row_dev = sqrt(s2 / my_nrows - row_avg * row_avg);

    uint32_t col_min = UINT32_MAX, col_max = 0;
    double col_avg, col_dev;
    s1 = s2 = 0;
    for (uint32_t j = 0; j < my_ncols; j++) {
	double d = col_weights[j];
	if (col_weights[j] > col_max)
	    col_max = col_weights[j];
	if (col_weights[j] < col_min)
	    col_min = col_weights[j];
	s1 += d;
	s2 += d * d;
    }
    col_avg = s1 / my_ncols;
    col_dev = sqrt(s2 / my_ncols - col_avg * col_avg);

    printf("[J%uT%u] N=%" PRIu32 " W=%" PRIu64 " "
	   "R[%" PRIu32 "..%" PRIu32 "],%.1f~%.1f "
	   "C[%" PRIu32 "..%" PRIu32 "],%.1f~%.1f.\n",
	   s->s->pi->m->jrank, s->s->pi->m->trank,
           my_nrows, tw,
	   row_min, row_max, row_avg, row_dev,
	   col_min, col_max, col_avg, col_dev);
}

data_dest_ptr slave_dest_alloc(slave_data s)
{
    slave_dest_ptr x = malloc(sizeof(slave_dest));
    memset(x, 0, sizeof(slave_dest));
    x->b->put = (int(*)(void *, uint32_t *, size_t)) slave_dest_put;
    x->s = s;
    x->current_row = 0;
    x->row_weights = malloc(s->my_nrows * sizeof(uint32_t));
    x->col_weights = malloc(s->my_ncols * sizeof(uint32_t));
    memset(x->row_weights, 0, s->my_nrows * sizeof(uint32_t));
    memset(x->col_weights, 0, s->my_ncols * sizeof(uint32_t));
    return (data_dest_ptr) x;
}

void slave_dest_free(data_dest_ptr xx)
{
    slave_dest_ptr x = (slave_dest_ptr) xx;
    free(x->row_weights);
    free(x->col_weights);
    free(xx);
}

void slave_dest_engage_final_pass(slave_dest_ptr D)
{
    slave_data_ptr s = D->s;
    if (D->s->mat->transpose) {
	s->mat->size = s->my_ncols + D->tw;
	s->mat->p = malloc(s->mat->size * sizeof(uint32_t));
	D->fill = malloc(s->my_ncols * sizeof(uint32_t *));
	uint32_t *lptr = s->mat->p;
	memset(D->fill, 0, s->my_ncols * sizeof(uint32_t));
	for (uint32_t j = 0; j < s->my_ncols; j++) {
	    *lptr++ = D->col_weights[j];
	    D->fill[j] = lptr;
	    lptr += D->col_weights[j];
	}
	ASSERT_ALWAYS((uint32_t) (lptr - s->mat->p) == s->my_ncols + D->tw);
    } else {
	s->mat->size = s->my_nrows + D->tw;
	s->mat->p = malloc(s->mat->size * sizeof(uint32_t));
	D->fill = malloc(s->my_nrows * sizeof(uint32_t *));
	uint32_t *lptr = s->mat->p;
	memset(D->fill, 0, s->my_nrows * sizeof(uint32_t));
	for (uint32_t u = 0; u < s->my_nrows; u++) {
	    *lptr++ = D->row_weights[u];
	    D->fill[u] = lptr;
	    lptr += D->row_weights[u];
	}
	ASSERT_ALWAYS((uint32_t) (lptr - s->mat->p) == s->my_nrows + D->tw);
    }
    D->b->r = 0;
}
/*}}}*/

void slave_loop(slave_data s)
{
    char name[80];
    // parallelizing_info_ptr pi = s->pi;
    // int gridpos = pi->m->jrank * pi->m->ncores + pi->m->trank;
    data_source_ptr input = (data_source_ptr) s->tp->src;

    data_dest_ptr output = slave_dest_alloc(s);
    snprintf(name, sizeof(name),
            "[J%uT%u] slave loop 1", s->pi->m->jrank, s->pi->m->trank);
    mf_pipe(input, output, name);

    thread_source_rewind((thread_source_ptr) input);
    slave_dest_stats((slave_dest_ptr) output);

    slave_dest_engage_final_pass((slave_dest_ptr) output);
    printf("[J%uT%u] building local matrix\n",
            s->pi->m->jrank, s->pi->m->trank);
    snprintf(name, sizeof(name),
            "[J%uT%u] slave loop 2", s->pi->m->jrank, s->pi->m->trank);
    mf_pipe(input, output, name);

    // the thread source should be freed elsewhere.
    // thread_source_free(input);
    slave_dest_free(output);
}

// this one is copied almost verbatim from mf_pipe.
void endpoint_loop(parallelizing_info_ptr pi, param_list_ptr pl, slave_data_ptr * slaves, int nslaves)
{
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;
    time_t t = t0;

    size_t queue_size = default_queue_size;
    if (param_list_parse_size_t(pl, "balancing_queue_size", &queue_size)) {
        queue_size /= sizeof(uint32_t);
    }

    // the master rank is zero.
    data_source_ptr input = mpi_source_alloc(0, queue_size);
    ((mpi_source_ptr)input)->nparallel = pi->m->ncores;

    for (;;) {
        uint32_t * ptr = NULL;  // we are not providing a buffer.
        size_t n = input->get(input, &ptr, 0);
        if (n == 0) break;
        int tag = ((mpi_source_ptr)input)->tag;
        ASSERT_ALWAYS(slaves[tag]->tp);
        data_dest_ptr output = (data_dest_ptr) slaves[tag]->tp->dst;
        ASSERT_ALWAYS(output);
        int r = output->put(output, ptr, n);
	t = time(NULL);
        if (r < 0) {
            fprintf(stderr, "endpoint leaves\n");
            break;
        }
	if (t >= t1) {
	    t1 = t + 1;
            // thread_dest don't return the row count anyway.
            // mf_progress(input, output, t-t0, name);
	}
    }
    // mf_progress(input, output, t-t0, name);
    // printf("\n");
    // we must flush all thread buffers.
    for(int tag = 0 ; tag < nslaves ; tag++) {
        if (!slaves[tag])
            continue;
        ASSERT_ALWAYS(slaves[tag]->tp);
        data_dest_ptr output = (data_dest_ptr) slaves[tag]->tp->dst;
        ASSERT_ALWAYS(output);
        output->put(output, NULL, 0);
    }
}

/* }}} */

/* {{{ master loop */

void set_master_variables(master_data m, parallelizing_info_ptr pi)/*{{{*/
{
    uint32_t quo_r = m->bal->trows / m->bal->h->nh;
    uint32_t quo_c = m->bal->tcols / m->bal->h->nv;
    int rem_r = m->bal->trows % m->bal->h->nh;
    int rem_c = m->bal->tcols % m->bal->h->nv;

    m->pi = pi;

    m->row_cellbase = quo_r;
    m->col_cellbase = quo_c;
    m->row_block0 = (quo_r + 1) * rem_r;
    m->col_block0 = (quo_c + 1) * rem_c;

    /* these two fields are in fact used only for debugging */
    m->sent_rows = malloc(pi->m->totalsize * sizeof(uint32_t));
    m->exp_rows = malloc(pi->m->totalsize * sizeof(uint32_t));

    for (int i = 0; i < (int) pi->wr[1]->totalsize; i++) {
        uint32_t e = quo_r + (i < rem_r);
        for(int j = 0 ; j < (int) pi->wr[0]->totalsize ; j++) {
            m->exp_rows[i*pi->wr[0]->totalsize+j] = e;
            m->sent_rows[i*pi->wr[0]->totalsize+j] = 0;
        }
    }

    if (m->bal->h->flags && FLAG_REPLICATE) {
        uint32_t maxdim = MAX(m->bal->h->nrows, m->bal->h->ncols);
	m->fw_rowperm = m->fw_colperm = malloc(maxdim * sizeof(uint32_t));
	memset(m->fw_rowperm, -1, maxdim * sizeof(uint32_t));
	uint32_t *x = (m->bal->colperm ? m->bal->colperm : m->bal->rowperm);
	for (uint32_t i = 0; i < maxdim; i++) {
	    uint32_t j = x[i];
	    ASSERT(m->fw_colperm[j] == UINT32_MAX);
	    m->fw_colperm[j] = i;
	}
    } else {
	ASSERT_ALWAYS(m->bal->h->flags & FLAG_COLPERM);
	m->fw_colperm = malloc(m->bal->tcols * sizeof(uint32_t));
	memset(m->fw_colperm, -1, m->bal->tcols * sizeof(uint32_t));
	for (uint32_t i = 0; i < m->bal->trows; i++) {
	    ASSERT(m->fw_colperm[m->bal->colperm[i]] == UINT32_MAX);
	    m->fw_colperm[m->bal->colperm[i]] = i;
	}

	ASSERT_ALWAYS(m->bal->h->flags & FLAG_ROWPERM);
	m->fw_rowperm = malloc(m->bal->trows * sizeof(uint32_t));
	memset(m->fw_rowperm, -1, m->bal->trows * sizeof(uint32_t));
	for (uint32_t i = 0; i < m->bal->tcols; i++) {
	    ASSERT(m->fw_rowperm[m->bal->colperm[i]] == UINT32_MAX);
	    m->fw_rowperm[m->bal->rowperm[i]] = i;
	}
    }

    /* one more check. The cost is tiny compared to what we do in other
     * parts of the code. */
    printf("Consistency check...");
    fflush(stdout);
    uint32_t *ttab = malloc(m->bal->h->nh * sizeof(uint32_t));
    memset(ttab, 0, m->bal->h->nh * sizeof(uint32_t));
    for (uint32_t j = 0; j < m->bal->h->nrows; j++) {
	ttab[who_has_row_bare(m, m->fw_rowperm[j])]++;
    }
    ASSERT_ALWAYS(m->bal->h->nh == pi->wr[1]->totalsize);
    for (uint32_t k = 0; k < m->bal->h->nh; k++) {
        ASSERT_ALWAYS(ttab[k] == m->exp_rows[k * m->bal->h->nv]);
    }
    printf("ok\n");
}/*}}}*/

void clear_master_variables(master_data m)/*{{{*/
{
    if (m->bal->h->flags && FLAG_REPLICATE) {
	free(m->fw_rowperm);
    } else {
	free(m->fw_rowperm);
	free(m->fw_colperm);
    }
    m->fw_rowperm = m->fw_colperm = NULL;
    free(m->sent_rows); m->sent_rows = NULL;
    free(m->exp_rows); m->exp_rows = NULL;
}/*}}}*/

struct master_dispatcher_s { /*{{{*/
// the dispatching task accomplished by
// the master job is stateful. The following struct defines the state.
    data_dest b;
    master_data_ptr m;
    uint32_t crow_togo;
    int noderow;
    int npeers;
    data_dest_ptr * x;
};
typedef struct master_dispatcher_s master_dispatcher[1];
typedef struct master_dispatcher_s *master_dispatcher_ptr;

int master_dispatcher_put(master_dispatcher_ptr d, uint32_t * p, size_t n)
{
    master_data_ptr m = d->m;
    uint32_t r = d->b->r;

    if (p == NULL && n == 0) {
        // special case: drain the output pipes (and do it in a blocking way)
        /* terminate everywhere */
        for (int i = 0; i < d->npeers; i++) {
            // printf("closing %d/%d\n", i, d->npeers);
            d->x[i]->put(d->x[i], NULL, 0);
        }

        return 0;
    }

    for (size_t s = 0; s < n; s++) {
        if (d->crow_togo == 0) {
            d->crow_togo = p[s];
            d->noderow = who_has_row(m, m->fw_rowperm[r]);
            /* Send a new row info _only_ to the nodes on the
             * corresponding row in the process grid ! */
            for (int i = 0; i < (int) m->bal->h->nv; i++) {
                uint32_t x[2] = {UINT32_MAX, m->fw_rowperm[r], };
                data_dest_ptr where = d->x[d->noderow + i];
                where->put(where, x, 2);
                /* debug only */
                m->sent_rows[d->noderow + i]++;
                if(m->sent_rows[d->noderow + i] >
                        m->exp_rows[d->noderow + i]) {
                    fprintf(stderr,
                            "d->noderow = %d\n", d->noderow);
                    fprintf(stderr,
                            "i=%d\n", i);
                    fprintf(stderr,
                            "m->sent_rows[d->noderow + i] = %"PRIu32"\n",
                            m->sent_rows[d->noderow + i]);
                    fprintf(stderr,
                            "m->exp_rows[d->noderow + i] = %"PRIu32"\n",
                            m->exp_rows[d->noderow + i]);
                }
                ASSERT_ALWAYS(m->sent_rows[d->noderow + i] <=
                        m->exp_rows[d->noderow + i]);
            }
            /* This advances, even if the row is not complete. */
            r++;
        } else {
            uint32_t c = m->fw_colperm[p[s]];
            int nodecol = who_has_col(m, c);
            ASSERT_ALWAYS(d->noderow + nodecol < (int) d->m->pi->m->totalsize);
            data_dest_ptr where = d->x[d->noderow + nodecol];
            where->put(where, &c, 1);
            d->crow_togo--;
        }
        d->b->pos++;
    }
    d->b->r = r;
    return 0;
}

void master_dispatcher_stats(master_dispatcher_ptr d)
{
    for (int i = 0; i < d->npeers; i++) {
        /* There's a -1 that comes from the end marker (which is in fact
         * the same as a new row marker)
         */
        char buf[16];
        printf("[%d] sent %zu rows, %zu coeffs (%s)\n",
                i, (size_t) d->m->sent_rows[i],
                d->x[i]->pos - 2 * d->m->sent_rows[i],
                size_disp(d->x[i]->pos * sizeof(uint32_t), buf));
    }
}

data_dest_ptr master_dispatcher_alloc(master_data m, parallelizing_info_ptr pi, slave_data_ptr * slaves, size_t queue_size)
{
    master_dispatcher_ptr d = malloc(sizeof(master_dispatcher));
    memset(d, 0, sizeof(master_dispatcher));
    d->b->put = (int(*)(void*,uint32_t *, size_t))master_dispatcher_put;
    d->m = m;
    int n = pi->m->totalsize;
    d->x = malloc(n * sizeof(data_dest_ptr));
    d->npeers = n;

    for(int j0 = 0 ; j0 < (int) pi->wr[0]->njobs ; j0++)
    for(int j1 = 0 ; j1 < (int) pi->wr[1]->njobs ; j1++)
    for(int t0 = 0 ; t0 < (int) pi->wr[0]->ncores ; t0++)
    for(int t1 = 0 ; t1 < (int) pi->wr[1]->ncores ; t1++) {
        int j = j0 + j1 * pi->wr[0]->njobs;
        // int t = t0 + t1 * pi->wr[0]->ncores;
        int i1 = t1 + j1 * pi->wr[1]->ncores;
        int i0 = t0 + j0 * pi->wr[0]->ncores;
        int i = i0 + i1 * pi->wr[0]->totalsize;
        if (j == 0) {
            // threads which are on the same node are special. We _have_ to
            // treat them in a different manner, for two reasons:
            // - there has to be a different control flow for _receiving_ the
            //   message. Unfortunately, by design, peers are not recognized
            //   as mpi jobs so we can't use the message passing interface
            //   for talking to these threads.
            // - this will also be the situation for smallish matrices, where
            //   requiring mpi would be undesirable.
            ASSERT_ALWAYS(slaves[i] != NULL);
            d->x[i] = (data_dest_ptr) slaves[i]->tp->dst;
        } else {
            ASSERT_ALWAYS(slaves[i] == NULL);
            d->x[i] = mpi_dest_alloc_partial();
            ((mpi_dest_ptr) d->x[i])->peer = j;
            ((mpi_dest_ptr) d->x[i])->tag = i;
            mpi_dest_alloc_one_queue((mpi_dest_ptr) d->x[i], queue_size);
        }
    }
    return (data_dest_ptr) d;
}

void master_dispatcher_free(data_dest_ptr xd, parallelizing_info_ptr pi, slave_data_ptr * slaves MAYBE_UNUSED)
{
    master_dispatcher_ptr d = (master_dispatcher_ptr) xd;
    // int n = pi->m->totalsize;
    // the thread_pipe structures are freed in another place.
    for(int j0 = 0 ; j0 < (int) pi->wr[0]->njobs ; j0++)
    for(int j1 = 0 ; j1 < (int) pi->wr[1]->njobs ; j1++)
    for(int t0 = 0 ; t0 < (int) pi->wr[0]->ncores ; t0++)
    for(int t1 = 0 ; t1 < (int) pi->wr[1]->ncores ; t1++) {
        int j = j0 + j1 * pi->wr[0]->njobs;
        // int t = t0 + t1 * pi->wr[0]->ncores;
        int i1 = t1 + j1 * pi->wr[1]->ncores;
        int i0 = t0 + j0 * pi->wr[0]->ncores;
        int i = i0 + i1 * pi->wr[0]->totalsize;
        if (j) {
            mpi_dest_free_one_queue((mpi_dest_ptr) d->x[i]);
            mpi_dest_free(d->x[i]);
        }
    }

    free(d->x);
    free(xd);
}
/*}}}*/

void master_loop_inner(master_data m, data_source_ptr input, data_dest_ptr output)/*{{{*/
{
    mf_pipe(input, output, "main");
    uint32_t r = output->r;

    ASSERT(r == m->bal->h->nrows);
    printf("Master loop finished\n");
    /* complete from nrows to maxdim if necessary ! */
    for (; r < m->bal->trows;) {
        uint32_t zero = 0;
        r += output->put(output, &zero, 1);
    }
    r += output->put(output, NULL, 0);
    printf("Final sends completed, everybody should finish.\n");
}

void master_loop(master_data m, parallelizing_info_ptr pi, param_list_ptr pl, slave_data_ptr * slaves, int nslaves MAYBE_UNUSED)
{
    size_t esz = (m->bal->h->ncoeffs + m->bal->trows) * sizeof(uint32_t);
    size_t queue_size = default_queue_size;
    if (param_list_parse_size_t(pl, "balancing_queue_size", &queue_size)) {
        queue_size /= sizeof(uint32_t);
    }
    data_source_ptr input;
#ifdef  HAVE_CURL
    int is_curl;
    if ((is_curl = (strstr(m->mfile, "://") != NULL)))
        input = curl_source_alloc(m->mfile, esz);
    else
#endif  /* HAVE_CURL */
    input = file_source_alloc(m->mfile, esz);
    data_dest_ptr output = master_dispatcher_alloc(m, pi, slaves, queue_size);
    master_loop_inner(m, input, output);
    master_dispatcher_stats((master_dispatcher_ptr) output);
    master_dispatcher_free(output, pi, slaves);
#ifdef  HAVE_CURL
    if (is_curl)
        curl_source_free(input);
    else
#endif
    file_source_free(input);
}

/*}}}*/

/* }}} */

struct do_loops_arg {
    parallelizing_info_ptr pi;
    param_list_ptr pl;
    slave_data_ptr * slaves;
    master_data_ptr m;
    int leader;
};

int gridpos(parallelizing_info_ptr pi)
{
    // everything in this program is stored in the row-major format which
    // correspond to the matrix's splitting.
    int j0 = pi->wr[0]->jrank;
    int j1 = pi->wr[1]->jrank;
    int t0 = pi->wr[0]->trank;
    int t1 = pi->wr[1]->trank;
    int i1 = t1 + j1 * pi->wr[1]->ncores;
    int i0 = t0 + j0 * pi->wr[0]->ncores;
    return i0 + i1 * pi->wr[0]->totalsize;
}

void * do_loops(struct do_loops_arg * arg)
{
    if (arg->leader) {
        ASSERT_ALWAYS(arg->pi->m->trank == 0);
        if (arg->pi->m->jrank == 0) {
            master_loop(arg->m, arg->pi, arg->pl, arg->slaves, arg->pi->m->totalsize);
        } else {
            endpoint_loop(arg->pi, arg->pl, arg->slaves, arg->pi->m->totalsize);
        }
    } else {
        parallelizing_info_ptr pi = arg->pi;
        slave_loop(arg->slaves[gridpos(pi)]);
    }
    return NULL;
}

void * balancing_get_matrix_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg)
{
    master_data m;
    slave_data_ptr * slaves;

    memset(m, 0, sizeof(master_data));

    m->mfile = arg->mfile;

    if (pi->m->trank == 0) {
        slaves = malloc(pi->m->totalsize * sizeof(slave_data_ptr));
        memset(slaves, 0, pi->m->totalsize * sizeof(slave_data_ptr));
    }
    thread_broadcast(pi->m, (void**) &slaves, 0);

    slave_data_ptr s = malloc(sizeof(slave_data));
    pthread_mutex_lock(pi->m->th->m);
    ASSERT_ALWAYS(slaves[gridpos(pi)] == NULL);
    slaves[gridpos(pi)] = s;
    pthread_mutex_unlock(pi->m->th->m);
    memset(s, 0, sizeof(slave_data));
    serialize(pi->m);

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        read_bfile(m, arg->bfile);

        int ok = 1;

        ok = ok && m->bal->h->nh == pi->wr[1]->ncores * pi->wr[1]->njobs;
        ok = ok && m->bal->h->nv == pi->wr[0]->ncores * pi->wr[0]->njobs;

        if (!ok) {
            fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                    "%ux%u splitting from %s\n",
                    pi->wr[1]->njobs,
                    pi->wr[0]->njobs,
                    pi->wr[1]->ncores,
                    pi->wr[0]->ncores,
                    m->bal->h->nh, m->bal->h->nv, arg->bfile);
            exit(1);
        }

	printf("Total %" PRIu32 " rows %" PRIu32 " cols %" PRIu64 " coeffs\n",
	       m->bal->h->nrows, m->bal->h->ncols, m->bal->h->ncoeffs);
	if (m->bal->h->flags & FLAG_PADDING) {
	    uint32_t maxdim = MAX(m->bal->h->nrows, m->bal->h->ncols);
	    printf("Padding to a matrix of size %" PRIu32 "x%" PRIu32 "\n",
		   maxdim, maxdim);
	}
        set_master_variables(m, pi);
        memcpy(s->bal->h, m->bal->h, sizeof(balancing_header));
        s->bal->trows = m->bal->trows;
        s->bal->tcols = m->bal->tcols;
    }

    set_slave_variables(s, pl, pi);
    memcpy(s->mat, arg, sizeof(matrix_u32));

    serialize(pi->m);

    //////////////////////////////////////////////////////////////////

    struct do_loops_arg s_arg[1] = {
        { .pi = pi, .pl = pl, .m = m, .slaves = slaves, .leader = 0, },
    };
    struct do_loops_arg l_arg[1] = {
        { .pi = pi, .pl = pl, .m = m, .slaves = slaves, .leader = 1, },
    };

    pthread_t leader_thread;

    if (pi->m->trank == 0)
        pthread_create(&leader_thread, NULL, (void*(*)(void*))do_loops, l_arg);

    do_loops(s_arg);

    if (pi->m->trank == 0)
        pthread_join(leader_thread, NULL);

    /* master must flush all pending send queues */

    serialize(pi->m);

    /* each local node has a matrix_u32 structure in s->mat ; copy this to a
     * safe place. */
    memcpy(arg, s->mat, sizeof(matrix_u32));
    memset(s->mat, 0, sizeof(matrix_u32));

    clear_slave_variables(s);
    clear_master_variables(m);

    return arg;
}
