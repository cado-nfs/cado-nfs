#include "cado.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <unistd.h>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include "portability.h"
#include "parallelizing_info.h"
#include "select_mpi.h"
#include "utils.h"
#include "balancing.h"
#include "balancing_workhorse.h"
#include "subdivision.hpp"

/* The entry point of this code is balancing_get_matrix_u32 ; called in a
 * parallel context, it fills the provided matrix_u32 parameter with the
 * sub-matrix that is relevant for the required balancing (also passed
 * from the matrix_u32_ptr parameter.
 *
 */
void balancing_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "sanity_check_vector",
            "while dispatching the matrix, store a fixed matrix times vector product in the given file");
}

void balancing_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "sanity_check_vector");
}

bool has_mpi_thread_multiple()
{
    int prov;
    MPI_Query_thread(&prov);
    return prov >= MPI_THREAD_MULTIPLE;
}

struct dispatcher {/*{{{*/
    parallelizing_info_ptr pi;
    std::string mfile;
    std::string bfile;
    std::string check_vector_filename;
    int withcoeffs;
    /* one matrix_u32_ptr per thread of this endpoint */
    matrix_u32_ptr * args_per_thread;
    std::vector<int> is_reader_map;
    std::vector<int> readers_index;
    int nreaders;
    balancing bal;
    unsigned int nhjobs;
    unsigned int nvjobs;
    uint32_t rows_chunk_big;
    uint32_t cols_chunk_big;
    uint32_t rows_chunk_small;
    uint32_t cols_chunk_small;

    MPI_Comm reader_comm;

    bool is_reader() const { return is_reader_map[pi->m->jrank]; }
    bool is_reader(int i) const { return is_reader_map[i]; }

    template<typename T> T integrate(std::vector<T> & v) const/*{{{*/
    {
        T s = 0;
        for(T & x : v) {
            T y = x;
            x = s;
            s += y;
        }
        return s;
    }/*}}}*/
    dispatcher(parallelizing_info_ptr pi,/*{{{*/
            param_list_ptr pl,
            matrix_u32_ptr * args_per_thread)
        : pi(pi)
          , mfile(args_per_thread[0]->mfile)
          , bfile(args_per_thread[0]->bfile)
          , check_vector_filename(param_list_lookup_string(pl, "sanity_check_vector"))
          , withcoeffs(args_per_thread[0]->withcoeffs)
          , args_per_thread(args_per_thread)
          , is_reader_map(pi->m->njobs, 0)
    {
        // Assume we are reading an N-rows matrix.  Assume we have n0*n1
        // nodes, t0*t1 threads.
        //
        // Determine the number of readers by finding out which of the n0*n1
        // nodes have read access from the matrix.

        is_reader_map[pi->m->jrank] = access(mfile.c_str(), R_OK) == 0;
        MPI_Allgather(MPI_IN_PLACE, 0, 0, &is_reader_map[0], 1, MPI_INT, pi->m->pals);
        readers_index = is_reader_map;
        nreaders = integrate(readers_index);

        if (pi->m->jrank == 0) {
            printf("Beginning balancing with %d readers for file %s\n",
                    nreaders, mfile.c_str());
            for(unsigned int i = 0 ; i < pi->m->njobs ; i++) {
                if (is_reader(i))
                    printf("Job %d is reader number %d\n", i, readers_index[i]);
            }
        }
        ASSERT_ALWAYS(nreaders);
        MPI_Comm_split(MPI_COMM_WORLD, is_reader(), pi->m->jrank, &reader_comm);

        balancing_init(bal);

        /* Reading the _full_ bfile is not useful. Readers do need the
         * colperm. As for the row perm, we only need a fraction of it,
         * but cuting hair might well end up costing more memory (if we
         * have FLAG_REPLICATE for instance).
         */
        if (pi->m->jrank == 0)
            balancing_read_header(bal, bfile.c_str());
        MPI_Bcast(bal, sizeof(balancing), MPI_BYTE, 0, pi->m->pals);

        nhjobs = pi->wr[1]->njobs;
        nvjobs = pi->wr[0]->njobs;
        rows_chunk_big = bal->trows / nhjobs;
        cols_chunk_big = bal->tcols / nvjobs;
        rows_chunk_small = bal->trows / bal->h->nh;
        cols_chunk_small = bal->tcols / bal->h->nv;

    }/*}}}*/
    ~dispatcher() {/*{{{*/
        MPI_Comm_free(&reader_comm);
        balancing_clear(bal);
    }/*}}}*/

    void main();

    int pass_number;

    /* reader stuff */
    std::vector<size_t> offset_per_reader;
    std::vector<uint32_t> fw_rowperm; // fragmented among all readers.
    std::vector<uint32_t> fw_colperm; // identical at all nodes
    void reader_compute_offsets();
    void reader_fill_index_maps();
    void reader_thread();

    /* MPI send from readers to endpoints */
    std::vector<MPI_Request> outstanding;
    std::vector<int> indices;
    // std::vector<MPI_Status> statuses;
    void post_send(std::vector<uint32_t> &, unsigned int);
    void progress(bool wait = false);
    void post_semaphore(unsigned int k);

    /* endpoint stuff */

    /* used in pass 1 for the weights per row */
    std::vector<std::vector<uint32_t>> thread_row_weights;

    /* used on pass 2 to give the row-beginning pointers. */
    std::vector<std::vector<size_t>> thread_row_positions;

    void prepare_pass();
    std::mutex incoming_mutex;
    void endpoint_handle_incoming(std::vector<uint32_t> & Q);
    void endpoint_thread();
};/*}}}*/

/* balancing_get_matrix_u32 -- This is our entry point. {{{
 * arg is thread private
 * We read:
 *      arg->mfile
 *      arg->bfile
 *      arg->withcoeffs
 * On output, we set:
 *      arg->size
 *      arg->p
 * Finally arg->transpose is currently ignored, maybe it's a bug.
 */
void balancing_get_matrix_u32(parallelizing_info_ptr pi, param_list pl,
        matrix_u32_ptr arg)
{
    // REQUIRED: arg->mfile      -- URLs no longer supported.
    ASSERT_ALWAYS(arg->mfile);
    ASSERT_ALWAYS(arg->bfile);
    ASSERT_ALWAYS(!arg->p);
    ASSERT_ALWAYS(!arg->size);

    matrix_u32_ptr * args_per_thread = (matrix_u32_ptr *) shared_malloc_set_zero(pi->m, pi->m->ncores * sizeof(matrix_u32_ptr));
    args_per_thread[pi->m->trank] = arg;
    serialize_threads(pi->m);

    if (pi->m->trank == 0) {
        dispatcher D(pi, pl, args_per_thread);
        D.main();
    }

    serialize_threads(pi->m);
}
/* }}} */

void dispatcher::post_send(std::vector<uint32_t> & Q, unsigned int k)/*{{{*/
{
    if (k == pi->m->jrank) {
        endpoint_handle_incoming(Q);
        return;
    }
    MPI_Request req;
    MPI_Isend(&Q[0], Q.size(), CADO_MPI_UINT32_T, k, 0, pi->m->pals, &req);
    outstanding.push_back(req);
}/*}}}*/

void dispatcher::post_semaphore(unsigned int k)/*{{{*/
{
    if (k != pi->m->jrank) {
        /* we might as well do it in a blocking way */
        uint32_t z = UINT32_MAX;
        MPI_Send(&z, 1, CADO_MPI_UINT32_T, k, 0, pi->m->pals);
    }
}/*}}}*/

void dispatcher::progress(bool wait)/*{{{*/
{
    int n_in, n_out;
    n_in = n_out = outstanding.size();
    if (!n_in) return;
    if (wait) {
        MPI_Waitall(n_in, &outstanding[0], MPI_STATUSES_IGNORE);
        return;
    }
    indices.assign(n_in, 0);
    // statuses.assign(n_in, 0);
    MPI_Testsome(n_in, &outstanding[0],
            &n_out, &indices[0],
            MPI_STATUSES_IGNORE);
    // &statuses[0]);
    int j = 0;
    for(int i = 0 ; i < n_out ; i++) {
        ASSERT_ALWAYS(indices[i] >= j);
        outstanding[j++]=outstanding[indices[i]];
    }
    outstanding.erase(outstanding.begin() + j, outstanding.end());
}/*}}}*/

void dispatcher::reader_compute_offsets()/*{{{*/
{
    subdivision readers_rows(bal->h->nrows, nreaders);
    unsigned int ridx = readers_index[pi->m->jrank];
    unsigned int row0 = readers_rows.nth_block_start(ridx);
    unsigned int row1 = readers_rows.nth_block_end(ridx);

    // Let R == nreaders.
    // All R nodes read from the rw file and deduce the byte size of the
    // orginal submatrix that has 1/R-th of the rows.
    std::string rwfile;
    {
        char * tmp = derived_filename(mfile.c_str(), "rw", ".bin");
        rwfile = tmp;
        free(tmp);
    }
    bool can_read_rw = access(rwfile.c_str(), R_OK) == 0;
    ASSERT_ALWAYS(!is_reader() || can_read_rw);
    FILE * frw = fopen(rwfile.c_str(), "rb");
    ASSERT_ALWAYS(frw);
    fseek(frw, row0 * sizeof(uint32_t), SEEK_SET);
    std::vector<uint32_t> rw(row1-row0,0);
    int r = fread(&rw[0], sizeof(uint32_t), row1-row0, frw);
    ASSERT_ALWAYS(r == (int) (row1 - row0));
    fclose(frw);
    std::vector<size_t> bytes_per_reader(pi->m->njobs, 0);
    for(unsigned int i = row0 ; i < row1 ; i++) {
        size_t coeff_size = (1 + withcoeffs) * sizeof(uint32_t);
        bytes_per_reader[pi->m->jrank] += (1 + rw[i - row0]) * coeff_size;
    }

    // This data is then allgathered into an array of R integers. Each
    // node thus determines at which offset it should read from the main
    // matrix.
    //
    MPI_Allgather(MPI_IN_PLACE, 0, 0,
            &bytes_per_reader[0], 1, CADO_MPI_SIZE_T,
            pi->m->pals);

    offset_per_reader = bytes_per_reader;
    integrate(offset_per_reader);
}/*}}}*/

void dispatcher::reader_thread()/*{{{*/
{
    if (!is_reader()) return;

    subdivision readers_rows(bal->h->nrows, nreaders);
    unsigned int ridx = readers_index[pi->m->jrank];
    unsigned int row0 = readers_rows.nth_block_start(ridx);
    unsigned int row1 = readers_rows.nth_block_end(ridx);

    FILE * f = fopen(mfile.c_str(), "rb");
    ASSERT_ALWAYS(f);
    int rc = fseek(f, offset_per_reader[pi->m->jrank], SEEK_SET);
    ASSERT_ALWAYS(rc == 0);

    std::vector<uint32_t> row;
    std::vector<std::vector<uint32_t>> noderows(nvjobs);
    std::vector<std::vector<uint32_t>> queues(pi->m->njobs);

    size_t queue_size_per_peer = 1 << 20;

    std::vector<uint64_t> check_vector;
    if (pass_number == 2 && !check_vector_filename.empty())
        check_vector.assign(row1 - row0, 0);

    for(unsigned int i = row0 ; i < row1 ; i++) {
        uint32_t rr = fw_rowperm[i - row0];
        // Readers read full lines from the matrix.
        uint32_t w;
        rc = fread(&w, sizeof(uint32_t), 1, f);
        ASSERT_ALWAYS(rc == 1);
        row.assign(w * (1 + withcoeffs), 0);
        int ww = w * (1 + withcoeffs);
        rc = fread(&row[0], sizeof(uint32_t), ww, f);
        ASSERT_ALWAYS(rc == ww);
        // Column indices are transformed.
        for(int j = 0 ; j < ww ; j += 1 + withcoeffs) {
            uint32_t c = row[j];
            uint32_t cc = fw_colperm[c];
            noderows[cc / cols_chunk_big].push_back(cc);
            if (pass_number == 2 && !check_vector_filename.empty())
                check_vector[i - row0] ^= DUMMY_VECTOR_COORD_VALUE(c);
        }

        // Queues to all nodes are filled
        for(unsigned int k = 0 ; k < nvjobs ; k++) {
            auto & C = noderows[k];
            if (C.empty()) continue;

            auto & Q = queues[rr / rows_chunk_big * nvjobs + k];
            Q.push_back(rr);
            Q.push_back(C.size());
            Q.insert(Q.end(), C.begin(), C.end());
            /* very important */
            C.clear();

            // and sends are posted every once in a while.
            if (Q.size() >= queue_size_per_peer)
                post_send(Q, rr / rows_chunk_big * nvjobs + k);
        }
        progress();
    }
    for(unsigned int kk = 0 ; kk < pi->m->njobs ; kk++) {
        auto & Q = queues[kk];
        if (!Q.empty())
            post_send(Q, kk);
    }
    progress(true);
    for(unsigned int kk = 0 ; kk < pi->m->njobs ; kk++) {
        post_semaphore(kk);
    }
    fclose(f);
    if (pass_number == 2 && !check_vector_filename.empty()) {
        /* Allocate a full vector on the leader node */
        std::vector<uint64_t> full;
        full.assign(bal->trows, 0);
        std::vector<int> sizes(nreaders, 0);
        std::vector<int> displs(nreaders, 0);
        for(int i = 0, d = 0 ; i < nreaders ; i++) {
            if (!is_reader(i)) continue;
            unsigned int ridx = readers_index[i];
            sizes[ridx] = readers_rows.nth_block_size(ridx);
            displs[ridx] = d;
            d += sizes[ridx];
        }
        std::copy(check_vector.begin(), check_vector.end(),
                full.begin() + displs[readers_index[pi->m->jrank]]);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                &full[0], &sizes[0], &displs[0],
                CADO_MPI_UINT64_T, reader_comm);

        if (readers_index[pi->m->jrank] == 0) {
            FILE * f = fopen(check_vector_filename.c_str(), "wb");
            int rc = fwrite(&full[0], sizeof(uint64_t), bal->h->nrows, f);
            ASSERT_ALWAYS(rc == (int) bal->h->nrows);
            fclose(f);
        }
    }
}/*}}}*/

void dispatcher::reader_fill_index_maps()/*{{{*/
{
    balancing xbal;
    balancing_init(xbal);
    balancing_read(xbal, bfile.c_str());
    uint32_t quo_r = xbal->trows / xbal->h->nh;
    ASSERT_ALWAYS(xbal->trows % xbal->h->nh == 0);

    fw_colperm.assign(xbal->tcols, -1);
    fw_rowperm.assign(xbal->trows, -1);

    if (xbal->h->flags & FLAG_REPLICATE) {
        uint32_t * xc = xbal->colperm;
        uint32_t * xr = xbal->rowperm;
        ASSERT_ALWAYS(xbal->tcols == xbal->trows);
        /* currently we seem to be supporting this only in case the
         * column permutation is authoritative. */
        if (!xc) {
            fprintf(stderr, "The current code expects a column permutation replicated on rows, not the converse. There is little adaptation work, but yet to be done. Maybe you could pass \"--reorder columns\" to mf_bal ?\n");
            abort();
        }
        ASSERT_ALWAYS(xc);
        ASSERT_ALWAYS(!xr);
        if (!xc) xc = xr;
        if (!xr) xr = xc;
        for (uint32_t i = 0; i < xbal->tcols; i++) {
            ASSERT_ALWAYS(xc[i] < xbal->tcols);
            uint32_t q = balancing_pre_unshuffle(bal, xc[i]);
            ASSERT_ALWAYS(fw_colperm[q] == UINT32_MAX);
            fw_colperm[q] = i;
        }
        /* In this case we arrange so that the replicated permutation is so
         * that eventually, we are still computing iterates of a matrix
         * which is conjugate to the one we're interested in */

        uint32_t nh = xbal->h->nh;
        uint32_t nv = xbal->h->nv;
        ASSERT_ALWAYS(xbal->trows % (nh * nv) == 0);
        uint32_t elem = xbal->trows / (nh * nv);
        uint32_t ix = 0;
        uint32_t iy = 0;
        for(uint32_t i = 0 ; i < nh ; i++) {
            for(uint32_t j = 0 ; j < nv ; j++) {
                ix = (i * nv + j) * elem;
                iy = (j * nh + i) * elem;
                for(uint32_t k = 0 ; k < elem ; k++) {
                    ASSERT(fw_rowperm[xr[iy+k]] == UINT32_MAX);
                    fw_rowperm[xr[iy+k]] = ix+k;
                }
            }
        }
    } else {
        /* In this case, because the row and column permutations depend
         * on the splitting, if we happen to be in block Wiedemann
         * context, and compute iterates of the form M^i, then we will
         * compute something which depends on the splitting. There is no
         * way we can reconcile what we're doing in a consistent way.
         * Therefore we don't bother trying to undo the effect of the
         * shuffled product.
         */

        uint32_t * xc = xbal->colperm;
        for (uint32_t i = 0; i < xbal->tcols; i++) {
            uint32_t j = xc ? xc[i] : i;
            ASSERT_ALWAYS(j < xbal->tcols);
            ASSERT_ALWAYS(fw_colperm[j] == UINT32_MAX);
            fw_colperm[j] = i;
        }

        uint32_t * xr = xbal->rowperm;
        for (uint32_t i = 0; i < xbal->trows; i++) {
            uint32_t j = xr ? xr[i] : i;
            ASSERT_ALWAYS(j < xbal->trows);
            ASSERT(fw_rowperm[j] == UINT32_MAX);
            fw_rowperm[j] = i;
        }
    }
    subdivision readers_rows(bal->h->nrows, nreaders);
    unsigned int row0 = readers_rows.nth_block_start(pi->m->jrank);
    unsigned int row1 = readers_rows.nth_block_end(pi->m->jrank);
    fw_rowperm.erase(fw_rowperm.begin() + row1, fw_rowperm.end());
    fw_rowperm.erase(fw_rowperm.begin(), fw_rowperm.begin() + row0);

    /* one more check. The cost is tiny compared to what we do in other
     * parts of the code. */
    printf("Consistency check ...");
    fflush(stdout);
    std::vector<uint32_t> ttab(xbal->h->nh, 0);
    for (uint32_t j = 0; j < xbal->trows; j++) {
        ttab[fw_rowperm[j] / rows_chunk_small]++;
    }
    ASSERT_ALWAYS(xbal->h->nh == pi->wr[1]->totalsize);
    for (uint32_t k = 0; k < xbal->h->nh; k++) {
        ASSERT_ALWAYS(ttab[k] == quo_r);
    }
    printf(" ok\n");

#if 0
    char buf[16];
    size_t sz = xbal->h->ncols * sizeof(*colmap);
    printf("Creating column map of size %s ...", size_disp(sz, buf));
    fflush(stdout);
    colmap = malloc(sz);
    for(uint32_t i = 0 ; i < xbal->h->ncols ; i++) {
        uint32_t q = balancing_pre_shuffle(xbal, i);
        ASSERT_ALWAYS(q < xbal->h->ncols);
        uint32_t c = fw_colperm[q];
        colmap[i].w = DUMMY_VECTOR_COORD_VALUE(q);
        colmap[i].c = fw_colperm[q];
        colmap[i].n = who_has_col(m, c);
    }
    printf(" done\n");
#endif
#if 0
    free(m->colmap);
#endif

    balancing_clear(xbal);
}/*}}}*/

void dispatcher::prepare_pass()/*{{{*/
{
    for(unsigned int i = 0 ; i < pi->m->ncores ; i++) {
        ASSERT_ALWAYS(args_per_thread[i]->p == NULL);
        ASSERT_ALWAYS(args_per_thread[i]->size == 0);
    }
    if (pass_number == 1) {
        // Each node allocates a local row weight info for each of its
        // theads.
        decltype(thread_row_weights)::value_type v(rows_chunk_small, 0);
        thread_row_weights = { pi->m->ncores, v };
    } else if (pass_number == 2) {
        for(unsigned int i = 0 ; i < pi->m->ncores ; i++) {
            auto & C = thread_row_weights[i];

            decltype(thread_row_positions)::value_type Cint(C.begin(), C.end());

            /* take into account the subrow length before taking the
             * integral */
            if (withcoeffs)
                for(auto & x : Cint) x = 1 + x * 2;
            else
                for(auto & x : Cint) x++;
            uint32_t w = integrate(Cint);

            args_per_thread[i]->size = w;
            args_per_thread[i]->p = (uint32_t *) malloc(w * sizeof(uint32_t));
            memset(args_per_thread[i]->p, 0, w * sizeof(uint32_t));
            /* place markers for the row weights.  */
            for(unsigned int j = 0 ; j < C.size() ; j++) {
                args_per_thread[i]->p[Cint[j]] = C[j];
                Cint[j]++;
            }

            thread_row_positions.emplace_back(std::move(Cint));
            C.clear();
        }
        thread_row_weights.clear();
    }
}/*}}}*/

void dispatcher::endpoint_handle_incoming(std::vector<uint32_t> & Q)/*{{{*/
{
    std::lock_guard<std::mutex> dummy(incoming_mutex);

    for(auto next = Q.begin() ; next != Q.end() ; ) {
        uint32_t rr = *next++;
        uint32_t rs = *next++;
        if (pass_number == 2) {
            if (!withcoeffs) {
                std::sort(next, next + rs);
            } else {
                /* ugly */
                struct cv {
                    uint32_t c;
                    int32_t v;
                    bool operator<(cv const & a) const { return c < a.c; }
                };
                cv * Q0 = (cv *) &next[0];
                cv * Q1 = (cv *) &next[2*rs];
                std::sort(Q0, Q1);
            }
        }

        unsigned int n_row_groups = pi->wr[1]->ncores;
        unsigned int n_col_groups = pi->wr[0]->ncores;
        unsigned int row_group = (rr / rows_chunk_small) % n_row_groups;
        unsigned int row_index = rr % rows_chunk_small;
        if (pass_number == 1) {
            for(unsigned int j = 0 ; j < rs ; j ++) {
                uint32_t c = *next++;
                unsigned int col_group = (c / cols_chunk_small) % n_col_groups;
                unsigned int group = row_group * n_col_groups + col_group;
                thread_row_weights[group][row_index]++;
                if (withcoeffs) next++;
            }
        } else if (pass_number == 2) {
            std::vector<uint32_t *> pointers;
            pointers.reserve(n_col_groups);
            for(unsigned int i = 0 ; i < n_col_groups ; i++) {
                unsigned int col_group = i;
                uint32_t group = row_group * n_col_groups + col_group;
                uint32_t * matrix = args_per_thread[group]->p;
                size_t pos0 = thread_row_positions[group][row_index];
                uint32_t * p0 = matrix + pos0;
                pointers.push_back(p0);
            }
            for(unsigned int j = 0 ; j < rs ; j ++) {
                uint32_t c = *next++;
                unsigned int col_group = (c / cols_chunk_small) % n_col_groups;
                uint32_t * & p = pointers[col_group];
                ASSERT(*p == 0);
                *p++ = c % cols_chunk_small;
                if (withcoeffs) {
                    ASSERT(*p == 0);
                    *p++ = *next++;
                }
            }
            for(unsigned int i = 0 ; i < n_col_groups ; i++) {
                unsigned int col_group = i;
                uint32_t group = row_group * n_col_groups + col_group;
                uint32_t * matrix = args_per_thread[group]->p;
                size_t pos0 = thread_row_positions[group][row_index];
                uint32_t * p0 = matrix + pos0;
                /* verify consistency with the first pass */
                ASSERT_ALWAYS((pointers[col_group] - p0) == (1 + withcoeffs) * p0[-1]);
            }
        }
    }
}/*}}}*/

void dispatcher::endpoint_thread()/*{{{*/
{
    int Qs;
    std::vector<uint32_t> Q;

    if (pi->m->njobs == 1) return;

    for(;;) {
        // On pass 1, endpoint threads do Recv from any source, and update the
        // local row weight for all threads.
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, 0, pi->m->pals, &status);
        MPI_Get_count(&status, CADO_MPI_UINT32_T, &Qs);
        Q.assign(Qs, 0);
        MPI_Recv(&Q[0], Qs, CADO_MPI_UINT32_T, status.MPI_SOURCE, 0, pi->m->pals, MPI_STATUS_IGNORE);
        if (Qs == 1 && Q[0] == UINT32_MAX)
            break;

        endpoint_handle_incoming(Q);
    }
}/*}}}*/


//
// Once this is all done, the reading threads adjust pointers for the
// row beginnings of each local matrix.
//
// A second pass is then done.
//
// The hard question is how we deal with MPI. What do we post in a
// blocking way, what do we post in a non-blockiing way, etc. It
// seems non trivial.
//

void dispatcher::main() {
    if (is_reader()) {
        reader_compute_offsets();
        // Each reader reads a fragment of N/R row index transformations.
        // Each reader reads the full set if N col index transformations.
        reader_fill_index_maps();
    }

    for(pass_number = 1 ; pass_number <= 2 ; pass_number++) {
        prepare_pass();
        // Each reader spawns a writing (to MPI) thread (reading from disk).
        // (non-readers exit immediately)
        std::thread reader([this] { reader_thread(); });

        // Each node spawns a reading (from MPI) thread.
        std::thread endpoint([this] { endpoint_thread(); });

        reader.join();
        endpoint.join();
    }
}
