#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cinttypes>
#include <ctime>
#include <cerrno>
#include <cmath>
#include <climits>

#include <memory>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>


#include "mod_ul.h"
#include "balancing.hpp"
#include "rowset_heap.h"
#include "mf_bal.hpp"
#include "portability.h" // strdup // IWYU pragma: keep
#include "timing.h"     // wct_seconds
#include "fix-endianness.h" // fread32_little
#include "macros.h"
#include "params.h"

#include "utils_cxx.hpp"        // for unique_ptr<FILE>

typedef int (*sortfunc_t) (const void *, const void *);

static int revcmp_2u32(const uint32_t * a, const uint32_t * b)
{
    if (a[0] > b[0]) return -1;
    if (b[0] > a[0]) return 1;
    if (a[1] > b[1]) return -1;
    if (b[1] > a[1]) return 1;
    return 0;
}

/* {{{ Datatype containing description for slices of the matrix */
struct slice {
    uint32_t * r;
    uint32_t nrows;
    uint32_t coeffs;
    uint32_t i0;
};

#if 0
static int slice_finding_helper(const uint32_t * key, const struct slice * elem)/*{{{*/
{
    if (*key < elem->i0) {
        return -1;
    } else if (*key >= elem->i0 + elem->nrows) {
        return 1;
    }
    return 0;
}

static unsigned int which_slice(const struct slice * slices, unsigned int ns, uint32_t k)
{
    const struct slice * found = (const struct slice *) bsearch(
            (const void *) &k, (const void *) slices,
            ns, sizeof(struct slice),
            (sortfunc_t) &slice_finding_helper);
    if (found == NULL) {
        return UINT_MAX;
    } else {
        return found - slices;
    }
}
#endif
/*}}}*/

static void free_slices(struct slice * slices, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i++) {
        free(slices[i].r);
    }
    free(slices);
}

static struct slice * alloc_slices(unsigned int water, unsigned int n)
{
    struct slice * res;
    unsigned int i;

    res = (struct slice*) malloc(n * sizeof(struct slice));

    uint32_t const common_size = water / n;

    // we now require equal-sized blocks.
    ASSERT_ALWAYS(water % n == 0);

    for(i = 0 ; i < n ; i++) {
        res[i].nrows = common_size; // min_size + (i < (water % n));
        res[i].r = (uint32_t *) malloc(res[i].nrows * sizeof(uint32_t));
        res[i].coeffs = 0;
    }

    return res;
}
/* }}} */
/* {{{ shuffle_rtable: This is the basic procedure which dispatches rows
 * (or columns) in buckets according to our preferred strategy
 */
static struct slice * shuffle_rtable(
        const char * text,
        uint32_t * rt,
        uint32_t n,
        unsigned int ns)
{
    uint32_t i;
    struct bucket * heap;
    struct slice * slices;
    clock_t t;

    t = -clock();
    qsort(rt, n, 2 * sizeof(uint32_t), (sortfunc_t) &revcmp_2u32);
    t += clock();
    printf("sort time %.1f s\n", (double) t / CLOCKS_PER_SEC);

    
    t = -clock();

    slices = alloc_slices(n, ns);

    heap = (struct bucket *) malloc(ns * sizeof(struct bucket));
    for(i = 0 ; i < ns ; i++) {
        heap[i].s = 0;
        heap[i].i = i;
        heap[i].room = slices[i].nrows;
    }
    make_heap(heap, heap + ns);

    /* Then we're putting each row through the appropriate bucket, using
     * the heap to constantly update.
     * Eventually, we have constructed the set of rows going in the
     * bucket into slices[0]...slices[nslices-1]
     */
    for(i = 0 ; i < n ; i++) {
        pop_heap(heap, heap + ns);
        int const j = heap[ns-1].i;
        int const pos = slices[j].nrows-heap[ns-1].room;
        ASSERT(heap[ns-1].room);
        slices[j].r[pos] = rt[2*i+1];
        heap[ns-1].s += rt[2*i];
        heap[ns-1].room--;
        push_heap(heap, heap + ns);
    }

    t += clock();

    printf("heap fill time %.1f\n", (double) t / CLOCKS_PER_SEC);

    qsort(heap, ns, sizeof(struct bucket), (sortfunc_t) &heap_index_compare);

    for(i = 0 ; i < ns ; i++) {
        int const j = heap[i].i;
        ASSERT(heap[i].i == (int) i);
        printf("%s slice %d, span=%ld, weight=%ld\n",
                text,
                i, slices[j].nrows - heap[i].room,
                heap[i].s);
        slices[j].coeffs = heap[i].s;
        if (heap[i].room != 0) {
            abort();
        }

#if 0
        // This is for giving the possibility of reading the matrix
        // linearly.
        if (sort_positionally) {
            qsort(slices[j].data, slices[j].nrows, sizeof(struct row),
                    (sortfunc_t) row_compare_index);
        }
#endif
    }
    free(heap);

    uint32_t i0 = 0;
    for(i = 0 ; i < ns ; i++) {
        slices[i].i0=i0;
        i0 += slices[i].nrows;
    }

    return slices;
}

/* }}} */

/* {{{ read the mfile header if we have it, deduce #coeffs */
static void read_mfile_header(balancing & bal, std::string const & mfile, int withcoeffs)
{
    struct stat sbuf_mat[1];
    int const rc = stat(mfile.c_str(), sbuf_mat);
    if (rc < 0) {
        fmt::print(stderr, "Reading {}: {} (not fatal)\n", mfile, strerror(errno));
        fmt::print("{}: {} rows {} cols\n",
                mfile, bal.nrows, bal.ncols);
        fmt::print("{}: main input file not present locally,"
                " total weight unknown\n", mfile);
        bal.ncoeffs = 0;
    } else {
        bal.ncoeffs = sbuf_mat->st_size / sizeof(uint32_t) - bal.nrows;
        if (withcoeffs) {
            if (bal.ncoeffs & 1) {
                fprintf(stderr, "Matrix with coefficient must have an even number of 32-bit entries for all (col index, coeff). Here, %" PRIu64 " is odd.\n", bal.ncoeffs);
                abort();
            }
            bal.ncoeffs /= 2;
        }

        int const extra = int(bal.ncols - bal.nrows);
        if (extra > 0) {
            fmt::print("{}: {} rows {} cols ({} extra cols) weight {}\n",
                    mfile, bal.nrows, bal.ncols,
                    extra,
                    bal.ncoeffs);
        } else if (extra < 0) {
            fmt::print("{}: {} rows {} cols ({} extra rows) weight {}\n",
                    mfile, bal.nrows, bal.ncols,
                    -extra,
                    bal.ncoeffs);
        } else {
            fmt::print("{}: {} rows {} cols weight {}\n",
                    mfile, bal.nrows, bal.ncols, bal.ncoeffs);
        }
    }
}
/* }}} */

void mf_bal_decl_usage(cxx_param_list & pl)
{
   param_list_decl_usage(pl, "mfile", "matrix file (can also be given freeform)");
   param_list_decl_usage(pl, "rwfile", "row weight file (defaults to <mfile>.rw)");
   param_list_decl_usage(pl, "cwfile", "col weight file (defaults to <mfile>.cw)");
   param_list_decl_usage(pl, "out", "output file name (defaults to stdout)");
   param_list_decl_usage(pl, "quiet", "be quiet");
   param_list_decl_usage(pl, "rectangular", "accept rectangular matrices (for block Lanczos)");
   param_list_decl_usage(pl, "withcoeffs", "expect a matrix with explicit coefficients (not just 1s)");
   param_list_decl_usage(pl, "rowperm", "permute rows in priority (defaults to auto)");
   param_list_decl_usage(pl, "colperm", "permute columns in priority (defaults to auto)");
   param_list_decl_usage(pl, "skip_decorrelating_permutation", "solve for the matrix M instead of the matrix P*M with P a fixed stirring matrix");
}

void mf_bal_configure_switches(cxx_param_list & pl, struct mf_bal_args * mba)
{
    param_list_configure_switch(pl, "--quiet", &mba->quiet);
    // param_list_configure_switch(pl, "--display-correlation", &display_correlation);
    param_list_configure_switch(pl, "--rectangular", &mba->rectangular);
    param_list_configure_switch(pl, "--withcoeffs", &mba->withcoeffs);
}


void mf_bal_parse_cmdline(struct mf_bal_args * mba, cxx_param_list & pl, int * p_argc, char const *** p_argv)
{
    unsigned int wild =  0;
    (*p_argv)++,(*p_argc)--;
    for(;(*p_argc);) {
        const char * q;
        if (param_list_update_cmdline(pl, p_argc, p_argv)) continue;

        if ((*p_argv)[0][0] != '-' && wild == 0 && (q = strchr((*p_argv)[0],'x')) != NULL) {
            mba->nh = atoi((*p_argv)[0]);
            mba->nv = atoi(q+1);
            wild+=2;
            (*p_argv)++,(*p_argc)--;
            continue;
        }

        if ((*p_argv)[0][0] != '-' && wild == 0) { mba->nh = atoi((*p_argv)[0]); wild++,(*p_argv)++,(*p_argc)--; continue; }
        if ((*p_argv)[0][0] != '-' && wild == 1) { mba->nv = atoi((*p_argv)[0]); wild++,(*p_argv)++,(*p_argc)--; continue; }
        if ((*p_argv)[0][0] != '-' && wild == 2) {
            mba->mfile = (*p_argv)[0];
            wild++;
            (*p_argv)++,(*p_argc)--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", (*p_argv)[0]);
        exit(1);
    }
}

void mf_bal_interpret_parameters(struct mf_bal_args * mba, cxx_param_list & pl)
{
    const char * tmp;

    if (!mba->nh || !mba->nv) {
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }

    if ((tmp = param_list_lookup_string(pl, "reorder")) != NULL) {
        if (strcmp(tmp, "auto") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_AUTO;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_AUTO;
        } else if (strcmp(tmp, "none") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_NO;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_NO;
        } else if (strcmp(tmp, "rows") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_YES;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_NO;
        } else if (strcmp(tmp, "columns") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_NO;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "rows,columns") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_YES;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "columns,rows") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_YES;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "both") == 0) {
            mba->do_perm[0] = mf_bal_args::MF_BAL_PERM_YES;
            mba->do_perm[1] = mf_bal_args::MF_BAL_PERM_YES;
        } else {
            fprintf(stderr, "Argument \"%s\" to the \"reorder\" parameter not understood\n"
                    "Supported values are:\n"
                    "\tauto (default)\n"
                    "\tnone\n"
                    "\trows\n"
                    "\tcolumns\n"
                    "\tboth (equivalent forms: \"rows,columns\" or \"columns,rows\"\n",
                    tmp);
            exit(EXIT_FAILURE);
        }
    }

    param_list_parse(pl, "mfile", mba->mfile);
    param_list_parse(pl, "rwfile", mba->rwfile);
    param_list_parse(pl, "cwfile", mba->cwfile);
    param_list_parse(pl, "out", mba->bfile);
    param_list_parse_int(pl, "skip_decorrelating_permutation", &mba->skip_decorrelating_permutation);
}

void mf_bal_adjust_from_option_string(struct mf_bal_args * mba, const char * opts)
{
    /* Take a comma-separated list of balancing_options, and use that to
     * adjust the mba structure.
     */
    if (!opts) return;
    ASSERT_ALWAYS(opts);

    /* Create a new param_list from opts {{{ */
    auto my_opts = split(opts, ",");
    int n_argc = 1 + int(my_opts.size());
    std::vector<const char *> n_argv_v;
    n_argv_v.reserve(n_argc);
    n_argv_v.push_back("mf_bal");
    for(auto const & s : my_opts) {
        n_argv_v.push_back(s.c_str());
    }
    const char ** n_argv = n_argv_v.data();
    /* }}} */

    cxx_param_list pl2;

    mf_bal_decl_usage(pl2);
    mf_bal_configure_switches(pl2, mba);
    mf_bal_parse_cmdline(mba, pl2, &n_argc, &n_argv);
    mf_bal_interpret_parameters(mba, pl2);
}

void mf_bal(struct mf_bal_args * mba)
{
    int rc;
    /* we sometimes allocate strings, which need to be freed eventually */
    if (!mba->mfile.empty() && mba->rwfile.empty()) {
        mba->rwfile = std::unique_ptr<char, free_delete<char>>(
                derived_filename(mba->mfile.c_str(), "rw", ".bin")).get();
    }

    if (!mba->mfile.empty() && mba->cwfile.empty()) {
        mba->cwfile = std::unique_ptr<char, free_delete<char>>(
                derived_filename(mba->mfile.c_str(), "cw", ".bin")).get();
    }

    if (mba->mfile.empty()) {
        fprintf(stderr, "Matrix file name (mfile) must be given, even though the file itself does not have to be present\n");
        exit(EXIT_FAILURE);
    }
    if (mba->rwfile.empty()) { fprintf(stderr, "No rwfile given\n"); exit(1); }
    if (mba->cwfile.empty()) { fprintf(stderr, "No cwfile given\n"); exit(1); }

    balancing bal;
    balancing_init(bal);

    bal.nh = mba->nh;
    bal.nv = mba->nv;

    struct stat sbuf[2][1];
    rc = stat(mba->rwfile.c_str(), sbuf[0]);
    if (rc < 0) { perror(mba->rwfile.c_str()); exit(1); }
    bal.nrows = sbuf[0]->st_size / sizeof(uint32_t);

    rc = stat(mba->cwfile.c_str(), sbuf[1]);
    if (rc < 0) { perror(mba->cwfile.c_str()); exit(1); }
    bal.ncols = sbuf[1]->st_size / sizeof(uint32_t);

    if (bal.ncols > bal.nrows) {
        fprintf(stderr, "Warning. More columns than rows. There could be bugs.\n");
    }

    read_mfile_header(bal, mba->mfile, mba->withcoeffs);

    /* {{{ Compute the de-correlating permutation */
    if (mba->skip_decorrelating_permutation) {
        /* internal, for debugging. This removes the de-correlating
         * permutation. Nothing to do with what is called
         * "shuffled-product" elsewhere, except that both are taken care
         * of within mf_bal. */
        bal.pshuf[0] = 1;
        bal.pshuf[1] = 0;
        bal.pshuf_inv[0] = 1;
        bal.pshuf_inv[1] = 0;
    } else {
        modulusul_t M;
        modul_initmod_ul(M, MIN(bal.nrows, bal.ncols));
        residueul_t a,b;
        residueul_t ai,bi;
        modul_init(a, M);
        modul_init(b, M);
        modul_init(ai, M);
        modul_init(bi, M);
        modul_set_ul(a, (unsigned long) sqrt(bal.ncols), M);
        modul_set_ul(b, 42, M);

        for( ; modul_inv(ai, a, M) == 0 ; modul_add_ul(a,a,1,M)) ;
        modul_mul(bi, ai, b, M);
        modul_neg(bi, bi, M);

        bal.pshuf[0] = modul_get_ul(a, M);
        bal.pshuf[1] = modul_get_ul(b, M);
        bal.pshuf_inv[0] = modul_get_ul(ai, M);
        bal.pshuf_inv[1] = modul_get_ul(bi, M);

        modul_clear(a, M);
        modul_clear(b, M);
        modul_clear(ai, M);
        modul_clear(bi, M);
        modul_clearmod(M);
    }
    /* }}} */

    const char * text[2] = { "row", "column", };
    unsigned int const matsize[2] = { bal.nrows, bal.ncols, };
    unsigned int const gridsize[2] = { (unsigned int) mba->nh, (unsigned int) mba->nv };
    unsigned int block[2];
    unsigned int padding[2];
    std::unique_ptr<uint32_t[]> weights[2];
    uint32_t nzero[2];
    double avg[2], sdev[2];

    for(int d = 0 ; d < 2 ; d++) {
        /* d == 0 : padding the rows */
        block[d] = iceildiv(matsize[d], mba->nh * mba->nv);
        /* We also want to enforce alignment of the block size with
         * respect to the SIMD things.
         * Given that mmt_vec_init provides 64-byte alignment of vector
         * areas, we may enforce the block size to be a multiple of
         * 8 in order to effectively guarantee 64-byte alignment for all
         * chunks. (Admittedly, this is a bit fragile; if we were to
         * possibly use smaller items, that would change stuff somewhat).
         */
        for ( ; block[d] % MINIMUM_ITEMS_IN_BWC_CHUNKS ; block[d]++);
    }

    if (!mba->rectangular) {
        printf("Padding to a square matrix\n");
        block[0] = block[1] = MAX(block[0], block[1]);
    }

    for(int d = 0 ; d < 2 ; d++) {
        padding[d] = mba->nv * mba->nh * block[d] - matsize[d];
        printf("Padding to %u+%u=%u %ss which is %u blocks of %u*%u=%u %ss\n",
                matsize[d], padding[d], mba->nv * mba->nh * block[d], text[d],
                gridsize[d],
                gridsize[!d], block[d], gridsize[!d] * block[d], text[d]);

        size_t const n = matsize[d] + padding[d];

        /* Read the file with row or column weights */
        std::string const filename = d == 0 ? mba->rwfile : mba->cwfile;
        std::unique_ptr<FILE> const fw(fopen(filename.c_str(), "rb"));
        if (!fw) { perror(filename.c_str()); exit(1); }
        weights[d].reset(new uint32_t[n]);
        std::fill_n(weights[d].get(), n, 0);
        /* Padding rows and cols have zero weight of course */
        double t_w;
        t_w = -wct_seconds();
        size_t const nr = fread32_little(weights[d].get(), matsize[d], fw.get());
        t_w += wct_seconds();

        if (nr < matsize[d]) {
            fmt::print(stderr, "{}: short {} count\n", filename, text[d]);
            exit(EXIT_FAILURE);
        }
        fmt::print("read {} in {:.1f} s ({:.1f} MB / s)\n",
                filename, t_w, 1.0e-6 * double(sbuf[d]->st_size) / t_w);
        /* count zero rows or columns */
        nzero[d] = 0;
        for(size_t r = 0 ; r < matsize[d] ; r++) {
            nzero[d] += (weights[d][r] == 0);
        }

        if (mba->do_perm[d] == mf_bal_args::MF_BAL_PERM_NO) {
            weights[d].reset();
            continue;
        }

        uint32_t ** pperm = d == 0 ? &bal.rowperm : &bal.colperm;
        *pperm = (uint32_t *) malloc(2 * n * sizeof(uint32_t));

        /* prepare for qsort */
        uint32_t * perm = *pperm;
        t_w = -wct_seconds();
        double s1 = 0, s2 = 0;
        uint64_t totalweight = 0;
        for(size_t r = 0 ; r < matsize[d] + padding[d] ; r++) {
            /* Column r in the matrix we work with is actually column rx in
             * the original matrix ! */
            size_t rx = r;
            if (r < matsize[d]) {
                /* The balancing_pre_* function satisfy
                 * f([0,bal.ncols[) \subset [0,bal.ncols[.
                 * and correspond to identity for x>=bal.ncols
                 */
                rx = balancing_pre_unshuffle(bal, r);
                ASSERT(balancing_pre_shuffle(bal, rx) == r);
                ASSERT(rx < matsize[d]);
            }
            uint32_t const w = weights[d][rx];
            totalweight += w;
            double const x = w;
            perm[2*r]=w;
            perm[2*r+1]=r;
            s1 += x;
            s2 += x * x;
        }
        avg[d] = s1 / matsize[d];
        sdev[d] = sqrt(s2 / matsize[d] - avg[d]*avg[d]);
        t_w += wct_seconds();

        if (bal.ncoeffs) {
            if (totalweight != bal.ncoeffs) {
                fmt::print(stderr, "Inconsistency in number of coefficients\n"
                        "From {}: {}, from file sizes; {}\n",
                        filename, totalweight, bal.ncoeffs);
                fmt::print(stderr, "Maybe use the --withcoeffs option for DL matrices ?\n");
                exit(1);
            }
        } else {
            bal.ncoeffs = totalweight;
            fmt::print("{} coefficients counted\n", totalweight);
        }
        fmt::print("{} {} ; avg {:.1f} sdev {:.1f} [scan time {:.1f} s]\n",
                matsize[d], text[d], avg[d], sdev[d], t_w);
    }

    bal.nzrows = nzero[0];
    bal.nzcols = nzero[1];

#if 0
    /* We used to examine the correlation between row and column
     * weight. This is important, as it accounts for some timing jitter
     * on the local matrix products.
     *
     * Unfortunately, in full generality, I'm having difficulties to make
     * sense out of this notion for rectangular matrices, so let's
     * comment it out for the moment.
     */

    if (display_correlation) {
        size_t n = bal.nrows + rpadding;
        FILE * frw = fopen(rwfile, "rb");
        if (frw == NULL) { perror(rwfile); exit(1); }
        uint32_t * rowweights = malloc(n * sizeof(uint32_t));
        memset(rowweights, 0, n * sizeof(uint32_t));
        double t_rw;
        t_rw = -wct_seconds();
        size_t nr = fread32_little(rowweights, bal.nrows, frw);
        t_rw += wct_seconds();
        fclose(frw);
        if (nr < bal.nrows) {
            fprintf(stderr, "%s: short row count\n", rwfile);
            exit(1);
        }
        double rs1 = 0;
        double rs2 = 0;
        double rc_plain = 0;
        double rc_decorr = 0;
        for(size_t r = 0 ; r < n ; r++) {
            double x = rowweights[r];
            rs1 += x;
            rs2 += x * x;
            rc_plain += x * colweights[r];
            rc_decorr += x * bal.colperm[2*r];
        }
        double ravg = rs1 / n;
        double rsdev = sqrt(rs2 / n - ravg*ravg);
        double cavg = s1 / n;
        double csdev = sqrt(s2 / n - cavg*cavg);
        double pcov = rc_plain / n - ravg * cavg;
        double pcorr = pcov / csdev / rsdev;
        double dcov = rc_decorr / n - ravg * cavg;
        double dcorr = dcov / csdev / rsdev;

        printf("%" PRIu32 " rows ; avg %.1f sdev %.1f [scan time %.1f s]\n",
                bal.nrows, ravg, rsdev, t_rw);
        printf("row-column correlation coefficient is %.4f\n",
                pcorr);
        printf("row-column correlation coefficient after decorrelation is %.4f\n",
                dcorr);
    }
#endif

    if (mba->do_perm[0] == mf_bal_args::MF_BAL_PERM_AUTO) {
        int const choose = sdev[1] > sdev[0];
        fmt::print("Choosing a {} permutation based on largest deviation"
                " ({:.2f} > {:.2f})\n",
                text[choose], sdev[choose], sdev[!choose]);
        mba->do_perm[choose] = mf_bal_args::MF_BAL_PERM_YES;
        mba->do_perm[!choose] = mf_bal_args::MF_BAL_PERM_NO;
    }

    bal.flags = 0;

    /* this is where the rowperm / colperm objects change from bein a
     * list of pairs (index, provenance tag) to just the indices
     */
    for(int d = 0 ; d < 2 ; d++) {
        if (mba->do_perm[d] == mf_bal_args::MF_BAL_PERM_NO) continue;
        bal.flags |= (d == 0) ? FLAG_ROWPERM : FLAG_COLPERM;
        uint32_t ** pperm = d == 0 ? &bal.rowperm : &bal.colperm;
        struct slice * h = shuffle_rtable(text[d], *pperm, matsize[d] + padding[d], gridsize[d]);
        /* we can now make the row or column permutation tidier */
        *pperm = (uint32_t *) realloc(*pperm, (matsize[d] + padding[d]) * sizeof(uint32_t));
        for(unsigned int ii = 0 ; ii < gridsize[d] ; ii++) {
            const struct slice * r = &(h[ii]);
            memcpy(*pperm + r->i0, r->r, r->nrows * sizeof(uint32_t));
        }
        free_slices(h, gridsize[d]);
    }

    if (!mba->rectangular) {
        bal.flags |= FLAG_REPLICATE;
    }

    balancing_finalize(bal);

    balancing_write(bal, mba->mfile, mba->bfile);
    balancing_clear(bal);


}

