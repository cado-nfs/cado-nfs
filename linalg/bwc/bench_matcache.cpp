#include "cado.h" // IWYU pragma: keep

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
 *
 * It must be used with the INNER matrix, not the .info one. So either use
 * the vanilla cado format (before balance), or one of the .h<i>.v<j>
 * matrices as output by balance.
 */

#include <stdio.h>
#include <stdint.h>             // for uint64_t
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>            // for pthread_mutex_lock, pthread_mutex_unlock
#include <gmp.h>
#include <memory>

#include "raw_matrix_u32.h"     // for matrix_u32
#include "crc.h"        // cado_crc_lfsr
#include "macros.h"
#include "matmul-mf.hpp"
#include "matmul.hpp"
#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "params.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "version_info.h" // cado_revision_string
#include "worker-threads.h"

void usage()
{
    fprintf(stderr,
            "Usage: ./bench [--impl <implementation>] [--tmax <time>] [--nmax <n_iter>] [--nchecks <number> | --nocheck] [-r|--rebuild] [-t|--transpose] [--nthreads <number>] [--cycles <frequency>] -- <file0> [<file1> ... ]\n");
    exit(1);
}


struct private_args {
    matmul_t mm;
    arith_generic::elt * colvec;
    arith_generic::elt * rowvec;
};

struct bench_args {
    arith_generic * xx;
    param_list pl;
    const char * impl;
    int nthreads;
    int transpose;
    int nchecks;
    int rebuild;
    mpz_t prime;
    double freq;
    char ** mfiles;
    const char * source_vec;
    struct private_args * p;
    struct worker_threads_group * tg;
};

void init_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    int t0 = time(NULL);
    struct private_args * p = ba->p + tnum;

    p->mm = matmul_init(ba->xx, 0, 0, ba->mfiles[tnum], ba->impl, ba->pl, !ba->transpose);

    fprintf(stderr, "Expect to find or create cache file %s\n", p->mm->cachefile_name);

    if (!ba->rebuild && matmul_reload_cache(p->mm)) {
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Reusing cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        fprintf(stderr, "T%d Cache load time %ds wct\n",
                tnum, (int) time(NULL) - t0);
        pthread_mutex_unlock(&tg->mu);
    } else {
        clock_t ct0 = clock();
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Building cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        pthread_mutex_unlock(&tg->mu);
        ASSERT_ALWAYS(p->mm->store_transposed == ba->transpose);
        matrix_u32 m;
        int withcoeffs = mpz_cmp_ui(ba->prime, 2) > 0;
        mf_prepare_matrix_u32(p->mm, m, ba->mfiles[tnum], withcoeffs);
        matmul_build_cache(p->mm, m);
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Cache build time %.2fs cpu\n",
                tnum, (double) (clock()-ct0) / CLOCKS_PER_SEC);
        fprintf(stderr, "T%d Saving cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        pthread_mutex_unlock(&tg->mu);
        t0 = time(NULL);
        matmul_save_cache(p->mm);
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Cache save time %ds wct\n",
                tnum, (int) time(NULL) - t0);
        pthread_mutex_unlock(&tg->mu);
        /* XXX ok, it's freakin ugly. We must really rethink this object. */
        if (m->mfile) free((void*) m->mfile);
    }

    pthread_mutex_lock(&tg->mu);
    fprintf(stderr, "T%u uses cache file %s\n",
            tnum,
            /* cache for mmt->locfile, */
            p->mm->cachefile_name);
    pthread_mutex_unlock(&tg->mu);

    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nr);
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nc);

    p->colvec = ba->xx->alloc(nc, ALIGNMENT_ON_ALL_BWC_VECTORS);
    p->rowvec = ba->xx->alloc(nr, ALIGNMENT_ON_ALL_BWC_VECTORS);
    ba->xx->vec_set_zero(p->colvec, nc);
    ba->xx->vec_set_zero(p->rowvec, nr);
}/*}}}*/

void check_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    arith_generic * A = ba->xx;

    std::unique_ptr<arith_cross_generic> AxA(arith_cross_generic::instance(A,A));

    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    unsigned int nr0 = nr;
    unsigned int nc0 = nc;
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nr);
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nc);

    arith_generic::elt * colvec_bis = A->alloc(nc, ALIGNMENT_ON_ALL_BWC_VECTORS);
    arith_generic::elt * rowvec_bis = A->alloc(nr, ALIGNMENT_ON_ALL_BWC_VECTORS);
    A->vec_set(colvec_bis, p->colvec, nc);
    A->vec_set(rowvec_bis, p->rowvec, nr);

    arith_generic::elt * check0 = A->alloc(A->simd_groupsize(), ALIGNMENT_ON_ALL_BWC_VECTORS);
    arith_generic::elt * check1 = A->alloc(A->simd_groupsize(), ALIGNMENT_ON_ALL_BWC_VECTORS);

    /* See the comment in matmul_mul about the direction argument and the
     * number of coordinates of source/destination vectors */

    printf("T%d colvec(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) p->colvec, A->vec_elt_stride(nc0)));
    matmul_mul(p->mm, p->rowvec, p->colvec, 1);
    printf("T%d rowvec(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) p->rowvec, A->vec_elt_stride(nr0)));

    A->vec_set_zero(check0, A->simd_groupsize());
    AxA->add_dotprod(check0, p->rowvec, rowvec_bis, nr0);

    printf("T%d rowvec_bis(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) rowvec_bis, A->vec_elt_stride(nr0)));
    matmul_mul(p->mm, colvec_bis, rowvec_bis, 0);
    printf("T%d colvec_bis(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) colvec_bis, A->vec_elt_stride(nc0)));

    A->vec_set_zero(check1, A->simd_groupsize());
    AxA->add_dotprod(check1, p->colvec, colvec_bis, nc0);

    if (A->vec_cmp(check0, check1, A->simd_groupsize()) != 0) {
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d : Check failed\n", tnum);
        pthread_mutex_unlock(&tg->mu);
        abort();
    }

    A->free(colvec_bis);
    A->free(rowvec_bis);
    A->free(check0);
    A->free(check1);


    matmul_aux(p->mm, MATMUL_AUX_ZERO_STATS);
}/*}}}*/

void mul_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    arith_generic::elt * dstvec = ba->transpose ? p->colvec : p->rowvec;
    arith_generic::elt * srcvec = ba->transpose ? p->rowvec : p->colvec;
    matmul_mul(p->mm, dstvec, srcvec, !ba->transpose);
}/*}}}*/

void clear_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    arith_generic * A = ba->xx;
    unsigned int nr MAYBE_UNUSED = p->mm->dim[0];
    unsigned int nc MAYBE_UNUSED = p->mm->dim[1];
    pthread_mutex_lock(&tg->mu);
    matmul_report(p->mm, ba->freq);
    printf("\n");
    pthread_mutex_unlock(&tg->mu);
    matmul_clear(p->mm);
    A->free(p->colvec); // nc
    A->free(p->rowvec); // nr
}/*}}}*/

void banner(int argc, char * argv[])
{
    /* print command line */
    fprintf (stderr, "# (%s) %s", cado_revision_string, (argv)[0]);
    for (int i = 1; i < (argc); i++)
        fprintf (stderr, " %s", (argv)[i]);
    fprintf (stderr, "\n");

#ifdef  __GNUC__
    fprintf(stderr, "# Compiled with gcc " __VERSION__ "\n");
#endif
    fprintf(stderr, "# Compilation flags " CFLAGS "\n");
}

int main(int argc, char * argv[])
{
    struct bench_args ba[1];

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    banner(argc, argv);

    memset(ba, 0, sizeof(ba));

    ba->impl = "bucket";
    char * file = NULL;
    int nocheck = 0;
    ba->nchecks = 4;
    ba->nthreads = 1;
    ba->freq = 1;
    mpz_init_set_ui(ba->prime, 2);

    /* {{{ */
    param_list_init(ba->pl);
    argv++,argc--;
    param_list_configure_switch(ba->pl, "--transpose", &ba->transpose);
    param_list_configure_switch(ba->pl, "--rebuild", &ba->rebuild);
    param_list_configure_switch(ba->pl, "--nocheck", &nocheck);
    param_list_configure_alias(ba->pl, "--transpose", "-t");
    param_list_configure_alias(ba->pl, "--rebuild", "-r");

    int wild = 0;
    for( ; argc ; ) {
        if (param_list_update_cmdline(ba->pl, &argc, &argv)) { continue; }
        if (strcmp(argv[0],"--") == 0) {
            argv++, argc--;
            break;
        }
        if (argv[0][0] != '-' && wild == 0) {
            file = argv[0];
            argv++,argc--;
            wild++;
            continue;
        }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    unsigned int nmax = UINT_MAX;
    param_list_parse_mpz(ba->pl, "prime", ba->prime);
    param_list_parse_uint(ba->pl, "nmax", &nmax);
    param_list_parse_int(ba->pl, "nthreads", &ba->nthreads);
    const char * unit;
    if (param_list_parse_double(ba->pl, "cycles", &ba->freq)) {
        unit = "cy/c";
    } else {
        unit = "ns/c";
    }

    if (ba->nthreads == 1) {
        if (file) {
            ba->mfiles = &file;
        } else {
            if (argc) {
                ba->mfiles = argv;
            } else {
                fprintf(stderr, "No file !\n");
                exit(1);
            }
        }
    } else {
        ba->mfiles = argv;
        if (file) {
            fprintf(stderr, "When using --nthreads, please specify files after a terminating --\n");
            exit(1);
        } else if (argc != ba->nthreads) {
            fprintf(stderr, "%u threads requested, but %u files given on the command line.\n", ba->nthreads, argc);
            exit(1);
        }
    }

    param_list_lookup_string(ba->pl, "matmul_bucket_methods");
    param_list_lookup_string(ba->pl, "l1_cache_size");
    param_list_lookup_string(ba->pl, "l2_cache_size");

    double tmax = 100.0;
    param_list_parse_double(ba->pl, "tmax", &tmax);
    param_list_parse_int(ba->pl, "nchecks", &ba->nchecks);
    if (nocheck) ba->nchecks = 0;

    const char * tmp;
    if ((tmp = param_list_lookup_string(ba->pl, "impl")) != NULL) {
        ba->impl = tmp;
    }


    /* The api mandates that we set the desired value for nbys. Here,
     * that's not really our intent, since we really want to bench
     * the layer in its favorite working context. Most of the time,
     * setting nbys is pointless.
     */
    unsigned int nbys = 64;
    std::unique_ptr<arith_generic> A(arith_generic::instance(ba->prime, nbys));
    ba->xx = A.get();

    /* may leave nbys at the default value ! */
    param_list_parse_uint(ba->pl, "nbys", &nbys);

    param_list_lookup_string(ba->pl, "srcvec");
    if (param_list_warn_unused(ba->pl)) {
        usage();
    }/*}}}*/

    ba->p = (private_args *) malloc(ba->nthreads * sizeof(struct private_args));
    ba->tg = worker_threads_init(ba->nthreads);
    fprintf(stderr, "Using implementation \"%s\"\n", ba->impl);
    worker_threads_do(ba->tg, (worker_func_t) &init_func, ba);

    /* {{{ display some info */
    uint64_t ncoeffs_total = 0;
    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        fprintf (stderr, "T%d %s: %u rows %u cols %" PRIu64 " coeffs\n",
                tnum, ba->mfiles[tnum],
                p->mm->dim[0],
                p->mm->dim[1],
                p->mm->ncoeffs);
        ncoeffs_total += p->mm->ncoeffs;
    }
    fprintf (stderr, "total %" PRIu64 " coeffs\n", ncoeffs_total);
    /* }}} */

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);

    /* {{{ Do checks if requested */
    for(int t = 0; t < ba->nchecks ; t++) {
        /* create deterministic test values */
        for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
            struct private_args * p = ba->p + tnum;
            A->vec_set_random(p->colvec, p->mm->dim[1], rstate);
            A->vec_set_random(p->rowvec, p->mm->dim[0], rstate);
            /* If we want shared vectors, this is the way to go. */
            /* Note that for such a test, the clear_func must be skipped
             * or improved, since we don't really want to free() the same
             * pointer twice */
            /*
            if (ba->transpose)
                p->rowvec = ba->p[0].rowvec;
            else
                p->colvec = ba->p[0].colvec;
                */
        }
        worker_threads_do(ba->tg, (worker_func_t) &check_func, ba);
        fprintf(stderr, "Check %d ok\n", t);
    }
    if (ba->nchecks)
        printf("All %d checks passed\n", ba->nchecks);
    /* }}} */

    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        A->vec_set_random(p->colvec, p->mm->dim[1], rstate);
        A->vec_set_random(p->rowvec, p->mm->dim[0], rstate);
    }

    if ((tmp = param_list_lookup_string(ba->pl, "srcvec")) != NULL) {
        struct private_args * p = ba->p;
        if (ba->nthreads > 1) {
            fprintf(stderr, "srcvec incompatible with multithread\n");
            exit(1);
        }
        A->vec_set_zero(p->colvec, p->mm->dim[1]);
        A->vec_set_zero(p->rowvec, p->mm->dim[0]);
        FILE * f = fopen(tmp, "rb");
        if (f == NULL) {
            fprintf(stderr, "fopen(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
        /* with ba->transpose, we're doing vector times matrix (rowvec
         * times matrix -> colvec) */
        arith_generic::elt * srcvec = ba->transpose ? p->rowvec : p->colvec;
        arith_generic::elt * dstvec = ba->transpose ? p->colvec : p->rowvec;
        int n = p->mm->dim[1 ^ ba->transpose];
        fprintf(stderr, "reading %zu bytes from %s\n",
                n * sizeof(uint64_t), tmp);
        int nread = fread(srcvec, sizeof(uint64_t), n, f);
        if (nread != n) {
            fprintf(stderr, "short read (%d < %d)\n", nread, n);
            exit(1);
        }
        fclose(f);
        worker_threads_do(ba->tg, (worker_func_t) &mul_func, ba);
        char * dstfile;
        int rc = asprintf(&dstfile, "%s.dst", tmp);
        ASSERT_ALWAYS(rc >= 0);
        f = fopen(dstfile, "wb");
        int nw = p->mm->dim[0 ^ ba->transpose];
        fprintf(stderr, "writing %zu bytes to %s\n",
                nw * sizeof(uint64_t), dstfile);
        int nwritten = fwrite(dstvec, sizeof(uint64_t), nw, f);
        if (nwritten != nw) {
            fprintf(stderr, "short write (%d < %d)\n", nwritten, nw);
            exit(1);
        }
        fclose(f);
        fprintf(stderr, "Saved [%s]%s * [%s] to [%s]\n",
                p->mm->cachefile_name,
                ba->transpose ? "^T" : "",
                tmp,
                dstfile);
        free(dstvec);
    }
    gmp_randclear(rstate);

#define NLAST   10
    double last[10]={0,};
    double sum_last = 0;

    clock_t t0 = clock();
    double next = 0.25 * CLOCKS_PER_SEC;
    double t1 = 0;
    double t, dt;
    if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
    unsigned int n;
    printf("Note: timings in seconds per iterations are given in cpu seconds ;\n");
    printf("Note: with %d threads, this means %.2f walltime seconds.\n", ba->nthreads, 1.0 / ba->nthreads);
    for(n = 0 ; n < nmax ; ) {
        worker_threads_do(ba->tg, (worker_func_t) &mul_func, ba);
        t = clock();
        dt = t - t0;
        sum_last -= last[n % NLAST];
        last[n % NLAST] = (t - t1) / CLOCKS_PER_SEC;
        sum_last += (t - t1) / CLOCKS_PER_SEC;
        t1 = t;
        unsigned int nlast = NLAST;
        n++;
        if (n < nlast) nlast = n;
        if (dt > next || n == nmax - 1) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
            dt /= CLOCKS_PER_SEC;
            printf("%d iters in %2.fs, %.3f/1, %.2f %s (last %u : %.3f/1, %.2f %s)        \r",
                    n, dt, dt/n, ba->freq * 1.0e9 * dt/n/ncoeffs_total, unit,
                    nlast, sum_last/nlast, ba->freq * 1.0e9 *sum_last/nlast/ncoeffs_total, unit);
            fflush(stdout);
            if (dt > tmax)
                break;
        }
    }
    printf("\n");
    // printf("Scanned %lu coeffs in total\n", n * ncoeffs_total);

    worker_threads_do(ba->tg, (worker_func_t) &clear_func, ba);
    worker_threads_clear(ba->tg);
    param_list_clear(ba->pl);
    free(ba->p);
    mpz_clear(ba->prime);

    return 0;
}
