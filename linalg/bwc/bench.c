#include "cado.h"

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
 *
 * It must be used with the INNER matrix, not the .info one. So either use
 * the vanilla cado format (before balance), or one of the .h<i>.v<j>
 * matrices as output by balance.
 */

#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <inttypes.h>
#include <unistd.h>
#include "bwc_config.h"
#include "matmul.h"
#include "abase.h"
#include "macros.h"
#include "params.h"
#include "worker-threads.h"
#include "utils.h"
// #include "debug.h"

void usage()
{
    fprintf(stderr,
            "Usage: ./bench [--impl <implementation>] [--tmax <time>] [--nmax <n_iter>] [--nchecks <number> | --nocheck] [-r|--rebuild] [-t|--transpose] [--nthreads <number>] [--cycles <frequency>] -- <file0> [<file1> ... ]\n");
    exit(1);
}


struct private_args {
    matmul_t mm;
    abt * src;
    abt * dst;
};

struct bench_args {
    abobj_t xx;
    param_list pl;
    const char * impl;
    int nthreads;
    int transpose;
    int nchecks;
    int rebuild;
    double freq;
    char ** mfiles;
    struct private_args * p;
    struct worker_threads_group * tg;
};

void init_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    int t0 = time(NULL);
    struct private_args * p = ba->p + tnum;

    p->mm = matmul_init(ba->xx, 0, 0, ba->mfiles[tnum], ba->impl, ba->pl, !ba->transpose);

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
        matmul_build_cache(p->mm, NULL);
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

    p->src = abinit(ba->xx, nc);
    p->dst = abinit(ba->xx, nr);

    abzero(ba->xx, p->dst, nr);
    abzero(ba->xx, p->src, nc);
}/*}}}*/

void check_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;

    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    unsigned int nr0 = nr;
    unsigned int nc0 = nc;
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nr);
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nc);

    abt * dstT = abinit(ba->xx, nc);
    abt * srcT = abinit(ba->xx, nr);
    abzero(ba->xx, dstT, nc);
    abzero(ba->xx, srcT, nr);

    abt * checkA = abinit(ba->xx, abnbits(ba->xx));
    abt * checkB = abinit(ba->xx, abnbits(ba->xx));

    abcopy(ba->xx, dstT, p->src, nc);
    abcopy(ba->xx, srcT, p->dst, nr);
    printf("T%d src(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) p->src, abbytes(ba->xx, nc0) / sizeof(unsigned long)));
    matmul_mul(p->mm, p->dst, p->src, 1);
    printf("T%d dst(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) p->dst, abbytes(ba->xx, nr0) / sizeof(unsigned long)));
    // debug_write(p->dst, abbytes(ba->xx, nr), "/tmp/Lmul");

    abdotprod(ba->xx, checkA, p->dst, srcT, nr0);

    printf("T%d srcT(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) srcT, abbytes(ba->xx, nr0) / sizeof(unsigned long)));
    matmul_mul(p->mm, dstT, srcT, 0);
    printf("T%d dstT(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) dstT, abbytes(ba->xx, nc0) / sizeof(unsigned long)));

    // debug_write(dstT, abbytes(ba->xx, nc), "/tmp/Rmul");

    abdotprod(ba->xx, checkB, p->src, dstT, nc0);

    if (memcmp(checkA, checkB, aboffset(ba->xx, abnbits(ba->xx)) * sizeof(abt)) != 0) {
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d : Check failed\n", tnum);
        pthread_mutex_unlock(&tg->mu);
        abort();
    }

    abclear(ba->xx, checkA, abnbits(ba->xx));
    abclear(ba->xx, checkB, abnbits(ba->xx));

    abclear(ba->xx, dstT, nc);
    abclear(ba->xx, srcT, nr);

    matmul_aux(p->mm, MATMUL_AUX_ZERO_STATS);
}/*}}}*/

void mul_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    matmul_mul(p->mm, ba->transpose ? p->src : p->dst,
            ba->transpose ? p->dst : p->src,
            !ba->transpose);
}/*}}}*/

void clear_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    pthread_mutex_lock(&tg->mu);
    matmul_report(p->mm, ba->freq);
    pthread_mutex_unlock(&tg->mu);
    matmul_clear(p->mm);
    abclear(ba->xx, p->dst, nr);
    abclear(ba->xx, p->src, nc);
}/*}}}*/

void banner(int argc, char * argv[])
{
    /* print command line */
    fprintf (stderr, "# (%s) %s", CADO_REV, (argv)[0]);
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

    banner(argc, argv);

    abobj_init(ba->xx);

    ba->impl = "bucket";
    char * file = NULL;
    int nocheck = 0;
    ba->nchecks = 4;
    ba->nthreads = 1;
    ba->freq = 1;

    /* {{{ */
    param_list_init(ba->pl);
    argv++,argc--;
    param_list_configure_knob(ba->pl, "--transpose", &ba->transpose);
    param_list_configure_knob(ba->pl, "--rebuild", &ba->rebuild);
    param_list_configure_knob(ba->pl, "--nocheck", &nocheck);
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

    double tmax = 100.0;
    param_list_parse_double(ba->pl, "tmax", &tmax);
    param_list_parse_int(ba->pl, "nchecks", &ba->nchecks);
    if (nocheck) ba->nchecks = 0;

    const char * tmp;
    if ((tmp = param_list_lookup_string(ba->pl, "impl")) != NULL) {
        ba->impl = tmp;
    }

    unsigned int nbys;

    if (param_list_parse_uint(ba->pl, "nbys", &nbys)) {
        /* The api mandates that we set the desired value for nbys. Here,
         * that's not really our intent, since we really want to bench
         * the layer in its favorite working context. Most of the time,
         * setting nbys is pointless.
         */
        abobj_set_nbys(ba->xx,nbys);
    }

    if (param_list_warn_unused(ba->pl)) {
        usage();
    }/*}}}*/

    ba->p = malloc(ba->nthreads * sizeof(struct private_args));
    ba->tg = worker_threads_init(ba->nthreads);
    fprintf(stderr, "Using implementation \"%s\"\n", ba->impl);
    worker_threads_do(ba->tg, (worker_func_t) &init_func, ba);

    /* {{{ display some info */
    uint64_t ncoeffs_total = 0;
    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        fprintf (stderr, "T%d %s: %u rows %u cols %" PRIu64 " coeffs\n",
                tnum, ba->mfiles[tnum], p->mm->dim[0], p->mm->dim[1], p->mm->ncoeffs);
        ncoeffs_total += p->mm->ncoeffs;
    }
    fprintf (stderr, "total %" PRIu64 " coeffs\n", ncoeffs_total);
    /* }}} */

    setup_seeding(1);
    /* {{{ Do checks if requested */
    for(int t = 0; t < ba->nchecks ; t++) {
        /* create deterministic test values */
        for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
            struct private_args * p = ba->p + tnum;
            abrandom(ba->xx, p->src, p->mm->dim[1]);
            abrandom(ba->xx, p->dst, p->mm->dim[0]);
        }
        worker_threads_do(ba->tg, (worker_func_t) &check_func, ba);
        fprintf(stderr, "Check %d ok\n", t);
    }
    if (ba->nchecks)
        printf("All %d checks passed\n", ba->nchecks);
    /* }}} */

    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        abrandom(ba->xx, p->src, p->mm->dim[1]);
        abrandom(ba->xx, p->dst, p->mm->dim[0]);
    }

#define NLAST   10
    double last[10]={0,};
    double sum_last = 0;

    clock_t t0 = clock();
    double next = 0.25 * CLOCKS_PER_SEC;
    double t1 = 0;
    double t, dt;
    if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
    unsigned int n;
    for(n = 0 ; n < nmax ; n++ ) {
        worker_threads_do(ba->tg, (worker_func_t) &mul_func, ba);
        t = clock();
        dt = t - t0;
        sum_last -= last[n % NLAST];
        last[n % NLAST] = (t - t1) / CLOCKS_PER_SEC;
        sum_last += (t - t1) / CLOCKS_PER_SEC;
        t1 = t;
        if (dt > next || n == nmax - 1) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
            dt /= CLOCKS_PER_SEC;
            printf("%d iters in %2.fs, %.2f/1, %.2f %s (last %u : %.2f/1, %.2f %s)        \r",
                    n, dt, dt/n, ba->freq * 1.0e9 * dt/n/ncoeffs_total, unit,
                    NLAST, sum_last/NLAST, ba->freq * 1.0e9 *sum_last/NLAST/ncoeffs_total, unit);
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

    return 0;
}
