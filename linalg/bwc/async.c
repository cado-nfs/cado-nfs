#include "cado.h"
#include <stddef.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "bwc_config.h"
#include "async.h"
#include "rusage.h"
#include "macros.h"
#include "bw-common.h"
#include "utils.h"

// int int_caught = 0;
int hup_caught = 0;

static void extract_interval(timing_interval_data * since_last_reset, timing_interval_data * since_beginning, struct timing_data * t)
{
    memcpy(since_last_reset, t->since_last_reset, 2 * sizeof(struct timing_interval_data_s));
    memcpy(since_beginning, t->since_beginning, 2 * sizeof(struct timing_interval_data_s));
    double d[2];
    seconds_user_sys(d);
    since_last_reset[t->which]->job[0] += d[0];
    since_last_reset[t->which]->job[1] += d[1];
    since_beginning[t->which]->job[0] += d[0];
    since_beginning[t->which]->job[1] += d[1];
    thread_seconds_user_sys(d);
    since_last_reset[t->which]->thread[0] += d[0];
    since_last_reset[t->which]->thread[1] += d[1];
    since_beginning[t->which]->thread[0] += d[0];
    since_beginning[t->which]->thread[1] += d[1];
    double w = wct_seconds();
    since_last_reset[t->which]->wct += w;
    since_beginning[t->which]->wct += w;
}

void timing_flip_timer(struct timing_data * t)
{
    double d[2];
    seconds_user_sys(d);
    t->since_last_reset[t->which]->job[0] += d[0];
    t->since_last_reset[t->which]->job[1] += d[1];
    t->since_beginning[t->which]->job[0] += d[0];
    t->since_beginning[t->which]->job[1] += d[1];
    t->since_last_reset[t->which^1]->job[0] -= d[0];
    t->since_last_reset[t->which^1]->job[1] -= d[1];
    t->since_beginning[t->which^1]->job[0] -= d[0];
    t->since_beginning[t->which^1]->job[1] -= d[1];
    thread_seconds_user_sys(d);
    t->since_last_reset[t->which]->thread[0] += d[0];
    t->since_last_reset[t->which]->thread[1] += d[1];
    t->since_beginning[t->which]->thread[0] += d[0];
    t->since_beginning[t->which]->thread[1] += d[1];
    t->since_last_reset[t->which^1]->thread[0] -= d[0];
    t->since_last_reset[t->which^1]->thread[1] -= d[1];
    t->since_beginning[t->which^1]->thread[0] -= d[0];
    t->since_beginning[t->which^1]->thread[1] -= d[1];
    double w = wct_seconds();
    t->since_last_reset[t->which]->wct += w;
    t->since_last_reset[t->which^1]->wct -= w;
    t->since_beginning[t->which]->wct += w;
    t->since_beginning[t->which^1]->wct -= w;
    t->which ^= 1;
}

static void timing_partial_init(struct timing_data * t, int iter)
{
    memset(t, 0, sizeof(struct timing_data));
    t->go_mark = iter;
    t->last_print = iter;
    t->next_print = iter + 1;
    t->next_async_check = iter + 1;
    t->async_check_period = 1;

    double d[2];
    seconds_user_sys(d);
    t->since_last_reset[t->which]->job[0] = -d[0];
    t->since_last_reset[t->which]->job[1] = -d[1];
    t->since_last_reset[t->which^1]->job[0] = 0;
    t->since_last_reset[t->which^1]->job[1] = 0;
    thread_seconds_user_sys(d);
    t->since_last_reset[t->which]->thread[0] = -d[0];
    t->since_last_reset[t->which]->thread[1] = -d[1];
    t->since_last_reset[t->which^1]->thread[0] = 0;
    t->since_last_reset[t->which^1]->thread[1] = 0;
    double w = wct_seconds();
    t->since_last_reset[t->which]->wct = -w;
    t->since_last_reset[t->which^1]->wct = 0;
}

void timing_init(struct timing_data * t, int start, int end)
{
    timing_partial_init(t, start);
    memcpy(t->since_beginning, t->since_last_reset, sizeof(t->since_last_reset));
    t->begin_mark = start;
    t->end_mark = end;
}

void timing_clear(struct timing_data * t MAYBE_UNUSED)
{
}

static void timing_rare_checks(pi_wiring_ptr wr, struct timing_data * t, int iter, int print)
{
    /* We've decided that it was time to check for asynchronous data.
     * Since it's an expensive operation, the whole point is to avoid
     * doing this check too often. */
    // timing_update_ticks(t, iter);

    timing_interval_data since_last_reset[2];
    timing_interval_data since_beginning[2];
    extract_interval(since_last_reset, since_beginning, t);

    /* First, re-evaluate the async checking period */
    double av;
    av = (since_last_reset[0]->wct + since_last_reset[1]->wct) / (iter - t->go_mark);

#ifndef MPI_LIBRARY_MT_CAPABLE
    /* Other threads might still be lingering in the matrix
     * multiplication routines. That's really not good, since we're also
     * going to do mpi ourselves ! */
    serialize_threads(wr);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    double good_period = PREFERRED_ASYNC_LAG / av;
    int guess = 1 + (int) good_period;
    global_broadcast(wr, &guess, sizeof(int), 0, 0);
    /* negative stuff is most probably caused by overflows in a fast
     * context */
    if (guess <= 0) guess = 10;
    if (iter == t->go_mark + 1 && print) {
        printf("Checking for asynchronous events every %d iterations\n", guess);
    }
    t->next_async_check += guess;
    t->async_check_period = guess;

    /* Now do some possibly expensive checks */

    /* This read is unsafe. We need to protect it. */
    unsigned int caught_something = hup_caught; // || int_caught;

    /* propagate the unsafe read. */
    if (wr->trank == 0) {
        int err = MPI_Allreduce(MPI_IN_PLACE, &caught_something, 1,
                MPI_UNSIGNED, MPI_MAX, wr->pals);
        ASSERT_ALWAYS(!err);
    }

    /* reconcile threads */
    // serialize_threads(wr);      // bug in thread_broadcast
    void * ptr = &caught_something;
    thread_broadcast(wr, &ptr, 0);
    caught_something = * (unsigned int *) ptr;

    /* ok, got it. Before we possibly leave, make sure everybody has read
     * the data from the pointer. */
    serialize(wr);

    if (!caught_something) {
        /* And the result of all this is... nothing :-(( */
        return;
    }

    /* We have something ! */
    if (print) {
        printf("Signal caught (iteration %d), restarting timings\n", iter);
    }
    serialize(wr);
    hup_caught = 0;
    // int_caught = 0;
    timing_partial_init(t, iter);
}


void timing_check(parallelizing_info pi, struct timing_data * timing, int iter, int print)
{
    pi_wiring_ptr wr = pi->m;

    // printf("timing_check %d timing%d\n", iter, wr->trank);

    if (iter == timing->next_async_check) {
        // printf("timing_check async %d timing%d\n", iter, wr->trank);
        timing_rare_checks(wr, timing, iter, print);
    }

    if (iter < timing->next_print)
        return;

    /* We're printing something, so we might as well check for signals
     * now */
    timing_rare_checks(wr, timing, iter, print);

    const int period[] = {1,5,20,100,1000,10000,100000,1000000,0};
    int i;
    for(i = 0 ; ; i++) {
        if (period[i] == 0 || iter % period[i])
            break;
    }
    i--;
    ASSERT(iter % period[i] == 0);
    timing->next_print = iter + period[i];

    // timing_update_ticks(timing, iter);

    char buf[20];

    timing_interval_data since_last_reset[2];
    timing_interval_data since_beginning[2];
    extract_interval(since_last_reset, since_beginning, timing);
    double di = iter - timing->go_mark;

    // (avg wct cpu)(cpu % cpu) + (avg comm)(cpu % comm)
    snprintf(buf, sizeof(buf), "%.2f@%.0f%%+%.2f@%.0f%%",
        since_last_reset[1]->wct / di,
        100.0 * since_last_reset[1]->thread[0] / since_last_reset[1]->wct,
        since_last_reset[0]->wct / di,
        100.0 * since_last_reset[0]->thread[0] / since_last_reset[0]->wct);

    if (print)
        printf("iteration %d\n", iter);

    grid_print(pi, buf, sizeof(buf), print);
}

/* stage=0 for krylov, 1 for mksol */
void timing_disp_collective_oneline(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print, int stage)
{
    timing_interval_data since_last_reset[2];
    timing_interval_data since_beginning[2];
    extract_interval(since_last_reset, since_beginning, timing);
    double di = iter - timing->go_mark;

    double dwct0 = since_last_reset[0]->wct;
    double dwct1 = since_last_reset[1]->wct;
    double dwct = dwct0 + dwct1;

    int ndoubles = sizeof(since_last_reset) / sizeof(double);

    // dt must be collected.
    int err;
    double ncoeffs_d = ncoeffs;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(pi->m) {
        err = MPI_Allreduce(MPI_IN_PLACE, since_last_reset, ndoubles, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);

        err = MPI_Allreduce(MPI_IN_PLACE, since_beginning, ndoubles, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);

        err = MPI_Allreduce(MPI_IN_PLACE, &ncoeffs_d, 1, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* dt_recent[] contains the JOB seconds. So we don't have to sum it up over
     * threads. However, we do for ncoeffs ! */
    double ncoeffs_total = 0;
    {
        void * ptr = &ncoeffs_total;
        thread_broadcast(pi->m, &ptr, 0);
        double * main_ncoeffs_total = ptr;
        my_pthread_mutex_lock(pi->m->th->m);
        * main_ncoeffs_total += ncoeffs_d;
        my_pthread_mutex_unlock(pi->m->th->m);
        serialize_threads(pi->m);
        ncoeffs_d = * main_ncoeffs_total;
    }

    /* wct intervals are not aggregated over threads either, so we must
     * do it by hand */
    double aggr_dwct[2] = {0,0};
    {
        void * ptr = aggr_dwct;
        thread_broadcast(pi->m, &ptr, 0);
        double * px = ptr;
        my_pthread_mutex_lock(pi->m->th->m);
        px[0] += since_last_reset[0]->wct;
        px[1] += since_last_reset[1]->wct;
        my_pthread_mutex_unlock(pi->m->th->m);
        serialize_threads(pi->m);
        since_last_reset[0]->wct = px[0];
        since_last_reset[1]->wct = px[1];
    }

    if (print) {
        {
            double puser = since_last_reset[1]->job[0] / dwct1;
            double psys = since_last_reset[1]->job[1] / dwct1;
            double pidle = aggr_dwct[1] / dwct1 - puser - psys;
            double avwct = dwct1 / di;
            double nsc = avwct / ncoeffs_d * 1.0e9;
            char * what_wct = "s";
            if (avwct < 0.1) { what_wct = "ms"; avwct *= 1000.0; }

            printf("N=%d ; CPU: %.2f [%.0f%%cpu, %.0f%%sys, %.0f%% idle]"
                    ", %.2f %s/iter"
                    ", %.2f ns/coeff"
                    "\n",
                    iter,
                    since_beginning[1]->job[0] + since_beginning[1]->job[1],
                    100.0 * puser,
                    100.0 * psys,
                    100.0 * pidle,
                    avwct, what_wct, nsc);
        }

        {
            double puser = since_last_reset[0]->job[0] / dwct0;
            double psys = since_last_reset[0]->job[1] / dwct0;
            double pidle = aggr_dwct[0] / dwct0 - puser - psys;
            double avwct = dwct0 / di;
            char * what_wct = "s";
            if (avwct < 0.1) { what_wct = "ms"; avwct *= 1000.0; }

            printf("N=%d ; COMM: %.2f [%.0f%%cpu, %.0f%%sys, %.0f%% idle]"
                    ", %.2f %s/iter"
                    "\n",
                    iter,
                    since_beginning[0]->job[0] + since_beginning[0]->job[1],
                    100.0 * puser,
                    100.0 * psys,
                    100.0 * pidle,
                    avwct, what_wct);
        }

        {
            /* still to go: timing->end_mark - iter */
            time_t now[1];
            time_t eta[1];
            char eta_string[32] = "not available yet\n";
            time(now);
            *eta = *now + (timing->end_mark - iter) * dwct / di;
            if (di)
                ctime_r(eta, eta_string);

            unsigned int s = strlen(eta_string);
            for( ; s && isspace(eta_string[s-1]) ; eta_string[--s]='\0') ;

            printf("%s: N=%d ; ETA (N=%d): %s [%.2f s/iter]\n",
                   (stage == 0) ? "krylov" : "mksol",
                   iter, timing->end_mark, eta_string, dwct / di);
        }
    }
}


void block_control_signals()
{
    /* Only the master thread receives control signals */
    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    // sigaddset(sset, SIGINT);
    my_pthread_sigmask(SIG_BLOCK, sset, NULL);
}

void sighandler(int sig)
{
    int caught = 0;
    if (sig == SIGHUP) hup_caught = caught = 1;
    // if (sig == SIGINT) int_caught = caught = 1;
    /* Of course, everybody is allowed to print this. */
    if (caught) printf("Signal caught, wait before it is acknowledged\n");
}

void catch_control_signals()
{
    struct sigaction sa[1];
    memset(sa, 0, sizeof(sa));
    sa->sa_handler = sighandler;
    sigaction(SIGHUP, sa, NULL);
    // sigaction(SIGINT, sa, NULL);

    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    // sigaddset(sset, SIGINT);
    my_pthread_sigmask(SIG_UNBLOCK, sset, NULL);
}
