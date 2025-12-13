#include "cado.h" // IWYU pragma: keep

#include <inttypes.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <gmp.h>
#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

#include "gmp_aux.h"
#include "macros.h"
#include "params.h"
#include "polyselect_arith.h"
#include "polyselect_data_series.h"
#include "polyselect_main_data.h"
#include "polyselect_main_queue.h"
#include "polyselect_poly_header.h"
#include "polyselect_qroots.h"
#include "polyselect_special_q.h"
#include "polyselect_stats.h"
#include "polyselect_thread.h"
#include "polyselect_thread_league.h"
#include "polyselect_thread_team.h"
#include "size_optimization.h"
#include "timing.h"

//#define DEBUG_POLYSELECT
//#define DEBUG_POLYSELECT2



void polyselect_main_data_init_defaults(polyselect_main_data_ptr main)
{
    main->incr = DEFAULT_INCR;
    mpz_init(main->admin);
    mpz_init(main->admax);
    mpz_init(main->N);
    pthread_mutex_init(&main->lock, NULL);
    main->verbose = 0;
    main->sopt_effort = SOPT_DEFAULT_EFFORT;
    main->maxtime = DBL_MAX;
    main->target_E = 0;
    main->keep = DEFAULT_POLYSELECT_KEEP;
    polyselect_stats_init(main->stats, main->keep);
    main->idx = 0;
#ifdef HAVE_HWLOC
    hwloc_topology_init(&main->topology);
    hwloc_topology_load(main->topology);
#endif
    main->nthreads = 0;
    main->finer_grain_threads = 1;
}

void polyselect_main_data_clear(polyselect_main_data_ptr main)
{
    mpz_clear(main->admin);
    mpz_clear(main->admax);
    mpz_clear(main->N);
    pthread_mutex_destroy(&main->lock);
    polyselect_stats_clear(main->stats);
#ifdef HAVE_HWLOC
    hwloc_topology_destroy(main->topology);
#endif
}

/* This returns the bound on i in the algorithm. We choose to set it to
 * the square of the largest prime (this is also a suggestion in the slides)
 * where we consider all primes in [P,2P].
 */
int64_t polyselect_main_data_get_M(polyselect_main_data_srcptr main)
{
    int64_t p = 2 * main->P; // main->Primes[main->lenPrimes - 1];
    return p * p;
}

size_t polyselect_main_data_expected_number_of_pairs(polyselect_main_data_srcptr main)
{
  uint32_t P = main->P; // main->Primes[0];
  int64_t M = polyselect_main_data_get_M(main);
  /* We add 2*M/p^2 entries to the hash table for each p that has
   * roots, and for each root. Since on average we have one root per p,
   * this means that we want the sum of 2M/p^2, for p prime ranging from P
   * to 2P. The sum of 1/(i^2*log(i)) over this same range is
   * 1/log(P)*(1/P-1/(2*P)), hence 1/(2*P*log(P)).
   *
   * The result is therefore 2M/(2*P*log(P)) = M/P/log(P).
   */
  return (size_t) ((uint64_t) M) / (double) P / log(P);
}

/* the number of expected collisions is 8*lenPrimes^2/2/(2P)^2 */
double polyselect_main_data_expected_collisions(polyselect_main_data_srcptr main)
{
#if 0
  uint32_t twoP = Primes[lenPrimes - 1];
  double m = (lenPrimes << 1) / (double) twoP;
  /* we multiply by 0.5 because we discard collisions with f[d] * f[d-2] > 0 */
  return 0.5 * m * m;
#else
  /* make it dependent on how M gets chosen above */
  double d = polyselect_main_data_expected_number_of_pairs(main);
  /* we multiply by 0.5 because we discard collisions with f[d] * f[d-2] > 0 */
  /* Gentle reminder: the domain on which collisions are being search is
   * [-M,+M], and so has size 2M */
  return d * d / 2 / (2*polyselect_main_data_get_M(main)) * 0.5;
#endif
}

/* check that l <= m0/P^2 where l = p1 * p2 * q with P <= p1, p2 <= 2P
   and q is the product of special-q primes.
   This will ensure that we can do rotation by x^(d-3)*g(x), since the
   expected value of a[d-2] is m0/P^2, and x^(d-3)*g(x) has coefficient
   l for degree d-2. */
int polyselect_main_data_check_parameters(polyselect_main_data_srcptr main, mpz_srcptr m0, double q)
{
    double p = 2 * main->P; // main->Primes[main->lenPrimes-1];
    return pow(p, 4) * q < mpz_get_d(m0);
}

/* find suitable values of lq and k, where the special-q part in the degree-1
   coefficient of the linear polynomial is made from k small primes among lq */
unsigned long
find_suitable_lq(polyselect_poly_header_srcptr header,
		 polyselect_qroots_srcptr SQ_R,
                 unsigned long *k,
                 polyselect_main_data_srcptr main)
{
  unsigned long prod = 1;
  unsigned int i;
  double sq = 1.0;
  unsigned long lq;

  for (i = 0, *k = 0; prod < main->nq && i < SQ_R->size; i++)
    {
      if (!polyselect_main_data_check_parameters(main, header->m0, sq * (double) SQ_R->q[i]))
	break;
      prod *= header->d;	/* We multiply by d instead of SQ_R->nr[i] to limit
				   the number of primes and thus the Y1 value. */
      sq *= (double) SQ_R->q[i];
      *k += 1;
    }

  /* We force k <= 4 on a 32-bit machine, and k <= 8 on a 64-bit machine,
     to ensure q fits on an "unsigned long". */
  if (*k > (sizeof(unsigned long) * CHAR_BIT) / 8)
    *k = (sizeof(unsigned long) * CHAR_BIT) / 8;
  if (*k < 1)
    *k = 1;

  /* If all factors in sq have d roots, then a single special-q is enough.
     Otherwise, we consider special-q's from combinations of k primes among lq,
     so that the total number of combinations is at least nq. */
  for (lq = *k; number_comb(SQ_R, *k, lq) < main->nq && lq < SQ_R->size; lq++);

  if (main->verbose > 0) {
    printf("# Info: nq=%lu leads to choosing k (k=%lu) primes among %lu (out of %u).\n",
            main->nq, *k, lq, SQ_R->size);
    printf("# Info: a combination of all %u possible primes leads to nq=%lu.\n",
            SQ_R->size, number_comb(SQ_R, *k, SQ_R->size));
  }

  return lq;
}

/* Given the current global state (on which we expect to have a lock!),
 * print the expected time to reach the given target E.
 *
 * This uses the Weibull moments functions. Note that the target_E
 * option is incompatible with maxtime (which also uses the same
 * things)
 */
static size_t snprintf_expected_time_target_E(
        char * buf, size_t size,
        polyselect_stats_ptr stats,
        double target_E)
{
    size_t np = 0;
    double beta, eta, prob;
    polyselect_data_series_srcptr exp_E = stats->exp_E;

    polyselect_data_series_estimate_weibull_moments2(&beta, &eta, exp_E);
    unsigned long collisions_good = stats->collisions_good;

    np += snprintf(buf + np, size - np, "# E:");
    ASSERT_ALWAYS(np <= size);

    np += polyselect_data_series_snprintf_summary(buf + np, size - np, exp_E);
    ASSERT_ALWAYS(np <= size);

    np += snprintf(buf + np, size - np, "\n");
    ASSERT_ALWAYS(np <= size);

    prob = 1.0 - exp(-pow(target_E / eta, beta));
    if (prob == 0)	/* for x small, exp(x) ~ 1+x */
        prob = pow(target_E / eta, beta);

    np += snprintf(buf + np, size - np,
            "# target_E=%.2f: collisions=%.2e, time=%.2e"
            " (beta %.2f,eta %.2f)\n", target_E, 1.0 / prob,
            seconds() / (prob * collisions_good), beta, eta);
    ASSERT_ALWAYS(np <= size);

    return np;
}

/* Given the current global state (on which we expect to have a lock!),
 * print what we can expect to reach after maxtime. We also take the ad
 * parameter, since it's handy to predict which ad we can expect to
 * reach.
 *
 * This uses the Weibull moments functions. Note that the maxtime
 * option is incompatible with target_E (which also uses the same
 * things)
 */
static size_t snprintf_expected_goal_maxtime(
        char * buf, size_t size,
        polyselect_stats_ptr stats,
        double maxtime,
        mpz_srcptr admin,
        mpz_srcptr ad)
{
    size_t np = 0;
    double beta, eta, prob;
    polyselect_data_series_srcptr exp_E = stats->exp_E;

    /* estimate the parameters of a Weibull distribution for E */
    polyselect_data_series_estimate_weibull_moments2(&beta, &eta, exp_E);
    unsigned long collisions_good = stats->collisions_good;

    np += snprintf(buf + np, size - np, "# E:");
    ASSERT_ALWAYS(np <= size);

    np += polyselect_data_series_snprintf_summary(buf + np, size - np, exp_E);
    ASSERT_ALWAYS(np <= size);

    np += snprintf(buf + np, size - np, "\n");
    ASSERT_ALWAYS(np <= size);

    /* time = seconds () / (prob * collisions_good)
       where  prob = 1 - exp (-(E/eta)^beta) */
    unsigned long n = collisions_good;	/* #polynomials found so far */
    double time_so_far = seconds();
    double time_per_poly = time_so_far / n;	/* average time per poly */
    double admin_d = mpz_get_d(admin);
    double ad_d = mpz_get_d(ad);
    double adrange = (ad_d - admin_d) * (maxtime / time_so_far);
    /* WTF ??? I'm commenting that line, but what does it mean ???
     * Is it a leftover from something?
     */
    // adrange = 2.00e+15 - 99900000000000.0;
    prob = time_so_far / (maxtime * n);
    double E = eta * pow(-log(1 - prob), 1.0 / beta);

    polyselect_data_series_ptr best_exp_E_Weibull = stats->best_exp_E_Weibull;
    polyselect_data_series_add(best_exp_E_Weibull, E);
    /* since the values of (eta,beta) fluctuate a lot, because
       they depend on the random samples in estimate_weibull_moments2,
       we take the average value for best_exp_E */
    E = polyselect_data_series_mean(best_exp_E_Weibull);

    np += snprintf(buf + np, size - np,
            "# %.2fs/poly, eta %.2f, beta %.3f, admax %.2e, best exp_E %.2f\n",
            time_per_poly, eta, beta, admin_d + adrange, E);
    return np;
}

void polyselect_main_data_commit_stats_unlocked(polyselect_main_data_ptr main, polyselect_stats_ptr stats, mpz_srcptr ad)
{
    polyselect_stats_accumulate(main->stats, stats);
    polyselect_stats_reset(stats);

    if (!ad) {
        /* if ad==NULL, we're committing stats for an asynchronous job,
         * meaning that we're not interested in the time projection stuff
         * below
         */
        return;
    }

    if (main->target_E != 0 && main->stats->exp_E->size) {
        char buf[1024];
        snprintf_expected_time_target_E(buf, sizeof(buf), main->stats, main->target_E);
        fputs(buf, stdout);
        fflush(stdout);
    }

    if (main->maxtime != DBL_MAX && main->stats->exp_E->size) {
        char buf[1024];
        snprintf_expected_goal_maxtime(buf, sizeof(buf), main->stats, main->maxtime, main->admin, ad);
        fputs(buf, stdout);
        fflush(stdout);
    }
}


void polyselect_main_data_commit_stats(polyselect_main_data_ptr main, polyselect_stats_ptr stats, mpz_srcptr ad)
{
    pthread_mutex_lock(&main->lock);
    polyselect_main_data_commit_stats_unlocked(main, stats, ad);
    pthread_mutex_unlock(&main->lock);
}

/* fetch N and d from the parameter list. It is of course necessary
 * before we can decide to setup many of the main_data parameters.
 */
void polyselect_main_data_parse_Nd(polyselect_main_data_ptr main, param_list_ptr pl)
{
    /* parse and check N in the first place */
    int have_n = param_list_parse_mpz(pl, "n", main->N);

    if (!have_n)
    {
        fprintf(stderr, "# Reading n from stdin\n");
        param_list_read_stream(pl, stdin, 0);
        have_n = param_list_parse_mpz(pl, "n", main->N);
    }

    if (!have_n)
    {
        fprintf(stderr, "No n defined ; sorry.\n");
        exit(EXIT_FAILURE);
    }

    if (mpz_cmp_ui(main->N, 0) <= 0)
        param_list_generic_failure(pl, "n");

    param_list_parse_uint(pl, "degree", &main->d);
    /* check degree */
    if (main->d <= 0)
        param_list_generic_failure(pl, "degree");
}

/*  parse incr, admin, admax */
void polyselect_main_data_parse_ad_range(polyselect_main_data_ptr main, param_list_ptr pl)
{
    param_list_parse_ulong(pl, "incr", &main->incr);
    if (main->incr <= 0)
    {
        fprintf(stderr, "Error, incr should be positive\n");
        exit(1);
    }

    /* if no -admin is given, mpz_init did set it to 0, which is exactly
       what we want */
    param_list_parse_mpz(pl, "admin", main->admin);
    /* admin should be nonnegative */
    if (mpz_cmp_ui(main->admin, 0) < 0)
    {
        fprintf(stderr, "Error, admin should be nonnegative\n");
        exit(1);
    }
    /* if admin = 0, start from incr */
    if (mpz_cmp_ui(main->admin, 0) == 0)
        mpz_set_ui(main->admin, main->incr);

    /* admin should be a non-zero multiple of 'incr', since when the global
       [admin, admax] range is cut by cado-nfs.py between different workunits,
       some bounds might no longer be multiple of 'incr'. */
    mpz_add_ui(main->admin, main->admin, (main->incr - mpz_fdiv_ui(main->admin, main->incr)) % main->incr);


    param_list_parse_mpz(pl, "admax", main->admax);
}

/* parse maxtime or target_E */
void polyselect_main_data_parse_maxtime_or_target(polyselect_main_data_ptr main, param_list_ptr pl)
{
    param_list_parse_double(pl, "maxtime", &main->maxtime);
    param_list_parse_double(pl, "target_E", &main->target_E);
    /* maxtime and target_E are incompatible */
    if (main->maxtime != DBL_MAX && main->target_E != 0.0)
    {
        fprintf(stderr, "Options -maxtime and -target_E are incompatible\n");
        exit(1);
    }
}

/* parse P, and check some conditions */
void polyselect_main_data_parse_P(polyselect_main_data_ptr main, param_list_ptr pl)
{
    unsigned long P = 0;

    param_list_parse_ulong(pl, "P", &P);
    if (P == 0)
        param_list_generic_failure(pl, "P");

    double Pd;
    Pd = (double) P;
    if (Pd > (double) UINT_MAX)
    {
        fprintf(stderr, "Error, too large value of P\n");
        exit(1);
    }

    if (4.0 * Pd * Pd >= (double) UINT64_MAX)
    {
        mpz_t tmp;
        mpz_init_set_uint64 (tmp, UINT64_MAX >> 2);
        mpz_sqrt (tmp, tmp);
        gmp_fprintf (stderr, "Error, too large value of P, maximum is %Zd\n",
                tmp);
        mpz_clear(tmp);
        exit (1);
    }

    if (P <= (unsigned long) SPECIAL_Q[LEN_SPECIAL_Q - 2])
    {
        fprintf(stderr, "Error, too small value of P, need P > %u\n",
                SPECIAL_Q[LEN_SPECIAL_Q - 2]);
        exit(1);
    }
    /* since for each prime p in [P,2P], we convert p^2 to int64_t, we need
       (2P)^2 < 2^63, thus P < 2^30.5 */
    /* XXX In fact there is a bug in collision_on_p which limit P to 2^30 */
    /* TODO: investigate, document. */
    if (P >= 1073741824UL)
    {
        fprintf(stderr, "Error, too large value of P\n");
        exit(1);
    }
    main->P = P;
}

unsigned long polyselect_main_data_number_of_ad_tasks(polyselect_main_data_srcptr main)
{

    if (mpz_cmp_ui(main->admax, 0) <= 0) {
        return ULONG_MAX;
    }

    unsigned long idx_max;
    /* ad = admin + idx * incr < admax
       thus idx_max = floor((admax-admin-1)/incr) + 1
       = floor((admax-admin+incr-1)/incr) */

    mpz_t t;
    mpz_init_set(t, main->admax);
    mpz_sub(t, t, main->admin);
    mpz_add_ui(t, t, main->incr - 1);
    mpz_fdiv_q_ui(t, t, main->incr);
    if (mpz_fits_ulong_p(t))
        idx_max = mpz_get_ui(t);
    else
        idx_max = ULONG_MAX;
    mpz_clear(t);

    return idx_max;
}


static void 
polyselect_main_data_auto_scale(polyselect_main_data_ptr main_data MAYBE_UNUSED)
{
#ifdef HAVE_HWLOC
    main_data->nthreads = hwloc_bitmap_weight(hwloc_get_root_obj(main_data->topology)->cpuset);
    const char * tmp = getenv("CADO_NFS_MAX_THREADS");
    if (tmp) {
        unsigned long const n = strtoul(tmp, NULL, 0);
        if (n < main_data->nthreads) {
            main_data->nthreads = n;
        }
    }
#endif
}

void polyselect_main_data_prepare_leagues(polyselect_main_data_ptr main_data)
{
    /* main_data->nnodes is already set by
     * polyselect_main_data_check_topology */

    /* prepare groups */
    main_data->leagues = malloc(main_data->nnodes * sizeof(polyselect_thread_league));
    printf("# %u nodes, 1 league on each\n", main_data->nnodes);

    for(unsigned int i = 0 ; i < main_data->nnodes ; i++) {
        /* This initializes the list of primes, one on each NUMA node */
        polyselect_thread_league_init(&(main_data->leagues[i]), main_data, i);
#ifdef HAVE_HWLOC
        if (main_data->bind && main_data->leagues[i].membind_set)
            printf("# node %u has %u cpus\n", i, hwloc_bitmap_weight(main_data->leagues[i].membind_set));
#endif
    }
}

void polyselect_main_data_dispose_leagues(polyselect_main_data_ptr main_data)
{
  for(unsigned int i = 0 ; i < main_data->nnodes ; i++) {
      polyselect_thread_league_clear(&main_data->leagues[i]);
  }
  free(main_data->leagues);
}

void polyselect_main_data_prepare_teams(polyselect_main_data_ptr main_data)
{
    unsigned int nteams = main_data->nthreads / main_data->finer_grain_threads;
    unsigned int w = nteams / main_data->nnodes;
    main_data->teams = malloc(nteams * sizeof(polyselect_thread_team));

    printf("# %u teams per league, %u team(s) in total\n", w, nteams);

    for(unsigned int i = 0 ; i < nteams ; i++) {
        polyselect_thread_team_ptr team = &main_data->teams[i];
        polyselect_thread_league_ptr league = &main_data->leagues[i / w];
        polyselect_thread_team_init(team, league, main_data, i);
    }
}

void polyselect_main_data_dispose_teams(polyselect_main_data_ptr main_data)
{
    unsigned int nteams = main_data->nthreads / main_data->finer_grain_threads;
    for(unsigned int i = 0 ; i < nteams ; i++) {
        polyselect_thread_team_ptr team = &main_data->teams[i];
        polyselect_thread_team_clear(team);
    }
    free(main_data->teams);
}

void polyselect_main_data_prepare_threads(polyselect_main_data_ptr main_data)
{
    main_data->threads = malloc(main_data->nthreads * sizeof(polyselect_thread));
    printf("# %u threads per team, %u threads in total\n", main_data->finer_grain_threads, main_data->nthreads);

    for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
        polyselect_thread_ptr thread = &main_data->threads[i];
        unsigned int ti = i / main_data->finer_grain_threads;
        polyselect_thread_team_ptr team = &(main_data->teams[ti]);
        polyselect_thread_init(thread, team, main_data, i);
    }

#ifdef HAVE_HWLOC
    /*  prepare the cpu sets for binding */
    if (main_data->bind) {
        if (hwloc_get_nbobjs_by_type(main_data->topology, HWLOC_OBJ_PU) % main_data->nthreads) {
            fprintf(stderr, "Warning, the number of threads is not compatible with this topology, we will not bind threads\n");
        } else {
            for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
                main_data->threads[i].cpubind_set = hwloc_bitmap_alloc();
            }
            int pu_depth = hwloc_get_type_depth(main_data->topology, HWLOC_OBJ_PU);
            unsigned int nk = hwloc_get_nbobjs_by_type(main_data->topology, HWLOC_OBJ_PU);
            int w = nk / main_data->nthreads;
            for(unsigned int i = 0 ; i < nk ; i++) {
                hwloc_bitmap_or(main_data->threads[i/w].cpubind_set, main_data->threads[i/w].cpubind_set, hwloc_get_obj_by_depth(main_data->topology,pu_depth, i)->cpuset);
            }
        }
    }
#endif
}

void polyselect_main_data_dispose_threads(polyselect_main_data_ptr main_data)
{
  for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
      polyselect_thread_clear(&main_data->threads[i]);
  }
#ifdef HAVE_HWLOC
  if (main_data->bind) {
    for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
       if(main_data->bind && main_data->threads[i].cpubind_set)
         hwloc_bitmap_free(main_data->threads[i].cpubind_set);
    }
  }
#endif
  free(main_data->threads);
}

void polyselect_main_data_go_parallel(polyselect_main_data_ptr main_data, void * (*thread_loop)(polyselect_thread_ptr))
{
    unsigned long idx_max = polyselect_main_data_number_of_ad_tasks(main_data);

    if (idx_max < main_data->nthreads)
    {
        fprintf(stderr,
                "# Warning: the current admin, admax, incr settings only make it possible to run %lu jobs in parallel, so that we won't be able to do %d-thread parallelism as requested\n",
                idx_max, main_data->nthreads);
    }

    for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
        pthread_create(&main_data->threads[i].tid, NULL, (void * (*)(void*)) thread_loop, &main_data->threads[i]);
    }
    for(unsigned int i = 0 ; i < main_data->nthreads ; i++) {
        pthread_join(main_data->threads[i].tid, NULL);
    }
}

void polyselect_main_data_check_topology(polyselect_main_data_ptr main_data)
{
    /* start with that as a default */
    main_data->nnodes = 1;
    if (main_data->nthreads) {
        /* if a number of threads was specified, this means that the job
         * placement is not our responsibility */
#ifdef HAVE_HWLOC
        main_data->bind = 0;
#endif
        return;
    }
#ifdef HAVE_HWLOC
    main_data->bind = 1;
    int mdepth = hwloc_get_type_depth(main_data->topology, HWLOC_OBJ_NUMANODE);
    /* hwloc 2 will answer HWLOC_TYPE_DEPTH_NUMANODE but hwloc 1.x
     * will return something in the topology, or possibly UNKNOWN */
    if (mdepth == HWLOC_TYPE_DEPTH_UNKNOWN)
        mdepth = 0;
    main_data->nnodes = hwloc_get_nbobjs_by_depth(main_data->topology, mdepth);
    /* mdepth == -1 can happen in some VM settings */
#endif

    polyselect_main_data_auto_scale(main_data);

    /*  sanity check nthreads / nnodes */
    if (main_data->nthreads % main_data->nnodes) {
        fprintf(stderr, "Error: the number of requested threads is incompatible with the number of nodes\n");
        /* what should we do here ? */
        /* perhaps restrict to one node and set loose binding ? */
        /* or just abort ? */
    }

    unsigned int nthreads_per_group = main_data->nthreads / main_data->nnodes;

    /*  sanity check nthreads_per_group / finer_grain_threads */
    if (nthreads_per_group % main_data->finer_grain_threads) {
        unsigned int f = main_data->finer_grain_threads;
        for( ; nthreads_per_group % f ; f--) ;
        main_data->finer_grain_threads = f;
        fprintf(stderr, "Warning, the number of finer-grain threads is incompatible with the current number of threads. Reducing finer-grain threads to %u\n", main_data->finer_grain_threads);
    }
}

