#include "cado.h"
#include <inttypes.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gmp.h>
#include "polyselect_main_data.h"
#include "polyselect_main_queue.h"      // some useful defaults
#include "polyselect_arith.h"
#include "polyselect_hash.h"
#include "size_optimization.h"          // SOPT_DEFAULT_EFFORT
#include "timing.h"     // milliseconds
#include "getprime.h"   // getprime
#include "misc.h"       // nprimes_interval
#include "auxiliary.h"       // ALG_SIDE RAT_SIDE       (TODO: get rid of that)
#include "params.h"

//#define DEBUG_POLYSELECT
//#define DEBUG_POLYSELECT2



void polyselect_main_data_init_defaults(polyselect_main_data_ptr main)
{
    main->incr = DEFAULT_INCR;
    mpz_init(main->admin);
    mpz_init(main->admax);
    mpz_init(main->N);
    cado_poly_init(main->best_poly);
    cado_poly_init(main->curr_poly);
    main->Primes = NULL;
    main->lenPrimes = 0;
    pthread_mutex_init(&main->stats_lock, NULL);
    polyselect_stats_init(main->stats);
    main->verbose = 0;
    main->sopt_effort = SOPT_DEFAULT_EFFORT;
    main->maxtime = DBL_MAX;
    main->target_E = 0;
}


void polyselect_main_data_clear(polyselect_main_data_ptr main)
{
    mpz_clear(main->admin);
    mpz_clear(main->admax);
    mpz_clear(main->N);
    cado_poly_clear(main->best_poly);
    cado_poly_clear(main->curr_poly);
    free(main->Primes);
    pthread_mutex_destroy(&main->stats_lock);
    polyselect_stats_clear(main->stats);
}

/* This returns the bound on i in the algorithm. We choose to set it to
 * the square of the max prime (this is also a suggestion in the slides)
 */
int64_t polyselect_main_data_get_M(polyselect_main_data_srcptr main)
{
    int64_t p = main->Primes[main->lenPrimes - 1];
    return p * p;
}

size_t polyselect_main_data_expected_number_of_pairs(polyselect_main_data_srcptr main)
{
  uint32_t P = main->Primes[0];
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
    double p = main->Primes[main->lenPrimes-1];
    return pow(p, 4) * q < mpz_get_d(m0);
}

/* find suitable values of lq and k, where the special-q part in the degree-1
   coefficient of the linear polynomial is made from k small primes among lq */
unsigned long
find_suitable_lq(polyselect_poly_header_srcptr header,
		 polyselect_qroots_srcptr SQ_R,
                 unsigned long *k,
                 polyselect_main_data_ptr main)
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

  if (main->verbose) {
    printf("# Info: nq=%lu leads to choosing k (k=%lu) primes among %lu (out of %u).\n",
            main->nq, *k, lq, SQ_R->size);
    printf("# Info: a combination of all %u possible primes leads to nq=%lu.\n",
            SQ_R->size, number_comb(SQ_R, *k, SQ_R->size));
  }

  return lq;
}

void polyselect_main_data_commit_stats(polyselect_main_data_ptr main, polyselect_stats_ptr stats)
{
    pthread_mutex_lock(&main->stats_lock);
    polyselect_stats_accumulate(main->stats, stats);
    pthread_mutex_unlock(&main->stats_lock);
}

/* init prime array */
/* initialize primes in [P,2*P] */
static unsigned long
initPrimes ( unsigned long P,
             uint32_t **primes )
{
  unsigned long p, nprimes = 0;
  unsigned long Pmax = 2*P;
#ifdef LESS_P // if impatient for root finding
  Pmax = P + P/2;
#endif
  unsigned long maxprimes = nprimes_interval(P, Pmax);

  *primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  prime_info pi;
  prime_info_init (pi);

  /* It's now fairly trivial to parallelize this prime search if we want
   * to, but I think that it's a trivial computation anyway, and most
   * probably not worth the work.
   */
  prime_info_seek(pi, P);

  for (p = P, nprimes = 0; (p = getprime_mt (pi)) <= Pmax; nprimes++) {
    if (nprimes + 1 >= maxprimes) {
      maxprimes += maxprimes / 10;
      *primes = (uint32_t*) realloc (*primes, maxprimes * sizeof (uint32_t));
      if ( (*primes) == NULL) {
        fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
        exit (1);
      }
    }
    (*primes)[nprimes] = p;
  }

  prime_info_clear (pi);

  *primes = (uint32_t*) realloc (*primes, (nprimes) * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  return nprimes;
}


/* clear prime array */
void
polyselect_main_data_print_primes (polyselect_main_data_srcptr main)
{
    uint32_t *primes = main->Primes;
    unsigned long size = main->lenPrimes;
    unsigned long i;
    for (i = 0; i < size; i++) {
        fprintf (stderr, "(%lu, %" PRIu32 ") ", i, primes[i]);
        if ((i+1) % 5 == 0)
            fprintf (stderr, "\n");
    }
    fprintf (stderr, "\n");
}

void polyselect_main_data_prepare_primes(polyselect_main_data_ptr main)
{
    unsigned long st = milliseconds();

    main->lenPrimes = initPrimes(main->P, &main->Primes);

    printf("# Info: initializing %lu P primes took %lums, nq=%lu\n",
            main->lenPrimes, milliseconds() - st, main->nq);
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

    /* set cpoly */
    mpz_set(main->best_poly->n, main->N);
    mpz_set(main->curr_poly->n, main->N);

    /* do we have cado_poly_* interface calls for these? */
    main->best_poly->pols[ALG_SIDE]->deg = main->d;
    main->best_poly->pols[RAT_SIDE]->deg = 1;
    main->curr_poly->pols[ALG_SIDE]->deg = main->d;
    main->curr_poly->pols[RAT_SIDE]->deg = 1;
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

unsigned long polyselect_main_data_number_of_ad_tasks(polyselect_main_data_ptr main)
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

