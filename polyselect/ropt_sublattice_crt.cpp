#include "cado.h"

#include <cmath>

#include "ropt_single_sublattice_priority_queue_impl.hpp"
#include "ropt_sublattice_crt.h"

struct ropt_crt_combiner_s {
    size_t nprimes;
    unsigned int * primes;
    unsigned int * e;
    float * logp;
    mpz_t modulus;
    mpz_t * multipliers;
};
typedef struct ropt_crt_combiner_s ropt_crt_combiner[1];
typedef struct ropt_crt_combiner_s * ropt_crt_combiner_ptr;
typedef const struct ropt_crt_combiner_s * ropt_crt_combiner_srctr;

void ropt_crt_combiner_init(ropt_crt_combiner_ptr R, size_t nprimes, const unsigned int * primes, unsigned int * e)
{
    R->nprimes = nprimes;
    R->primes = (unsigned int *) malloc(nprimes * sizeof(unsigned int));
    R->logp = (float *) malloc(nprimes * sizeof(float));
    R->multipliers = (mpz_t *) malloc(nprimes * sizeof(mpz_t));
    memcpy(R->primes, primes, nprimes * sizeof(unsigned int));
    R->e = (unsigned int *) malloc(nprimes * sizeof(unsigned int));
    memcpy(R->e, e, nprimes * sizeof(unsigned int));
    mpz_init(R->modulus);
    for(size_t i = 0 ; i < nprimes ; i++)
        mpz_init(R->multipliers[i]);

    mpz_set_ui(R->modulus, 1);

    mpz_t tmp;
    unsigned int * pe = (unsigned int *) malloc(nprimes * sizeof(unsigned int));
    mpz_init(tmp);

    for(size_t i = 0 ; i < nprimes ; i++) {
        R->logp[i] = log( (double) primes[i] );
        pe[i] = primes[i];
        if (e) {
            for(unsigned int j = 1 ; j < e[i] ; j++) {
                pe[i] *= primes[i];
            }
        }
        mpz_mul_ui(R->modulus, R->modulus, pe[i]);
    }
    /* It's fine with very small sets of primes, but this approach is
     * sub-par when nprimes becomes large */
    for (size_t i = 0; i < nprimes; i ++) {
        mpz_divexact_ui (R->multipliers[i], R->modulus, pe[i]);
        mpz_set_ui (tmp, pe[i]);
        mpz_invert (tmp, R->multipliers[i], tmp);
        mpz_mul (R->multipliers[i], R->multipliers[i], tmp);
    }

    mpz_clear(tmp);
    free(pe);
}

void ropt_crt_combiner_clear(ropt_crt_combiner_ptr R)
{
    for(size_t i = 0 ; i < R->nprimes ; i++)
        mpz_clear(R->multipliers[i]);
    mpz_clear(R->modulus);
    free(R->e);
    free(R->primes);
    free(R->logp);
    free(R->multipliers);
}

/**
 * Compute crt and add (u, v) to queue.
 *
 * infos[i]: lattice modulo tprimes[i], to power infos[i]->e
 */
static unsigned int
return_combined_sublattice_crt ( unsigned int nprimes,
                                 unsigned int const *primes,
                                 ropt_bound_srcptr bound,
                                 single_sublattice_info const * const * infos,
                                 sublattice_priority_queue_ptr pqueue )
{
  mpz_t sum, tmpu1, tmpu2;
  ropt_crt_combiner R;
  unsigned int * e = (unsigned int *) malloc (nprimes * sizeof(unsigned int));

  for(size_t i = 0 ; i < nprimes ; i++)
      e[i] = infos[i]->e;

  ropt_crt_combiner_init(R, nprimes, primes, e);

  /* old code was setting this, is it useful ? */
  // mpz_set(s1param->modulus, R->modulus);
  
  mpz_init (sum);
  mpz_init (tmpu1);
  mpz_init (tmpu2);

  /* compute u */
  mpz_set_ui (sum, 0);
  for (size_t i = 0; i < nprimes; i ++) {
      mpz_addmul_ui(sum, R->multipliers[i], infos[i]->u);
  }
  mpz_mod (tmpu1, sum, R->modulus);
  mpz_sub (tmpu2, R->modulus, tmpu1); // tmpu2 > 0
  
  unsigned int count = 0;

  /* check if a signed, centered representative of u falls within
   * global_u_boundl .. global_u_boundr
   */
  /* if u is good, compute v */
  if ( mpz_cmp_si (tmpu1, bound->global_u_boundr) <= 0 ||
       mpz_cmp_si (tmpu2, -bound->global_u_boundl) <= 0 ) {

      /* compute v */
      float val = 0.0;
      mpz_set_ui (sum, 0);
      for (size_t i = 0; i < nprimes; i ++) {
          mpz_addmul_ui(sum, R->multipliers[i], infos[i]->v);
          val += R->logp[i] * infos[i]->val;
      }
      mpz_mod (tmpu2, sum, R->modulus);

      /* insert this node */
      sublattice_priority_queue_push(pqueue, tmpu1, tmpu2, R->modulus, val );
      count = 1;
  }

  ropt_crt_combiner_clear(R);
  free(e);
  mpz_clear (sum);
  mpz_clear (tmpu1);
  mpz_clear (tmpu2);

  return count;
}

unsigned int ropt_sublattice_combine_all_crt_inner(unsigned int nprimes, unsigned int * primes, single_sublattice_priority_queue * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue,
        single_sublattice_info const ** infos, unsigned int i)
{
    if (i == nprimes) {
        return return_combined_sublattice_crt (nprimes, primes, bound, infos, pqueue);
    }
    /* iterate through all elements of tops[i], put them one by one in
     * infos[i], and recurse on that.
     */
    auto qi = single_sublattice_priority_queue_impl::cast(tops[i]);
    unsigned int count = 0;
    for(auto const & q : *qi) {
        infos[i] = &q;
        count += ropt_sublattice_combine_all_crt_inner(nprimes, primes, tops, bound, pqueue, infos, i+1);
    }
    return count;
}

unsigned int ropt_sublattice_combine_all_crt(unsigned int nprimes, unsigned int * primes, single_sublattice_priority_queue * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue)
{
    std::vector<single_sublattice_info const *> infos(nprimes, nullptr);
    return ropt_sublattice_combine_all_crt_inner(nprimes, primes, tops, bound, pqueue, infos.data(), 0);
}
