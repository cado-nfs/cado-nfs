#include "polyselect2l_arith.h"
#include "polyselect2l_str.h"

/* return 1/a mod p, assume 0 <= a < p */
unsigned long
invert (unsigned long a, unsigned long p)
{
  modulusul_t q;
  residueul_t b;

  modul_initmod_ul (q, p);
  modul_init (b, q);
  assert (a < p);
  modul_set_ul_reduced (b, a, q);
  modul_inv (b, b, q);
  a = modul_get_ul (b, q);
  modul_clear (b, q);
  modul_clearmod (q);
  return a;
}


/* lift the n roots r[0..n-1] of N = x^d (mod p) to roots of
   N = (m0 + r)^d (mod p^2) */
void
roots_lift (uint64_t *r, mpz_t N, unsigned long d, mpz_t m0,
            unsigned long p, unsigned long n)
{
  uint64_t pp;
  unsigned long j, inv;
  mpz_t tmp, lambda;
  mpz_init (tmp);
  mpz_init (lambda);
  pp = (uint64_t) p;
  pp *= (uint64_t) p;

  if (sizeof (unsigned long) == 8) {
    for (j = 0; j < n; j++) {
	/* we have for r=r[j]: r^d = N (mod p), lift mod p^2:
	   (r+lambda*p)^d = N (mod p^2) implies
	   r^d + d*lambda*p*r^(d-1) = N (mod p^2)
	   lambda = (N - r^d)/(p*d*r^(d-1)) mod p */
	mpz_ui_pow_ui (tmp, r[j], d - 1);
	mpz_mul_ui (lambda, tmp, r[j]);    /* lambda = r^d */
	mpz_sub (lambda, N, lambda);
	mpz_divexact_ui (lambda, lambda, p);
	mpz_mul_ui (tmp, tmp, d);         /* tmp = d*r^(d-1) */
	inv = invert (mpz_fdiv_ui (tmp, p), p);
	mpz_mul_ui (lambda, lambda, inv * p); /* inv * p fits in 64 bits if
						 p < 2^32 */
	mpz_add_ui (lambda, lambda, r[j]); /* now lambda^d = N (mod p^2) */

	/* subtract m0 to get roots of (m0+r)^d = N (mod p^2) */
	mpz_sub (lambda, lambda, m0);
	r[j] = mpz_fdiv_ui (lambda, pp);
      }
  }
  else {
#if 0   
    printf ("p: %lu, ppl %"PRId64": ", p, pp);
#endif
    uint64_t tmp1;
    mpz_t ppz, *rz, tmpz;
    rz = (mpz_t*) malloc (n * sizeof (mpz_t));
    mpz_init (ppz);
    mpz_init (tmpz);
    for (j = 0; j < n; j++) {
      mpz_init_set_ui (rz[j], 0UL);
      mpz_set_uint64 (rz[j], r[j]);
#if 0   
      printf (" %"PRIu64"", r[j]);
#endif
    }

    for (j = 0; j < n; j++) {
	mpz_pow_ui (tmp, rz[j], d - 1);
	mpz_mul (lambda, tmp, rz[j]);    /* lambda = r^d */
	mpz_sub (lambda, N, lambda);
	mpz_divexact_ui (lambda, lambda, p);
	mpz_mul_ui (tmp, tmp, d);         /* tmp = d*r^(d-1) */
	inv = invert (mpz_fdiv_ui (tmp, p), p);
	tmp1 = (uint64_t) inv;
	tmp1 *= (uint64_t) p;
	mpz_set_uint64 (tmpz, tmp1);
	mpz_mul (lambda, lambda, tmpz); 
	mpz_add (lambda, lambda, rz[j]); /* now lambda^d = N (mod p^2) */
	/* subtract m0 to get roots of (m0+r)^d = N (mod p^2) */
	mpz_sub (lambda, lambda, m0);
	mpz_set_uint64 (tmpz, pp);
	mpz_fdiv_r (rz[j], lambda, tmpz);
	r[j] = mpz_get_uint64 (rz[j]);
      }

    for (j = 0; j < n; j++)
      mpz_clear (rz[j]);
    free (rz);
    mpz_clear (ppz);
    mpz_clear (tmpz);
  }

  mpz_clear (tmp);
  mpz_clear (lambda);
}


/* first combination: 0, 1, 2, 3, \cdots */
void
first_comb ( unsigned long k,
             unsigned long *r )
{
  unsigned long i;
  for (i = 0; i < k; ++i)
    r[i] = i;
}


/* next combination */
unsigned long
next_comb ( unsigned long n,
            unsigned long k,
            unsigned long *r )
{
  unsigned long j;

  /* if the last combination */
  if (r[0] == (n - k)) {
    return k;
  }

  /* r[k-1] doesnot equal to the n-1, just increase it */
  j = k - 1;
  if (r[j] < (n-1)) {
    r[j] ++;
    return j;
  }

  /* find which one we should increase */
  while ( (r[j] - r[j-1]) == 1)
    j --;

  unsigned long ret = j - 1;
  unsigned long z = ++r[j-1];

  while (j < k) {
    r[j] = ++z;
    j ++;
  }
  return ret;
}


/* debug */
void
print_comb ( unsigned long k,
             unsigned long *r )
{
  unsigned long i;
  for (i = 0; i < k; i ++)
    fprintf (stderr, "%lu ", r[i]);
  fprintf (stderr, "\n");
}


/* return number of n choose k */
unsigned long
binom ( unsigned long n,
        unsigned long k )
{
  if (k > n)
    return 0;
  if (k == 0 || k == n)
    return 1;
  if (2*k > n)
    k = n - k;

  unsigned long tot = n - k + 1, f = tot, i;
  for (i = 2; i <= k; i++) {
    f ++;
    tot *= f;
    tot /= i;
  }
  return tot;
}

/* sort r[0..n-1] in increasing order */
static void
roots_sort (uint64_t *r, unsigned long n)
{
  unsigned long i, j;
  uint64_t v;

  for (i = 1; i < n; i++)
    {
      v = r[i];
      for (j = i; j > 0 && v < r[j - 1]; j--)
        r[j] = r[j - 1];
      r[j] = v;
    }
}

/* prepare special-q's roots */
void
comp_sq_roots ( header_t header,
                qroots_t SQ_R )
{
  unsigned long i, q, nrq;
  uint64_t *rq;
  modulusul_t qq;
  mpz_t *f;

  f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
  if (f == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in comp_sq_r\n");
    exit (1);
  }
  for (i = 0; i <= header->d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[header->d], 1);

  rq = (uint64_t*) malloc (header->d * sizeof (uint64_t));
  if (rq == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in comp_sq_q\n");
    exit (1);
  }

  /* prepare the special-q's */
  for (i = 1; (q = SPECIAL_Q[i]) != 0 ; i++)
  {
    if ((header->d * header->ad) % q == 0)
      continue;

    if ( mpz_fdiv_ui (header->Ntilde, q) == 0 )
      continue;

    modul_initmod_ul (qq, q * q);
    mpz_mod_ui (f[0], header->Ntilde, q);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrq = poly_roots_uint64 (rq, f, header->d, q);
    /* sort roots to get reproducible results with the same seed
       (since poly_roots_uint64 also uses the random generator) */
    roots_sort (rq, nrq);
    roots_lift (rq, header->Ntilde, header->d, header->m0, q, nrq);

#if 0
    unsigned int j = 0;
    mpz_t r1, r2, *rqz;
    mpz_init (r1);
    mpz_init (r2);
    rqz = (mpz_t*) malloc (header->d * sizeof (mpz_t));
    mpz_set (r1, header->Ntilde);
    mpz_mod_ui (r1, r1, q*q);

    gmp_fprintf (stderr, "Ntilde: %Zd, Ntilde (mod %u): %Zd\n", 
                header->Ntilde, q*q, r1);

    for (j = 0; j < nrq; j ++) {
      mpz_set (r2, header->m0);

      mpz_add_ui (r2, r2, rq[j]);
      mpz_pow_ui (r2, r2, header->d);
      mpz_mod_ui (r2, r2, q*q);

      if (mpz_cmp (r1, r2) != 0) {
        fprintf (stderr, "Root computation wrong in comp_sq_roots().\n");
        fprintf (stderr, "q: %lu, rq: %lu\n", q, rq[j]);
        exit(1);
      }
    }

    mpz_clear (r1);
    mpz_clear (r2);
#endif

    qroots_add (SQ_R, q, nrq, rq);
    modul_clearmod (qq);
  }

  /* Reorder R entries by nr. It is safe to comment it. */
  qroots_rearrange (SQ_R);

  free(rq);
  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  qroots_realloc (SQ_R, SQ_R->size); /* free unused space */
}


/* given individual q's, return crted rq */
uint64_t
return_q_rq ( qroots_t SQ_R,
              unsigned long *idx_q,
              unsigned long k,
              mpz_t qqz,
              mpz_t rqqz,
              unsigned long lq )
{
  unsigned long i, j, idv_q[k], idv_rq[k];
  uint64_t q = 1;

  /* q and roots */
  for (i = 0; i < k; i ++) {
    idv_q[i] = SQ_R->q[idx_q[i]];
    q = q * idv_q[i];
    j = rand() % SQ_R->nr[idx_q[i]];
    idv_rq[i] = SQ_R->roots[idx_q[i]][j];
  }

#if 0
  for (i = 0; i < k; i ++) {
    fprintf (stderr, "(%lu:%lu) ", idv_q[i], idv_rq[i]);
  }
  //gmp_fprintf (stderr, "%Zd\n", rqqz);
#endif

  /* crt roots */
  crt_sq (qqz, rqqz, idv_q, idv_rq, lq);

  return q;
}


/* given individual q's, return \product q and q^2 only, no rq */
uint64_t
return_q_norq ( qroots_t SQ_R,
                unsigned long *idx_q,
                unsigned long k,
                mpz_t qqz )
{
  unsigned long i;
  uint64_t q = 1;

  for (i = 0; i < k; i ++)
    q = q * SQ_R->q[idx_q[i]];
  mpz_set_uint64 (qqz, q);
  mpz_mul (qqz, qqz, qqz);
  return q;
}

