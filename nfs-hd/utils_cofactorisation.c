#include "cado.h"
#include "utils_cofactorisation.h"
#include "utils_norm.h"

#include "ecm/facul.h"
#include "ecm/facul_doit.h"

void factor_init(factor_ptr factor, unsigned int alloc)
{
  ASSERT(alloc > 0);

  factor->alloc = alloc;
  factor->number = 0;
  factor->factorization = (mpz_t * ) malloc(sizeof(mpz_t) * alloc);
  for (unsigned int i = 0; i < alloc; i++) {
    mpz_init(factor->factorization[i]);
  }
}

void factor_clear(factor_ptr factor)
{
  for (unsigned int i = 0; i < factor->alloc; i++) {
    mpz_clear(factor->factorization[i]);
  }
  free(factor->factorization);
  factor->alloc = 0;
  factor->number = 0;
}

void factor_append(factor_ptr factor, mpz_srcptr z)
{
  if (factor->alloc == factor->number) {
    unsigned int newalloc = factor->alloc + 10;
    factor->factorization = (mpz_t * ) realloc(factor->factorization,
        sizeof(mpz_t) * newalloc);
    for (unsigned int i = factor->alloc; i < newalloc; ++i)
      mpz_init(factor->factorization[i]);
    factor->alloc = newalloc;
  }

  mpz_set(factor->factorization[factor->number], z);
  factor->number++;
}

unsigned int factor_remove(factor_ptr factor, mpz_srcptr z)
{
  unsigned int found = 0;
  unsigned int i = 0;
  for ( ; i < factor->number; i++) {
    if (mpz_cmp(z, factor->factorization[i]) == 0) {
      found = 1;
      break;
    }
  }
  if (found == 0) {
    return 0;
  }

  for ( ; i < factor->number - 1; i++) {
    mpz_set(factor->factorization[i], factor->factorization[i + 1]);
  }
  factor->number--;
  return 1;
}

void factor_fprintf(FILE * file, factor_srcptr factor)
{
  fprintf(file, "[");
  if (factor->number != 0) {
    for (unsigned int i = 0; i < factor->number - 1; i++) {
      gmp_fprintf(file, "%Zd, ", factor->factorization[i]);
    }
    gmp_fprintf(file, "%Zd]\n", factor->factorization[factor->number - 1]);
  } else {
    fprintf(file, "]\n");
  }
}

unsigned int factor_is_smooth(factor_srcptr factor, mpz_t B,
    unsigned int sort)
{
  if (sort) {
    ASSERT(sort == 1);

    if (mpz_cmp(factor->factorization[factor->number - 1], B) > 0) {
      return 0;
    }
  } else {
    ASSERT(sort == 0);

    for (unsigned int i = 0; i < factor->number; i++) {
      if (mpz_cmp(factor->factorization[i], B) > 0) {
        return 0;
      }
    }
  }

  return 1;
}

unsigned int factor_assert(factor_srcptr factor, mpz_srcptr z)
{
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(tmp, 1);
  for (unsigned int i = 0; i < factor->number; i++) {
    mpz_mul(tmp, tmp, factor->factorization[i]);
  }
  unsigned int assert_facto = 0;
  if (!mpz_cmp(tmp, z)) {
    ASSERT(mpz_cmp(tmp, z) == 0);

    assert_facto = 1;
  }
  mpz_clear(tmp);

  return assert_facto;
}

int brute_force_factorize_ul(factor_ptr factor, mpz_ptr z,
    mpz_srcptr z_root, unsigned long bound)
{
  prime_info pi;
  prime_info_init (pi);

  mpz_set(z, z_root);
  mpz_t prime_Z;
  mpz_init(prime_Z);
  unsigned long prime = 2;

  for (prime = 2; prime <= bound; prime = getprime_mt(pi)) {
    if (mpz_cmp_ui(z, 1) != 0) {
      mpz_set_ui(prime_Z, prime);
      mpz_t q;
      mpz_init(q);
      mpz_t r;
      mpz_init(r);
      mpz_fdiv_qr(q, r, z, prime_Z);
      while (mpz_cmp_ui(r, 0) == 0) {
        mpz_set(z, q);
        factor_append(factor, prime_Z);
        mpz_fdiv_qr(q, r, z, prime_Z);
      }
      mpz_clear(q);
      mpz_clear(r);
    }
  }

  mpz_clear(prime_Z);
  prime_info_clear (pi);

  if (mpz_cmp_ui(z, 1) == 0) {
    return 1;
  }
  return 0;
}

/*
 * To sort the factors by ascending values.
 */
static int compare_factor(const void * p0, const void * p1)
{
  return(mpz_cmp((mpz_srcptr) p0, (mpz_srcptr) p1));
}

/*
 * To sort by acending value the element of the factor array.
 *
 * factor: array of factors.
 */
static void sort_factor(factor_ptr factor)
{
  qsort(factor->factorization, factor->number,
      sizeof(factor->factorization[0]), compare_factor);
}

#ifndef NOT_PRINT_RELATION
/*
 * Print a relation.
 *
 * factor: factorisation the smooth norms.
 * I:
 * L:
 * a: the polynomial in the original lattice.
 * t: dimension of the lattice.
 * V: number of number fields.
 * size:
 * assert_facto: 1 if the factorisation is correct, 0 otherwise.
 * c: coordinate of a in the q-lattice.
 */
static void printf_relation(factor_t * factor,
    mpz_poly_srcptr a, unsigned int t, unsigned int V,
    FILE * outstd, MAYBE_UNUSED unsigned int * assert_facto,
    MAYBE_UNUSED mpz_vector_srcptr c, int * smooth)
{
  //Print a.
#ifdef PRINT_POLY_RELATION
  fprintf(outstd, "# ");
  mpz_poly_fprintf(outstd, a);
#endif // PRINT_POLY_RELATION

#ifdef PRINT_VECTOR_RELATION
  fprintf(outstd, "# ");
  mpz_vector_fprintf(outstd, c);
#endif // PRINT_VECTOR_RELATION

  //Print if the factorisation is good.
#ifdef ASSERT_FACTO
  fprintf(outstd, "# ");
  for (unsigned int i = 0; i < V - 1; i++) {
    ASSERT(assert_facto[i] < 3);
    fprintf(outstd, "%u:", assert_facto[i]);
  }
  ASSERT(assert_facto[V - 1] < 3);
  fprintf(outstd, "%u\n", assert_facto[V - 1]);
#endif // ASSERT_FACTO

  //Print the coefficient of a.
  for (int i = 0; i < a->deg; i++) {
    gmp_fprintf(outstd, "%Zd,", a->coeff[i]);
  }
  if ((int)t - 1 == a->deg) {
    gmp_fprintf(outstd, "%Zd:", a->coeff[a->deg]);
  } else {
    gmp_fprintf(outstd, "%Zd,", a->coeff[a->deg]);
    for (int i = a->deg + 1; i < (int)t - 1; i++) {
      fprintf(outstd, "0,");
    }
    fprintf(outstd, "0:");
  }

  //Print the factorisation in the different number field.
  for (unsigned int i = 0; i < V - 1; i++) {
    if (factor[i]->number != 0 && smooth[i]) {
      for (unsigned int j = 0; j < factor[i]->number - 1; j++) {
#ifdef PRINT_DECIMAL
        gmp_fprintf(outstd, "%Zd,", factor[i]->factorization[j]);
#else // PRINT_DECIMAL
        mpz_out_str(outstd, 16, factor[i]->factorization[j]);
        fprintf(outstd, ",");
#endif // PRINT_DECIMAL
      }
#ifdef PRINT_DECIMAL
      gmp_fprintf(outstd, "%Zd:", factor[i]->factorization[
          factor[i]->number - 1]);
#else // PRINT_DECIMAL
      mpz_out_str(outstd, 16,
          factor[i]->factorization[factor[i]->number - 1]);
      fprintf(outstd, ":");
#endif // PRINT_DECIMAL
    } else {
      fprintf(outstd, ":");
    }
  }

  if (factor[V - 1]->number != 0 && smooth[V - 1]) {
    for (unsigned int j = 0; j < factor[V - 1]->number - 1; j++) {
#ifdef PRINT_DECIMAL
      gmp_fprintf(outstd, "%Zd,", factor[V - 1]->factorization[j]);
#else // PRINT_DECIMAL
      mpz_out_str(outstd, 16, factor[V - 1]->factorization[j]);
      fprintf(outstd, ",");
#endif // PRINT_DECIMAL
    }
#ifdef PRINT_DECIMAL
    gmp_fprintf(outstd, "%Zd:", factor[V - 1]->factorization[
        factor[V - 1]->number - 1]);
#else // PRINT_DECIMAL
    mpz_out_str(outstd, 16,
        factor[V - 1]->factorization[factor[V - 1]->number - 1]);
#endif // PRINT_DECIMAL
  }
  fprintf(outstd, "\n");
}
#endif // NOT_PRINT_RELATION


typedef struct {
  unsigned long lpb;            // large prime bound = 2^lpb
  unsigned long fbb;            // fbb (the real bound, not its log)
  double BB;                    // square of fbb
  double BBB;                   // cube of fbb
  facul_method_t * methods;     // list of ECMs (null-terminated)
} facul_aux_data;

// FIXME: this does not belong to here!
static int
nb_curves95 (const unsigned int lpb)
{
    /* same, but with target probability 95% */
    /* do_table(10,64,ntries=500,target_prob=0.95)
     */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
/* lpb=10 */ 0, /* 0:1.000000 */
/* lpb=11 */ 0, /* 0:1.000000 */
/* lpb=12 */ 0, /* 0:1.000000 */
/* lpb=13 */ 0, /* 0:1.000000 */
/* lpb=14 */ 0, /* 0:1.000000 */
/* lpb=15 */ 0, /* 0:0.998000 */
/* lpb=16 */ 0, /* 0:0.986000 */
/* lpb=17 */ 0, /* 0:0.972000 */
/* lpb=18 */ 0, /* 0:0.964000 */
/* lpb=19 */ 1, /* 0:0.926000, 1:0.986000 */
/* lpb=20 */ 1, /* 1:0.986000 */
/* lpb=21 */ 2, /* 1:0.946000, 2:0.976000 */
/* lpb=22 */ 3, /* 2:0.940000, 3:0.966000 */
/* lpb=23 */ 3, /* 3:0.952000 */
/* lpb=24 */ 5, /* 4:0.926000, 5:0.962000 */
/* lpb=25 */ 7, /* 6:0.934000, 7:0.956000 */
/* lpb=26 */ 8, /* 7:0.934000, 8:0.956000 */
/* lpb=27 */ 10, /* 9:0.936000, 10:0.956000 */
/* lpb=28 */ 13, /* 12:0.938000, 13:0.954000 */
/* lpb=29 */ 16, /* 15:0.946000, 16:0.956000 */
/* lpb=30 */ 18, /* 17:0.936000, 18:0.950000 */
/* lpb=31 */ 22, /* 21:0.934000, 22:0.958000 */
/* lpb=32 */ 26, /* 25:0.940000, 26:0.950000 */
/* lpb=33 */ 29, /* 28:0.942000, 29:0.952000 */
/* lpb=34 */ 33, /* 32:0.948000, 33:0.956000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}

// FIXME: This has been duplicated from facul.cpp
// (September 2015).
// Maybe merge again, at some point.
// It returns -1 if the factor is not smooth, otherwise the number of
// factors.
// Remark: FACUL_NOT_SMOOTH is just -1.
static int
facul_aux (mpz_t *factors, const struct modset_t m,
    const facul_aux_data *data, int method_start)
{
  int found = 0;
  facul_method_t* methods = data->methods;
  if (methods == NULL)
    return found;

  int i = 0;
  for (i = method_start; methods[i].method != 0; i++)
  {
    struct modset_t fm, cfm;
    int res_fac = 0;

    switch (m.arith) {
      case CHOOSE_UL:
        res_fac = facul_doit_onefm_ul(factors, m.m_ul,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_15UL:
        res_fac = facul_doit_onefm_15ul(factors, m.m_15ul,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_2UL2:
        res_fac = facul_doit_onefm_2ul2 (factors, m.m_2ul2,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_MPZ:
        res_fac = facul_doit_onefm_mpz (factors, m.m_mpz,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      default: abort();
    }
    // check our result
    // res_fac contains the number of factors found
    if (res_fac == -1)
    {
      /*
         The cofactor m is not smooth. So, one stops the
         cofactorization.
         */
      found = FACUL_NOT_SMOOTH;
      break;
    }
    if (res_fac == 0)
    {
      /* Zero factor found. If it was the last method for this
         side, then one stops the cofactorization. Otherwise, one
         tries with an other method */
      continue;
    }

    found += res_fac;
    if (res_fac == 2)
      break;

    /*
       res_fac == 1  Only one factor has been found. Hence, our
       factorization is not finished.
       */
    if (fm.arith != CHOOSE_NONE)
    {
      int found2 = facul_aux(factors+res_fac, fm, data, i+1);
      if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
      {
        found = FACUL_NOT_SMOOTH;
        modset_clear(&cfm);
        modset_clear(&fm);
        break;
      }
      else
        found += found2;
      modset_clear(&fm);
    }
    if (cfm.arith != CHOOSE_NONE)
    {
      int found2 = facul_aux(factors+res_fac, cfm, data, i+1);
      if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
      {
        found = FACUL_NOT_SMOOTH;
      }
      else
        found += found2;
      modset_clear(&cfm);
      break;
    }
    break;
  }
  return found;
}


// Returns 1 if norm is smooth, 0 otherwise.
// Factorization is given in factors.
static int call_facul(factor_ptr factors, mpz_srcptr norm_r,
    facul_aux_data * data) {
  unsigned long B = data->fbb;

  mpz_t norm;
  mpz_init(norm);
  mpz_set(norm, norm_r);

  // Trial divide all the small factors.
  int success = brute_force_factorize_ul(factors, norm, norm, B);
  if (success) {
    mpz_clear(norm);
    return 1;
  }

  // TODO: move this test after the next block
  // and do the primality test with fast routines.
  if (mpz_probab_prime_p(norm, 1)) {
    if (mpz_sizeinbase(norm, 2) <= data->lpb) {
      factor_append(factors, norm);
      mpz_clear(norm);
      return 1;
    } else {
      mpz_clear(norm);
      return 0;
    }
  }

  // Prepare stuff for calling facul() machinery.
  // Choose appropriate arithmetic according to size of norm
  // Result is stored in n, which is of type struct modset_t (see facul.h)
  struct modset_t n;
  n.arith = CHOOSE_NONE;
  size_t bits = mpz_sizeinbase(norm, 2);
  if (bits <= MODREDCUL_MAXBITS) {
    ASSERT(mpz_fits_ulong_p(norm));
    modredcul_initmod_ul(n.m_ul, mpz_get_ui(norm));
    n.arith = CHOOSE_UL;
  }
  else if (bits <= MODREDC15UL_MAXBITS)
  {
    unsigned long t[2];
    modintredc15ul_t m;
    size_t written;
    mpz_export(t, &written, -1, sizeof(unsigned long), 0, 0, norm);
    ASSERT_ALWAYS(written <= 2);
    modredc15ul_intset_uls(m, t, written);
    modredc15ul_initmod_int(n.m_15ul, m);
    n.arith = CHOOSE_15UL;
  }
  else if (bits <= MODREDC2UL2_MAXBITS)
  {
    unsigned long t[2];
    modintredc2ul2_t m;
    size_t written;
    mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, norm);
    ASSERT_ALWAYS(written <= 2);
    modredc2ul2_intset_uls (m, t, written);
    modredc2ul2_initmod_int (n.m_2ul2, m);
    n.arith = CHOOSE_2UL2;
  }
  else
  {
    modmpz_initmod_int(n.m_mpz, norm);
    n.arith = CHOOSE_MPZ;
  }
  ASSERT_ALWAYS(n.arith != CHOOSE_NONE);

  // Call the facul machinery.
  // TODO: think about this hard-coded 16...
  mpz_t * fac = (mpz_t *) malloc(sizeof(mpz_t)*16);
  for (int i = 0; i < 16; ++i)
    mpz_init(fac[i]);
  int found = facul_aux(fac, n, data, 0);
  if (found > 0) {
    for (int i = 0; i < found; ++i)
      factor_append(factors, fac[i]);
  }

  for (int i = 0; i < 16; ++i)
    mpz_clear(fac[i]);
  free(fac);
  mpz_clear(norm);
  return found > 0;
}

static unsigned int is_irreducible(mpz_poly_srcptr a)
{
  if (a->deg == 1) {
    return 1;
  }
  else if (a->deg == 2) {
    mpz_t delta;
    mpz_init(delta);
    mpz_submul_ui(delta, a->coeff[0], 4);
    mpz_mul(delta, delta, a->coeff[2]);
    mpz_addmul(delta, a->coeff[1], a->coeff[1]);

    if (mpz_perfect_square_p(delta) == 0) {
      mpz_clear(delta);
      return 1;
    }

    mpz_clear(delta);
    return 0;
  } else {
    //TODO: wrong.
    return 1;
  }
}

#ifndef NOT_PRINT_RELATION
static void automorphism_6_1_1(mpz_poly_ptr b, mpz_poly_srcptr a)
{
  mpz_t * c = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (int i = 0; i < 2; i++) {
    mpz_init(c[i]);
  }

  mpz_add(c[1], a->coeff[0], a->coeff[1]);

  mpz_addmul_ui(c[0], a->coeff[0], 2);
  mpz_sub(c[0], c[0], a->coeff[1]);


  mpz_poly_setcoeffs(b, c, 1);
  for (int i = 0; i < 2; i++) {
    mpz_clear(c[i]);
  }
  free(c);
}

static void automorphism_6_1_2(mpz_poly_ptr b, mpz_poly_srcptr a)
{
  mpz_t * c = (mpz_t *) malloc(sizeof(mpz_t) * 3);
  for (int i = 0; i < 3; i++) {
    mpz_init(c[i]);
  }

  mpz_add(c[2], a->coeff[0], a->coeff[1]);
  mpz_add(c[2], c[2], a->coeff[2]);

  mpz_addmul_ui(c[1], a->coeff[0], 4);
  mpz_add(c[1], c[1], a->coeff[1]);
  mpz_submul_ui(c[1], a->coeff[2], 2);

  mpz_set(c[0], a->coeff[2]);
  mpz_submul_ui(c[0], a->coeff[1], 2);
  mpz_addmul_ui(c[0], a->coeff[0], 4);


  mpz_poly_setcoeffs(b, c, 2);
  for (int i = 0; i < 3; i++) {
    mpz_clear(c[i]);
  }
  free(c);
}

static void rewrite_poly_6_1(mpz_poly_ptr b, factor_t * fac, unsigned int V,
    int * smooth, factor_t * gal_norm_denom, const mpz_poly * f)
{
  for (unsigned int i = 0; i < V; i++) {
    if (fac[i]->number != 0 && smooth[i]) {
      for (unsigned int j = 0; j < (unsigned int)b->deg; j++) {
        for (unsigned int k = 0; k < gal_norm_denom[i]->number; k++) {
          factor_append(fac[i], gal_norm_denom[i]->factorization[k]);
        }
      }
      sort_factor(fac[i]);
    }
  }

  mpz_t c;
  mpz_init(c);
  mpz_poly_content(c, b);

  mpz_poly_divexact_mpz(b, b, c);

  if (mpz_cmp_ui(c, 1) > 0) {
    factor_t fac_tmp;
    factor_init(fac_tmp, 10);
    brute_force_factorize_ul(fac_tmp, c, c, 1000);

    ASSERT(mpz_cmp_ui(c, 1) == 0);

    for (unsigned int j = 0; j < fac_tmp->number; j++) {
      for (unsigned int i = 0; i < V; i++) {
        if (fac[i]->number != 0 && smooth[i]) {
          for (unsigned int k = 0; k < (unsigned int)f[i]->deg; k++) {
            MAYBE_UNUSED unsigned int found = factor_remove(fac[i],
                fac_tmp->factorization[j]);
            ASSERT(found == 1);
          }
        }
      }
    }

    factor_clear(fac_tmp);
  }

  if (mpz_cmp_ui(b->coeff[b->deg], 0) < 0) {
    mpz_poly_neg(b, b);
  }

  mpz_clear(c);
}

//TODO. assert facto is false: need to recompute.
//TODO: galois is probably not compatible with MNFS.
static void printf_relation_galois_6_1(factor_t * factor,
    mpz_poly_srcptr a, unsigned int t, unsigned int V,
    FILE * outstd, MAYBE_UNUSED unsigned int * assert_facto,
    MAYBE_UNUSED mpz_vector_srcptr c, int * smooth,
    const mpz_poly * f, factor_t * gal_norm_denom)
{
  if (a->deg > 2) {
    fprintf(outstd, "# Can not use Galois 6.1.\n");
  } else {
    fprintf(outstd, "# Use Galois 6.1.\n");
    mpz_poly b;
    mpz_poly_init(b, a->deg);
    mpz_poly_set(b, a);

    for (unsigned int i = 0; i < 5; i++) {
      if (a->deg == 2) {
        automorphism_6_1_2(b, b);
      } else if (a->deg == 1) {
        automorphism_6_1_1(b, b);
      }
      rewrite_poly_6_1(b, factor, V, smooth, gal_norm_denom, f);
#ifdef ASSERT_FACTO
      mpz_t res;
      mpz_init(res);
      for (unsigned int j = 0; j < V; j++) {
        if (factor[j]->number != 0 && smooth[j]) {
          norm_poly(res, f[j], b);
          assert_facto[j] = factor_assert(factor[j], res);
        }
      }
      mpz_clear(res);
#endif // ASSERT_FACTO
      printf_relation(factor, b, t, V, outstd, assert_facto, c, smooth);
    }
    mpz_poly_clear(b);
    fprintf(outstd, "# ----------\n");
  }
}
#endif // NOT_PRINT_RELATION

/*
 * A potential polynomial is found. Factorise it to verify if it gives a
 *  relation.
 *
 * a: the polynomial in the original lattice.
 * f: the polynomials that define the number fields.
 * lpb: the large prime bounds.
 * L:
 * size:
 * t: dimension of the lattice.
 * V: number of number fields.
 */
static void good_polynomial(mpz_poly_srcptr a, const mpz_poly * f, unsigned int
    * L, MAYBE_UNUSED unsigned int t, unsigned int V, int main, facul_aux_data
    *data, unsigned int * nb_rel_found, ideal_spq_srcptr special_q, unsigned int
    q_side, unsigned int size, MAYBE_UNUSED FILE * outstd, MAYBE_UNUSED unsigned
    int * number_factorisation, MAYBE_UNUSED mpz_vector_srcptr c,
    unsigned int gal, unsigned int * nb_rel_gal, unsigned int gal_version,
    factor_t * gal_norm_denom)
{
  mpz_t res;
  mpz_init(res);

  factor_t * factor = (factor_t * ) malloc(sizeof(factor_t) * V);
  int * smooth = (int *) malloc(sizeof(int) * V);

  MAYBE_UNUSED unsigned int * assert_facto;
#ifdef ASSERT_FACTO
  assert_facto = (unsigned int * ) malloc(sizeof(unsigned int) * V);
  mpz_t res_tmp;
  mpz_init(res_tmp);
#else // ASSERT_FACTO
  assert_facto = NULL;
#endif // ASSERT_FACTO

  unsigned int find = 0;
  if (main == -1) {
    unsigned int explore = 0;
    for (unsigned int i = 0; i < V; i++) {
      factor_init(factor[i], 10);
      if (L[i] && (explore < size - 1 || find != 0)) {
        explore++;
        norm_poly(res, f[i], a);

#ifdef ASSERT_FACTO
        mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

        if (i == q_side) {
          mpz_t spq;
          mpz_init(spq);
          mpz_set_uint64(spq, ideal_spq_get_q(special_q));
          ASSERT(ideal_spq_get_deg_g(special_q) >= 1);
          for (int deg = 0; deg < ideal_spq_get_deg_g(special_q); deg++) {
            mpz_divexact(res, res, spq);
            factor_append(factor[i], spq);
          }
          mpz_clear(spq);
        }

        smooth[i] = call_facul(factor[i], res, &data[i]);

#ifdef ASSERT_FACTO
        assert_facto[i] = factor_assert(factor[i], res_tmp);
        if (assert_facto[i] == 0) {
          number_factorisation[0] = number_factorisation[0] + 1;
#ifdef PRINT_ASSERT_FACTO
          gmp_fprintf(outstd, "# Incomplete: %Zd\n", res_tmp);
#endif // PRINT_ASSERT_FACTO
        } else {
          number_factorisation[1] = number_factorisation[1] + 1;
        }
#endif // ASSERT_FACTO

        if (smooth[i]) {
          find++;
          sort_factor(factor[i]);
        }
      } else {
#ifdef ASSERT_FACTO
        assert_facto[i] = 2;
#endif // ASSERT_FACTO
      }
    }
  } else {
    // TODO: this part of code is false for cofactorisation at least.
    ASSERT(main >= 0);

    factor_init(factor[main], 10);
    norm_poly(res, f[main], a);

#ifdef ASSERT_FACTO
    mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

    smooth[main] = call_facul(factor[main], res, &data[main]);

#ifdef ASSERT_FACTO
    assert_facto[main] = factor_assert(factor[main], res);
    if (assert_facto[main] == 0) {
      number_factorisation[0] = number_factorisation[0] + 1;
#ifdef PRINT_ASSERT_FACTO
      gmp_fprintf(outstd, "# Incomplete: %Zd\n", res_tmp);
#endif // PRINT_ASSERT_FACTO
    } else {
      number_factorisation[1] = number_factorisation[1] + 1;
    }
#endif // ASSERT_FACTO

    if (smooth[main]) {
      find = 1;

      for (unsigned int i = 0; i < V; i++) {
        if (i != (unsigned int) main) {
          factor_init(factor[i], 10);
          if (L[i]) {
            norm_poly(res, f[i], a);

#ifdef ASSERT_FACTO
            mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

            smooth[i] = call_facul(factor[i], res, &data[i]);

#ifdef ASSERT_FACTO
            assert_facto[i] = factor_assert(factor[i], res);
            if (assert_facto[i] == 0) {
              number_factorisation[0] = number_factorisation[0] + 1;
#ifdef PRINT_ASSERT_FACTO
              gmp_fprintf(outstd, "# Incomplete: %Zd\n", res_tmp);
#endif // PRINT_ASSERT_FACTO
            } else {
              number_factorisation[1] = number_factorisation[1] + 1;
            }
#endif // ASSERT_FACTO

            if (smooth[i]) {
              find++;
              sort_factor(factor[i]);
            }
          } else {
#ifdef ASSERT_FACTO
            assert_facto[i] = 2;
#endif // ASSERT_FACTO
          }
        } 
      }
    }
  }

  if (find >= 2) {
#ifndef NOT_PRINT_RELATION
    printf_relation(factor, a, t, V, outstd, assert_facto, c, smooth);
#endif // NOT_PRINT_RELATION
    (* nb_rel_found)++;

    if (gal == 6) {
      if (gal_version == 1) {
#ifndef NOT_PRINT_RELATION
        printf_relation_galois_6_1(factor, a, t, V, outstd, assert_facto, c,
            smooth, f, gal_norm_denom);
#endif // NOT_PRINT_RELATION
      }
      (* nb_rel_gal) += 5;
    }
  }

  mpz_clear(res);
  if (main >= 0 && find == 0) {
    factor_clear(factor[main]);
  } else {
    for (unsigned int i = 0; i < V; i++) {
      factor_clear(factor[i]);
    }
  }
  free(factor);
  free(smooth);

#ifdef ASSERT_FACTO
  mpz_clear(res_tmp);
  free(assert_facto);
#endif
}

/*
 * Return the sum of the element in index. index has size V.
 *
 * index: array of index.
 * V: number of element, number of number fields.
 */
static uint64_t sum_index(uint64_t * index, unsigned int V, int main)
{
  ASSERT(main >= -1);
  ASSERT(V >= 2);

  if (main != -1) {
    return index[main];
  }
  uint64_t sum = 0;
  for (unsigned int i = 0; i < V; i++) {
    sum += index[i];
  }
  return sum;
}

/*
 * Find in all the indices which one has the smallest value.
 *
 *
 *
 */
static unsigned int find_indices(unsigned int ** L,
    uint64_array_t * indices, uint64_t * index, unsigned int V,
    uint64_t max_indices, unsigned int * pos)
{
  unsigned int size = 0;
  uint64_t min = max_indices;
  for (unsigned int i = 0; i < V; i++) {
    if (indices[i]->length != 0 && index[i] < indices[i]->length) {
      if (indices[i]->array[index[i]] < min) {
        min = indices[i]->array[index[i]];
        * pos = i;
        (*L)[i] = 1;
        size = 1;
        for (unsigned int j = 0; j < i; j++) {
          (*L)[i] = 0;
        }
      } else if (min == indices[i]->array[index[i]]) {
        (*L)[i] = 1;
        size++;
      }
    }
  }
  return size;
}

static unsigned int find_indices_main(unsigned int ** L,
    uint64_array_t * indices, uint64_t * index, unsigned int V, int main)
{
  ASSERT(main >= 0);

  (*L)[main] = 1;
  unsigned int size = 1;
  uint64_t target = indices[main]->array[index[main]];
  for (unsigned int i = 0; i < V; i++) {
    if (i != (unsigned int) main && indices[i]->length != 0 && index[i] <
        indices[i]->length) {
      while (index[i] < indices[i]->length - 1 &&
          target > indices[i]->array[index[i]]) {
        index[i] = index[i] + 1;
      }
      if (target == indices[i]->array[index[i]]) {
        (*L)[i] = 1;
        size++;
      } else {
        (*L)[i] = 0;
      }
    }
  }
  return size;
}

/*
 * For
 */
static void find_relation(uint64_array_t * indices, uint64_t * index,
    uint64_t number_element, mat_Z_srcptr matrix, const mpz_poly * f,
    sieving_bound_srcptr H, unsigned int V, int main, uint64_t max_indices,
    facul_aux_data *data, unsigned int * nb_rel_found,
    ideal_spq_srcptr special_q, unsigned int q_side, FILE * outstd,
    MAYBE_UNUSED unsigned int * number_factorisation, unsigned int gal,
    unsigned int * nb_rel_gal, unsigned int gal_version,
    factor_t * gal_norm_denom)
{
  unsigned int * L = (unsigned int *) malloc(V * sizeof(unsigned int));
  unsigned int size = 0;
  unsigned int pos = 0;
  if (main == -1) {
    size = find_indices(&L, indices, index, V, max_indices, &pos);
  } else {
    pos = main;
    size = find_indices_main(&L, indices, index, V, main);
  }

  ASSERT(size >= 1);

  if (size >= 2) {
    mpz_vector_t c;
    mpz_t gcd;
    mpz_init(gcd);
    mpz_vector_init(c, H->t);
    array_index_mpz_vector(c, indices[pos]->array[index[pos]], H,
        number_element);

    mpz_poly a;
    mpz_poly_init(a, 0);
    mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, c);
    mpz_poly_content(gcd, a);

    //a must be irreducible.
    if (mpz_cmp_ui(gcd, 1) == 0 && a->deg > 0 &&
        mpz_cmp_ui(mpz_poly_lc_const(a), 0) > 0 && is_irreducible(a)) {

      good_polynomial(a, f, L, H->t, V, main, data, nb_rel_found,
          special_q, q_side, size, outstd, number_factorisation, c, gal,
          nb_rel_gal, gal_version, gal_norm_denom);
    }

    mpz_poly_clear(a);

    mpz_clear(gcd);
    mpz_vector_clear(c);

    for (unsigned int i = 0; i < V; i++) {
      if (L[i] == 1) {
        index[i] = index[i] + 1;
      }
    }
  } else {
    index[pos] = index[pos] + 1;
  }

  free(L);
}

unsigned int find_relations(uint64_array_t * indices, uint64_t number_element,
    unsigned int * lpb, mat_Z_srcptr matrix, const mpz_poly * f,
    sieving_bound_srcptr H, unsigned int V, ideal_spq_srcptr special_q,
    unsigned int q_side, int main, FILE * outstd, unsigned int gal,
    unsigned int gal_version, factor_t * gal_norm_denom)
{
  //index[i] is the current index of indices[i].
  uint64_t * index = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  uint64_t length_tot = 0;
  //Maximum of all the indices.
  uint64_t max_indices = 0;
  for (unsigned int i = 0; i < V; i++) {
    index[i] = 0;
    if (i == (unsigned int) main) {
      ASSERT(main >= 0);

      length_tot = indices[main]->length;
    } else if (indices[i]->length != 0 && main == -1) {
      length_tot += indices[i]->length;
      if (max_indices < indices[i]->array[indices[i]->length - 1]) {
        max_indices = indices[i]->array[indices[i]->length - 1];
      }
    }
  }

  //  prepare cofactorization strategy
  facul_aux_data * data;
  data = (facul_aux_data *)malloc(V*sizeof(facul_aux_data));
  ASSERT_ALWAYS(data != NULL);
  for (unsigned int i = 0; i < V; ++i) {
    unsigned int B = 1000; // FIXME: should be a parameter.
    /*unsigned int B = (unsigned int) (1 << (lpb[i] - 1));*/
    data[i].lpb = (unsigned long) lpb[i];
    data[i].fbb = B;
    data[i].BB = ((double)B) * ((double)B);
    data[i].BBB = ((double)B) * data[i].BB;
    data[i].methods = facul_make_default_strategy(nb_curves95(lpb[i]), 0);
  }

  MAYBE_UNUSED unsigned int * number_factorisation =
    (unsigned int *) malloc(sizeof(unsigned int) * 2);
  //Bad factorisation.
  number_factorisation[0] = 0;
  //Good factorisation.
  number_factorisation[1] = 0;
  unsigned int nb_rel_found = 0;
  unsigned int nb_rel_gal = 0;
  if (0 != length_tot) {
    while(sum_index(index, V, main) < length_tot) {
      find_relation(indices, index, number_element, matrix, f, H, V, main,
          max_indices, data, &nb_rel_found, special_q, q_side, outstd,
          number_factorisation, gal, &nb_rel_gal, gal_version, gal_norm_denom);
    }
  }
  if (nb_rel_found != 0) {
    if (gal == 6) {
      fprintf(outstd, "# Number of found relations: %u (with %u by Galois \
6).\n", nb_rel_found + nb_rel_gal, nb_rel_gal);
    } else {
      fprintf(outstd, "# Number of found relations: %u.\n", nb_rel_found);
    }
  } else {
    fprintf(outstd, "# No relations\n");
  }

#ifdef ASSERT_FACTO
  fprintf(outstd, "# Number of complete factorisations: %u.\n",
      number_factorisation[1]);
  fprintf(outstd, "# Number of incomplete factorisations: %u.\n",
      number_factorisation[0]);
#endif // ASSERT_FACTO
  free(number_factorisation);

  for (unsigned int i = 0; i < V; ++i)
    facul_clear_methods(data[i].methods);
  free(data);
  free(index);

  return nb_rel_found + nb_rel_gal;
}
