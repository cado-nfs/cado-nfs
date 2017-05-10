#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "portability.h"
#include "utils.h"

/********************** internal functions *****************************/

static const int ugly[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};


static uint64_t next_prime_of_powers_of_2[64] = { 0x2, 0x3, 0x5, 0xb, 0x11,
        0x25, 0x43, 0x83, 0x101, 0x209, 0x407, 0x805, 0x1003, 0x2011, 0x401b,
        0x8003, 0x10001, 0x2001d, 0x40003, 0x80015, 0x100007, 0x200011,
        0x40000f, 0x800009, 0x100002b, 0x2000023, 0x400000f, 0x800001d,
        0x10000003, 0x2000000b, 0x40000003, 0x8000000b, 0x10000000f,
        0x200000011, 0x400000019, 0x800000035, 0x100000001f, 0x2000000009,
        0x4000000007, 0x8000000017, 0x1000000000f, 0x2000000001b, 0x4000000000f,
        0x8000000001d, 0x100000000007, 0x20000000003b, 0x40000000000f,
        0x800000000005, 0x1000000000015, 0x2000000000045, 0x4000000000037,
        0x8000000000015, 0x10000000000015, 0x20000000000005, 0x4000000000009f,
        0x80000000000003, 0x100000000000051, 0x200000000000009,
        0x400000000000045, 0x800000000000083, 0x1000000000000021,
        0x200000000000000f, 0x4000000000000087, 0x800000000000001d };

static uint64_t previous_prime_of_powers_of_2[65] = { 0x0, 0x0, 0x3, 0x7, 0xd,
        0x1f, 0x3d, 0x7f, 0xfb, 0x1fd, 0x3fd, 0x7f7, 0xffd, 0x1fff, 0x3ffd,
        0x7fed, 0xfff1, 0x1ffff, 0x3fffb, 0x7ffff, 0xffffd, 0x1ffff7, 0x3ffffd,
        0x7ffff1, 0xfffffd, 0x1ffffd9, 0x3fffffb, 0x7ffffd9, 0xfffffc7,
        0x1ffffffd, 0x3fffffdd, 0x7fffffff, 0xfffffffb, 0x1fffffff7,
        0x3ffffffd7, 0x7ffffffe1, 0xffffffffb, 0x1fffffffe7, 0x3fffffffd3,
        0x7ffffffff9, 0xffffffffa9, 0x1ffffffffeb, 0x3fffffffff5, 0x7ffffffffc7,
        0xfffffffffef, 0x1fffffffffc9, 0x3fffffffffeb, 0x7fffffffff8d,
        0xffffffffffc5, 0x1ffffffffffaf, 0x3ffffffffffe5, 0x7ffffffffff7f,
        0xfffffffffffd1, 0x1fffffffffff91, 0x3fffffffffffdf, 0x7fffffffffffc9,
        0xfffffffffffffb, 0x1fffffffffffff3, 0x3ffffffffffffe5, 0x7ffffffffffffc9,
        0xfffffffffffffa3, 0x1fffffffffffffff, 0x3fffffffffffffc7,
        0x7fffffffffffffe7, 0xffffffffffffffc5 };



/* Skip line beginning with '#'. Return 0 if fgets return NULL.*/
static inline size_t
get_one_line (FILE *f, char *s)
{
  char *rets;
  size_t n;
  do
  {
    rets = fgets(s, RENUMBER_MAXLINE, f);
    if (rets == NULL)
    {
      n = 0;
      break;
    }
    else if (rets[0] != '#')
    {
      n = strnlen(s, RENUMBER_MAXLINE);
      ASSERT_ALWAYS(n != RENUMBER_MAXLINE);
      break;
    }
    // else we skip the line
  } while (1);
  return n;
}

// p,r values can be in decimal (in the badideals file) or in hexa (in
// the renumber file). This is the purpose of the BASE parameter.
static void
parse_one_line_bad_ideals(struct bad_ideals_s * bad, const char * str, int k,
    int BASE)
{
  int t;
  const char *ptr = str;
  p_r_values_t p = 0, r = 0;
  int side = 0, nb = 0;

  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    p = (p * BASE) + t;
  ASSERT_ALWAYS(*ptr == ',');
  ptr++;

  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    r = (r * BASE) + t;
  ASSERT_ALWAYS(*ptr == ':');
  ptr++;

  side = ugly[(unsigned char) *ptr];
  ptr++;
  ASSERT_ALWAYS(side == 0 || side == 1); /* FIXME: for MNFS */
  ASSERT_ALWAYS(*ptr == ':');

  ptr++;
  while (isspace(ptr[0]))
    ptr++;
  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    nb = (nb * 10) + t;
  ASSERT_ALWAYS(*ptr == '\n');
  ASSERT_ALWAYS(nb <= RENUMBER_MAX_ABOVE_BADIDEALS);

  bad->p[k] = p;
  bad->r[k] = r;
  bad->nb[k] = nb;
  bad->side[k] = side;
}

static inline p_r_values_t
parse_one_line (char * str)
{
  p_r_values_t v = 0;
  int t;
  char *p;

  for (p = str; (t=ugly[(unsigned char) *p]) >= 0; p++)
    v = (v << 4) + t;

  ASSERT_ALWAYS(*p == '\n');
  return v;
}

static void
parse_bad_ideals_file (FILE *badidealsfile, renumber_ptr renum)
{
  renum->bad_ideals.n = 0;
  char s[RENUMBER_MAXLINE];

  /* first count the number of lines, i.e., the number of bad ideals */
  while (get_one_line (badidealsfile, s) != 0)
    renum->bad_ideals.n++;
  ASSERT_ALWAYS (feof (badidealsfile));

  /* Allocate memory to store badideals information */
  size_t badideals_pr_size = renum->bad_ideals.n * sizeof (p_r_values_t);
  size_t badideals_int_size = renum->bad_ideals.n * sizeof (int);
  renum->bad_ideals.p = (p_r_values_t *) malloc (badideals_pr_size);
  renum->bad_ideals.r = (p_r_values_t *) malloc (badideals_pr_size);
  renum->bad_ideals.side = (int *) malloc (badideals_int_size);
  renum->bad_ideals.nb = (int *) malloc (badideals_int_size);
  ASSERT_ALWAYS (renum->bad_ideals.p != NULL);
  ASSERT_ALWAYS (renum->bad_ideals.r != NULL);
  ASSERT_ALWAYS (renum->bad_ideals.nb != NULL);
  ASSERT_ALWAYS (renum->bad_ideals.side != NULL);

  /* Finally, read the information and store it */
  rewind (badidealsfile);
  int k = 0;
  renum->bad_ideals.max_p = 0;
  while (get_one_line (badidealsfile, s) != 0)
  {
    parse_one_line_bad_ideals (&(renum->bad_ideals), s, k, 10);
    renum->size += renum->bad_ideals.nb[k];
    renum->bad_ideals.max_p=MAX(renum->bad_ideals.max_p,renum->bad_ideals.p[k]);
    k++;
  }
  ASSERT_ALWAYS (k == renum->bad_ideals.n);
  ASSERT_ALWAYS (feof (badidealsfile));
}

static void
print_info (FILE * f, renumber_srcptr r, int after_reading)
{
  char pre[9] = "# INFO: ";
  fprintf (f, "# Information on renumber struct:\n%ssizeof(p_r_values_t) = %zu\n"
              "%snb_bits = %" PRIu8 "\n%snumber of polynomials = %u\n", pre,
              sizeof(p_r_values_t), pre, r->nb_bits, pre, r->nb_polys);
  if (r->rat == -1)
    fprintf (f, "%sThere is no rational side\n", pre);
  else
    fprintf (f, "%sPolynomial on side %d is rational\n", pre, r->rat);
  fprintf (f, "%s#badideals = %d [max_p = %" PRpr "]\n", pre, r->bad_ideals.n,
                                                         r->bad_ideals.max_p);
  fprintf (f, "%s#additional columns = %u\n", pre, r->naddcols);
  if (r->nonmonic)
  {
    fprintf (f, "%sNon monic polynomial on side:", pre);
    for (uint64_t t = r->nonmonic, s = 0; t != 0; t>>=1, s++)
      if (t & ((uint64_t) 1))
        fprintf (f, " %" PRIu64 "", s);
    fprintf (f, "\n");
  }
  else
    fprintf (f, "%sAll polynomials are monic\n", pre);
  for (unsigned int i = 0; i < r->nb_polys; i++)
    fprintf (f, "%slpb%d = %lu\n", pre, i, r->lpb[i]);

  if (after_reading) /* there is more info to print*/
  {
    fprintf (f, "%ssize = %" PRIu64 "\n%ssmallest prime not cached = "
             "0x%" PRpr " at index 0x%" PRid "\n", pre, r->size, pre,
             r->smallest_prime_not_cached, r->index_smallest_prime_not_cached);
    for (unsigned int i = 0; i < r->nb_polys; i++)
      fprintf (f, "%sbiggest prime below lbp%d is 0x%" PRpr " at index "
                  "0x%" PRid "\n", pre, i, r->biggest_prime_below_lpb[i],
                  r->index_biggest_prime_below_lpb[i]);
  }
  fflush (stdout);
}

/* sort in decreasing order. Fastest for ~ < 15 values in r[] vs qsort */
static inline void
renumber_sort_ul (unsigned long *r, size_t n)
{
  unsigned long rmin;

  if (UNLIKELY (n < 2))
    return;

  if (UNLIKELY (n == 2)) {
    if (r[0] < r[1]) {
      rmin = r[0];
      r[0] = r[1];
      r[1] = rmin;
    }
    return;
  }

  for (size_t i = n; --i;) {
    size_t min = i;
    rmin = r[min];
    for (size_t j = i; j--;) {
      unsigned long rj = r[j];
      if (UNLIKELY (rj < rmin)) {
	min = j;
	rmin = rj;
      }
    }
    if (LIKELY (min != i)) {
      r[min] = r[i];
      r[i] = rmin;
    }
  }
}

/* return zero if no roots mod p, else non-zero */
static int
get_largest_root_mod_p (p_r_values_t *r, mpz_poly_srcptr f, p_r_values_t p)
{
  int deg = f->deg;
  /* If there is a projective root, it is the largest (r = p by convention) */
  if (mpz_divisible_ui_p (f->coeff[deg], p))
  {
    *r = p;
    return 1;
  }

  unsigned long roots[deg];
  size_t k = (size_t) mpz_poly_roots_ulong (roots, (mpz_poly_ptr) f, p);
  if (k)
  {
    unsigned long max = roots[--k];
    while (k--)
      if (UNLIKELY (roots[k] > max)) max = roots[k];
    *r = max;
    return 2;
  }

  return 0;
}

static void
renumber_write_first_line (renumber_ptr renum)
{
  fprintf (renum->file, "%" PRIu8 " %d %d %d %" PRIx64 " %u", renum->nb_bits,
           renum->rat, renum->bad_ideals.n, renum->naddcols, renum->nonmonic,
           renum->nb_polys);
  for (unsigned int i = 0; i < renum->nb_polys; i++)
    fprintf (renum->file, " %lu", renum->lpb[i]);
  fprintf (renum->file, "\n");
}

static void
renumber_read_first_line (renumber_ptr renum)
{
  int ret;
  ret = fscanf (renum->file, "%" SCNu8 " %d %d %" SCNu8 " %" SCNx64 " %u",
                &(renum->nb_bits), &(renum->rat), &(renum->bad_ideals.n),
                &(renum->naddcols), & (renum->nonmonic), &(renum->nb_polys));
  ASSERT_ALWAYS (ret == 6);
  ASSERT_ALWAYS (renum->nb_polys >= 1);
  ASSERT_ALWAYS (-1 <= renum->rat && renum->rat < (int) renum->nb_polys);
  ASSERT_ALWAYS (renum->nb_bits <= 8 * sizeof(p_r_values_t));
  ASSERT_ALWAYS (renum->nb_bits == 32 || renum->nb_bits == 64);
  renum->lpb = (unsigned long *) malloc (renum->nb_polys*sizeof(unsigned long));
  ASSERT_ALWAYS (renum->lpb != NULL);

  for (unsigned int i = 0; i < renum->nb_polys; i++)
  {
    if (i == renum->nb_polys - 1)
      ret = fscanf (renum->file, " %lu\n", &(renum->lpb[i]));
    else
      ret = fscanf (renum->file, " %lu", &(renum->lpb[i]));
    ASSERT_ALWAYS (ret == 1);
  }
}

static inline p_r_values_t
compute_vp_from_p (renumber_srcptr tab, p_r_values_t p)
{
  if (tab->nb_polys == 2)
  {
    if (tab->rat >= 0) /* If there is a rational side */
      return p + 1;
    else
      return (p << 1) + 1;
  }
  else
  {
    if (tab->rat >= 0) /* If there is a rational side */
      return (tab->nb_polys - 1) * (p + 1);
    else
      return tab->nb_polys * (p + 1) - 1;
  }
}

static inline p_r_values_t
compute_vr_from_p_r (renumber_srcptr tab, p_r_values_t p, p_r_values_t r, int side)
{
  if (tab->nb_polys == 2)
  {
    if (tab->rat >= 0) /* If there is a rational side */
      return (side == tab->rat) ? (p + 1) : r;
    else
      return (side == 1) ? (p + 1 + r) : r;
  }
  else
  {
    if (tab->rat >= 0) /* If there is a rational side */
    {
      if (side == tab->rat)
        return (tab->nb_polys - 1) * (p + 1);
      else if (side < tab->rat)
        return side * (p + 1) + r;
      else
        return (side-1) * (p + 1) + r;
    }
    else
      return side * (p + 1) + r;
  }
}

/* Inverse function of compute_vp_from_p */
static inline p_r_values_t
compute_p_from_vp (renumber_srcptr tab, p_r_values_t vp)
{
  if (tab->nb_polys == 2)
  {
    if (tab->rat >= 0) /* One alg and one rat side -> p is v-1 */
      return (vp - 1);
    else               /* Two alg sides -> p is (v-1)/2 */
      return (vp >> 1);
  }
  else
  {
    if (tab->rat >= 0) /* One rat side and (nb_polys-1) alg side */
      return (vp/(tab->nb_polys - 1)) - 1;
    else               /* nb_polys alg sides -> p = (vp+1)/nb_polys - 1 */
      return ((vp + 1)/tab->nb_polys) - 1;
  }
}

static inline void
compute_r_side_from_p_vr (p_r_values_t *r, int *side, renumber_srcptr tab,
                          p_r_values_t p, p_r_values_t vr)
{
  *side = 0;
  *r = vr;
  while (*r > p)
  {
    *r -= (p + 1);
    *side += 1;
  }

  if (tab->rat >= 0) /* If there is a rational side */
    if (*side >= tab->rat)
      *side += 1;
}

/*********************** End internal functions  ******************************/

void
renumber_init_for_reading (renumber_ptr renumber_info)
{
  memset(renumber_info, 0, sizeof(renumber_t));
  /* Will be set later, by renumber_read_table, with the values of the first
   * line of the renumber file */
}

/* nb_polys contains the number of polynomials (must be at least 2).
 * rat contains the rational side or -1 if no rational side.
 * add_full_col is non-zero if we need to add a column of 1 in the matrix, 0
 * otherwise (for factorization, always 0, for DL 1 if one of the polynomials is
 * not monic). */
void
renumber_init_for_writing (renumber_ptr renumber_info, unsigned int nb_polys,
                           int rat, int lcideals, uint64_t nonmonic,
                           unsigned long *lpb)
{
  memset(renumber_info, 0, sizeof(renumber_t));

  ASSERT_ALWAYS (nb_polys >= 1);
  ASSERT_ALWAYS (-1 <= rat && rat < (int) nb_polys);
  ASSERT_ALWAYS (lpb != NULL);
  renumber_info->nb_polys = nb_polys;
  renumber_info->rat = rat;
  
  /* Set lpb table */
  size_t size_ul = nb_polys * sizeof (unsigned long);
  renumber_info->lpb = (unsigned long *) malloc (size_ul);
  ASSERT_ALWAYS (renumber_info->lpb != NULL);
  memcpy (renumber_info->lpb, lpb, size_ul);

  /* Set max_lpb */
  renumber_info->max_lpb = lpb[0]; /* There are at least 1 sides */
  for (unsigned int i = 1; i < nb_polys; i++)
    renumber_info->max_lpb = MAX(renumber_info->max_lpb, lpb[i]);

  /* Set max_nb_bits */
  uint64_t max_p = previous_prime_of_powers_of_2[renumber_info->max_lpb];
  uint64_t max_vp; /* same as compute_vp_from_p but with uint64_t */
  if (rat >= 0)
    max_vp = (nb_polys - 1) * (max_p + 1);
  else
    max_vp = nb_polys * (max_p + 1) - 1;

  if (nbits (max_vp) <= 32)
    renumber_info->nb_bits = 32;
  else
    renumber_info->nb_bits = 64;
  if (renumber_info->nb_bits > 8 * sizeof(p_r_values_t))
    {
      fprintf (stderr, "Error, p_r_values_t is too small to store ideals, "
               "recompile with FLAGS_SIZE=\"-D__SIZEOF_P_R_VALUES__=8\"\n");
      exit (1);
    }

  /* Compute the number of additional columns needed in the renumbering table
   * due to the non monic polynomials.
   * If lcideals == 0 (as it should be for factorization):
   *    --> no column is added (not needed).
   * If lcideals != 0:
   *    if nb_polys == 2:
   *        if nonmonic == 0:
   *          --> no column is added (all polynomials are monic).
   *        if nonmonic != 0:
   *          --> 1 column is added.
   *    if nb_polys != 2:
   *          --> add as many columns as the number of set bits in nonmonic
   */
  renumber_info->nonmonic = nonmonic;
  if (lcideals == 0)
    renumber_info->naddcols = 0;
  else
  {
    if (nb_polys == 2)
      renumber_info->naddcols = (nonmonic != 0);
    else
      for (renumber_info->naddcols = 0; nonmonic ; nonmonic>>=1)
        renumber_info->naddcols += nonmonic & ((uint64_t) 1);
  }
}

void
renumber_clear (renumber_ptr renumber_info)
{
  if (renumber_info->table != NULL)
    free (renumber_info->table);
  if (renumber_info->cached != NULL)
    free (renumber_info->cached);

  renumber_info->table = NULL;
  renumber_info->cached = NULL;

 if (renumber_info->lpb != NULL)
  free (renumber_info->lpb);
 if (renumber_info->biggest_prime_below_lpb != NULL)
  free (renumber_info->biggest_prime_below_lpb);
 if (renumber_info->index_biggest_prime_below_lpb != NULL)
  free (renumber_info->index_biggest_prime_below_lpb);

 renumber_info->lpb = NULL;
 renumber_info->biggest_prime_below_lpb = NULL;
 renumber_info->index_biggest_prime_below_lpb = NULL;

  if (renumber_info->bad_ideals.p != NULL)
    free (renumber_info->bad_ideals.p);
  if (renumber_info->bad_ideals.r != NULL)
    free (renumber_info->bad_ideals.r);
  if (renumber_info->bad_ideals.nb != NULL)
    free (renumber_info->bad_ideals.nb);
  if (renumber_info->bad_ideals.side != NULL)
    free (renumber_info->bad_ideals.side);

  renumber_info->bad_ideals.r = NULL;
  renumber_info->bad_ideals.p = NULL;
  renumber_info->bad_ideals.nb = NULL;
  renumber_info->bad_ideals.side = NULL;
}

/* The renumber_t struct _must_ have been initialized before
 * poly = NULL is accepted. It will not print the polynomials on the file */
void
renumber_write_open (renumber_ptr tab, const char *tablefile, 
		     const char *badidealsfile, cado_poly poly)
{
  printf ("# Opening %s to write the renumbering table\n", tablefile);
  fflush (stdout);
  tab->file = fopen_maybe_compressed (tablefile, "w");
  ASSERT_ALWAYS(tab->file != NULL);
  FILE *fbad = NULL;
  if (badidealsfile != NULL)
  {
    printf ("# Opening %s to read the bad ideals\n", badidealsfile);
    fbad = fopen(badidealsfile, "r"); /* never compressed, always small */
    ASSERT_ALWAYS (fbad != NULL);
  }

  /* additional columns are always at the beginning of the renumbering table. */
  tab->size = tab->naddcols;

  /* Read bad ideals files */
  if (badidealsfile != NULL)
    parse_bad_ideals_file (fbad, tab); /* update size et bad_ideals.n */

  /* Write the first line */
  renumber_write_first_line (tab);

  /* Print info on stdout (~ what is written on the first line of the file) */
  print_info (stdout, tab, 0);

  /* Write the two polynomials on a line beginning by #, if given */
  if (poly != NULL)
  {
    ASSERT_ALWAYS (poly->nb_polys == (int)tab->nb_polys);
    for (int i = 0; i < poly->nb_polys; i++)
    {
      fprintf (tab->file, "# pol%d: ", i);
      mpz_poly_fprintf_coeffs (tab->file, poly->pols[i], ',');
    }
  }

  /* Write first the bad ideals information at the beginning of file */
  if (badidealsfile != NULL)
  {
    for (int i = 0; i < tab->bad_ideals.n; ++i) {
      fprintf(tab->file, "%" PRpr ",%" PRpr ":%d: %d\n",
          tab->bad_ideals.p[i],
          tab->bad_ideals.r[i],
          tab->bad_ideals.side[i],
          tab->bad_ideals.nb[i]);
    }
  }

  if (fbad != NULL)
    fclose (fbad);
}

void
renumber_write_close (renumber_ptr tab, const char *tablefile)
{
  fclose_maybe_compressed (tab->file, tablefile);
}


/* The renumber_t struct _must_ have been initialized before */
void
renumber_read_table (renumber_ptr tab, const char * filename)
{
  char s[RENUMBER_MAXLINE];
  size_t bytes_read = 0, bytes_line;
  stats_data_t infostats;  /* for displaying progress */

  /* open file for reading */
  printf ("# Opening %s to read the renumbering table\n", filename);
  fflush (stdout);
  tab->file = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK (tab->file == NULL, "Cannot open file for reading");

  /* read size of renumbering table */
  renumber_read_first_line (tab);

  /* Allocating memory */
  size_t badideals_pr_size = tab->bad_ideals.n * sizeof (p_r_values_t);
  size_t badideals_int_size = tab->bad_ideals.n * sizeof (int);
  size_t cached_table_size = (2 << MAX_LOG_CACHED) * sizeof (index_t);
  size_t primes_nb_polys_size = tab->nb_polys * sizeof (p_r_values_t);
  size_t index_nb_polys_size = tab->nb_polys * sizeof (index_t);

  /* Do not know the size yet. Reallocating while reading.
     We assume that RENUMBER_DEFAULT_SIZE is enough to hold at least the bad
     ideals and the added column (if add_full_col is set) */
  uint64_t allocated = RENUMBER_DEFAULT_SIZE;
  size_t default_size = RENUMBER_DEFAULT_SIZE * sizeof (p_r_values_t);

  tab->table           = (p_r_values_t *)  malloc (default_size);
  tab->cached          = (index_t *)       malloc (cached_table_size);
  tab->biggest_prime_below_lpb = (p_r_values_t *) malloc (primes_nb_polys_size);
  tab->index_biggest_prime_below_lpb = (index_t *) malloc (index_nb_polys_size);
  tab->bad_ideals.p    = (p_r_values_t *)  malloc (badideals_pr_size);
  tab->bad_ideals.r    = (p_r_values_t *)  malloc (badideals_pr_size);
  tab->bad_ideals.side = (int *)           malloc (badideals_int_size);
  tab->bad_ideals.nb   = (int *)           malloc (badideals_int_size);

  ASSERT_ALWAYS (tab->table != NULL);
  ASSERT_ALWAYS (tab->cached != NULL);
  ASSERT_ALWAYS (tab->biggest_prime_below_lpb != NULL);
  ASSERT_ALWAYS (tab->index_biggest_prime_below_lpb != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.p != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.r != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.nb != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.side != NULL);

  memset (tab->cached, 0, cached_table_size);

  /* additional columns are always at the beginning of the renumbering table. */
  tab->size = tab->naddcols;
  for (uint8_t i = 0; i < tab->naddcols; i++)
    tab->table[i] = RENUMBER_SPECIAL_VALUE;

  /* Reading the bad ideals at the top of the renumbering file */
  tab->bad_ideals.max_p = 0;
  for (int k = 0; k < tab->bad_ideals.n; k++)
  {
    bytes_read += get_one_line(tab->file, s);
    parse_one_line_bad_ideals (&tab->bad_ideals, s, k, 16);
    tab->bad_ideals.max_p = MAX (tab->bad_ideals.max_p, tab->bad_ideals.p[k]);
    for (int j = 0; j < tab->bad_ideals.nb[k]; j++)
    {
      tab->table[tab->size] = RENUMBER_SPECIAL_VALUE;
      tab->size++;
    }
  }

  p_r_values_t prime_cache_limit = next_prime_of_powers_of_2[MAX_LOG_CACHED];
  p_r_values_t * expected_biggest_prime_lpb =
                                 (p_r_values_t *) malloc (primes_nb_polys_size);
  for (unsigned int i = 0; i < tab->nb_polys; i++)
    expected_biggest_prime_lpb[i] = previous_prime_of_powers_of_2[tab->lpb[i]];
  int has_smallest = 0;

  /* Reading the renumbering table */
  stats_init (infostats, stdout, &(tab->size), 24, "Read", "elements", "",
              "elts");

  while ((bytes_line = get_one_line(tab->file, s)) > 0)
  {
    bytes_read += bytes_line;
    if (tab->size >= allocated) /* Not enough space, reallocated tab->table */
    {
      allocated += RENUMBER_DEFAULT_SIZE;
      size_t new_size = allocated * sizeof (p_r_values_t);
      tab->table = (p_r_values_t *) realloc (tab->table, new_size);
      ASSERT_ALWAYS (tab->table != NULL);
    }
    tab->table[tab->size] = parse_one_line(s);

    if (tab->size == 0 || tab->table[tab->size-1] == RENUMBER_SPECIAL_VALUE
                       || tab->table[tab->size] > tab->table[tab->size-1])
    {
      /* We just switch to a new prime in the renumbering table, see if we need
       * to cache it (we cached primes below 2^MAX_LOG_CACHED)
       */
      p_r_values_t p = compute_p_from_vp (tab, tab->table[tab->size]);
      if (p < prime_cache_limit) /* p < 2^MAX_LOG_CACHED */
        tab->cached[p] = tab->size;
      else if (!has_smallest)
      {
        has_smallest = 1;
        tab->index_smallest_prime_not_cached = tab->size;
        tab->smallest_prime_not_cached = p;
      }

      for (unsigned int i = 0; i < tab->nb_polys; i++)
      {
        if (p <= expected_biggest_prime_lpb[i])
        {
          tab->index_biggest_prime_below_lpb[i] = tab->size;
          tab->biggest_prime_below_lpb[i] = p;
        }
      }
    }
    tab->size++;

    if (stats_test_progress(infostats))
      stats_print_progress (infostats, tab->size, 0, bytes_read, 0);
  }

  if (!has_smallest) /* Every prime is cached. */
  {
    tab->index_smallest_prime_not_cached = tab->size;
    tab->smallest_prime_not_cached =
           next_prime_of_powers_of_2[MAX(tab->lpb[0],tab->lpb[1])];
  }

  stats_print_progress (infostats, tab->size, 0, bytes_read, 1);
  size_t final_size = tab->size * sizeof (p_r_values_t);
  tab->table = (p_r_values_t *) realloc (tab->table, final_size);


  ASSERT_ALWAYS (feof (tab->file));

  print_info (stdout, tab, 1);

  free (expected_biggest_prime_lpb);
  fclose_maybe_compressed (tab->file, filename);
}

int renumber_is_bad (int *nb, index_t *first, renumber_srcptr rn, p_r_values_t p,
                     p_r_values_t r, int side)
{
  /* bad ideals start after the 'naddcols' additional columns in the table. */
  if (first != NULL)
    *first = rn->naddcols;
  int bad = 0;
  for (int i = 0; i < rn->bad_ideals.n; ++i)
  {
    if (p == rn->bad_ideals.p[i] && r == rn->bad_ideals.r[i]
                                 && side == rn->bad_ideals.side[i])
    {
      if (nb) *nb = rn->bad_ideals.nb[i];
      bad = 1;
      break;
    }
    else
      if (first) *first += rn->bad_ideals.nb[i];
  }
  return bad;
}

/* Return 1 if i does not correspond to a bad ideal, else return 0 */
int
renumber_badideal_get_p_r_below (renumber_srcptr renumber_info, p_r_values_t *p,
                                 p_r_values_t * r, int *side, index_t i)
{
  index_t bad = renumber_info->naddcols;
  if (i < bad)
    return 1; /* i is not a bad ideal, it is a additional column */
  int k;
  for (k = 0; k < renumber_info->bad_ideals.n; k++)
  {
    if (i < bad + renumber_info->bad_ideals.nb[k])
      break;
    bad += renumber_info->bad_ideals.nb[k];
  }
  if (k == renumber_info->bad_ideals.n)
    return 1; /* i is not a bad ideal, it is a "normal" ideal */
  *p = renumber_info->bad_ideals.p[k];
  *r = renumber_info->bad_ideals.r[k];
  *side = renumber_info->bad_ideals.side[k];
  return 0;
}

int
renumber_is_additional_column (renumber_srcptr tab, index_t h)
{
  return (h < tab->naddcols);
}

/* This function writes in a string the part of the renumbering table
 * corresponding to the prime p. It returns the number of chars that were
 * written.
 * /!\ roots[i] are sorted and are modified by this function.
 * /!\ This function assumes that the string out has enough space.
 */
size_t
renumber_write_buffer_p (char *out, renumber_ptr tab, unsigned long p,
                         unsigned long roots[][MAX_DEGREE], int nroots[])
{
  size_t n = 0;
  /* The root on rat side becomes vp. If there is no rat side or not root on
   * rat side for this prime (i.e., lpb < p), then the largest root becomes vp.
   */
  p_r_values_t vp = compute_vp_from_p (tab, (p_r_values_t) p);
  unsigned int replace_first = 0;

  for (unsigned int i = 0; i < tab->nb_polys; i++)
    renumber_sort_ul (roots[i], nroots[i]);

  /* Only one poly */
  if (tab->nb_polys == 1)
  {
    roots[0][0] = vp;
    for (int i = 0; i < nroots[0]; i++)
      n += sprintf (out + n, "%lx\n", roots[0][i]);
  }
  /* Two polys: alg and rat */
  else if (tab->nb_polys == 2 && tab->rat != -1)
  {
    int algside = 1-tab->rat;
    if (LIKELY(nroots[tab->rat])) /* There is at most 1 rational root. */
      n = sprintf (out, "%" PRpr "\n", vp);
    else
      roots[algside][0] = vp; /* No rational root */
    for (int i = 0; i < nroots[algside]; i++)
      n += sprintf (out + n, "%lx\n", roots[algside][i]);
  }
  /* Two alg polys */
  else if (tab->nb_polys == 2 && tab->rat == -1)
  {
    if (LIKELY (nroots[1]))
    {
      n = sprintf (out, "%" PRpr "\n", vp);
      for (int i = 1; i < nroots[1]; i++)
        n += sprintf (out + n, "%lx\n", roots[1][i] + p + 1);
    }
    else
      roots[0][0] = vp;

    for (int i = 0; i < nroots[0]; i++)
      n += sprintf (out + n, "%lx\n", roots[0][i]);
  }
  /* More than two polys (with or without rat side). */
  else
  {
    if (tab->rat == -1 || nroots[tab->rat] == 0) /* The largest root becomes vp */
      replace_first = 1;
    else
      n = sprintf (out, "%" PRpr "\n", vp);

    for (int i = tab->nb_polys - 1; i >= 0; i--)
    {
      if (i != tab->rat)
      {
        for (int j = 0; j < nroots[i]; j++)
        {
          if (UNLIKELY(replace_first))
          {
            n = sprintf (out, "%" PRpr "\n", vp);
            replace_first = 0;
          }
          else
          {
            p_r_values_t vr = compute_vr_from_p_r (tab, p, roots[i][j], i);
            n += sprintf (out + n, "%" PRpr "\n", vr);
          }
        }
      }
    }
  }
  return n;
}

void
renumber_write_p (renumber_ptr tab, unsigned long p,
                  unsigned long r[][MAX_DEGREE], int k[])
{
  size_t size_buffer;
  char buffer[2048];

  size_buffer = renumber_write_buffer_p (buffer, tab, p, r, k);

  fwrite ((void *) buffer, size_buffer, 1, tab->file);
  for (unsigned int i = 0; i < tab->nb_polys; i++)
    tab->size += k[i];
}

index_t
renumber_get_first_index_from_p(renumber_srcptr renumber_info,
    p_r_values_t p, int side)
{
  index_t i;
  p_r_values_t *tab = renumber_info->table;
  p_r_values_t vp = compute_vp_from_p (renumber_info,  p);

  /**************************************************************************/
  /* Search for i such that
        renumber_info->table[i] = vp
        this is the beginning of a decreasing sequence
  */

  /* For small value of p, the corresponding value of i is cached */
  if (p < renumber_info->smallest_prime_not_cached)
  {
    i = renumber_info->cached[p];
    if (UNLIKELY(tab[i] != vp))
    {
      /* There is a problem, most probably p is not prime. */
      fprintf(stderr, "Fatal error in %s at %s:%d\nError with the cached part of"
              " the renumbering table\n  p = 0x%" PRpr "\n  vp = 0x%" PRpr "\n"
              "  i = cached[p] = 0x%" PRid "\n  tab[i] = 0x%" PRpr "\n",
              __func__, __FILE__, __LINE__, p, vp, i, tab[i]);
      abort();
    }
  }
  /* p is not cached and below the lpb[side] */
  else if (p <= renumber_info->biggest_prime_below_lpb[side])
  {
    index_t max = renumber_info->index_biggest_prime_below_lpb[side];
    index_t min = renumber_info->index_smallest_prime_not_cached;

    if (UNLIKELY(p == renumber_info->biggest_prime_below_lpb[side]))
      i = max;
    else if (UNLIKELY(p == renumber_info->smallest_prime_not_cached))
      i = min;
    else /* We have to look for i such that tab[i] == vp between min and max. */
    {
      i = min + (max - min) / 2; /* avoids overflow when
                                    min + max >= UMAX(index_t) */

      /* Looking for vp: the values of vp are ordered in increasing order and are
        always at the beginning of a decreasing sequence */
      while (1)
      {
        index_t i_bak = i;
        while (i > 0 && tab[i-1] > tab[i])
          i--;

        if (i != min)
        {
          if (vp == tab[i])
            break;
          else if (vp < tab[i])
            max = i;
          else /* vp > tab[i] */
            min = i;
        }
        else
        {
          i = i_bak + 1;
          while (i < max && tab[i-1] > tab[i])
            i++;

          if (i != max)
          {
            if (vp == tab[i])
              break;
            else if (vp < tab[i])
              max = i;
            else /* vp > tab[i] */
              min = i;
          }
          else
          {
            /* prime p is not in the table => Fatal error */
            fprintf(stderr, "Fatal error in %s at %s:%d\nIdeals above p = 0x"
                            "%" PRpr " are not in the renumbering table\n"
                            "Maybe p is not prime ?\n", __func__, __FILE__,
                            __LINE__, p);
            abort();
          }
        }

        i = min + (max - min) / 2; /* avoids overflow (see above) */
      }
    }
  }
  else /* Error */
  {
    /* prime p is bigger than lpb[side] => Fatal error */
    fprintf(stderr, "Fatal error in %s at %s:%d\nIdeal (p, side) = (0x%" PRpr
                    ", %d) is bigger that large prime bound "
                    "2^%ld\n", __func__, __FILE__, __LINE__, p, side,
                    renumber_info->lpb[side]);
    abort();
  }
  return i;
}

index_t
renumber_get_random_index_from_p_side(renumber_srcptr renumber_info,
    p_r_values_t p, int side)
{
  // Get the first index corresponding to p
  index_t i = renumber_get_first_index_from_p(renumber_info, p, side);

  p_r_values_t *tab = renumber_info->table;
  // Some cases where there is no choice
  if (side == renumber_info->rat || i == renumber_info->size - 1
                                 || tab[i] <= tab[i+1]) {
    return i;
  }

  // Find list of valid indices
  index_t valid_i[MAX_DEGREE];
  int n = 0;
  do {
    p_r_values_t r;
    int s;
    compute_r_side_from_p_vr (&r, &s, renumber_info, p, tab[i]);
    if (s == side) {
      ASSERT_ALWAYS(n < MAX_DEGREE);
      valid_i[n] = i;
      n++;
    }
    i++;
  } while (i != renumber_info->size - 1 && tab[i-1] > tab[i]);
  ASSERT_ALWAYS(n != 0);

  long rnd = random()%n;
  return valid_i[rnd];
}


/* side is 0 if (p,r) corresponds to the left part in the relation,
 * side is 1 if (p,r) corresponds to the right part.
 * If side corresponds to the rational side (if it exists), the value of r is
 * meaningless. */
index_t
renumber_get_index_from_p_r (renumber_srcptr renumber_info, p_r_values_t p,
                             p_r_values_t r, int side)
{
  index_t i;
  p_r_values_t *tab = renumber_info->table;
  p_r_values_t vr; /* values of r as it is stored in the table*/

  vr = compute_vr_from_p_r (renumber_info,  p, r, side);

  i = renumber_get_first_index_from_p(renumber_info, p, side);

  /**************************************************************************/
  /* Now i points at the beginning of a decreasing sequence of values of vr */

  /* Return i in 4 cases:
    first case: an ideal on rational side always corresponds to the first
                element of a sequence
    second case: i is the last index of the table
    third case: the sequence contains only one value
    fourth case: next element of the sequence is too small to correspond to vr
  */
  if (side == renumber_info->rat || i == renumber_info->size - 1
                                 || tab[i] <= tab[i+1]
                                 || vr > tab[i+1])
    return i;
  else /* else we go through the sequence until we find vr */
  {
    while(i != renumber_info->size - 1 && tab[i] > tab[i+1])
    {
      i++;
      if (vr == tab[i])
        return i;
    }
    /* if we arrive here, there is a problem, the ideal was not found in the
       renumbering table */
    fprintf(stderr, "Fatal error in %s at %s:%d\nIdeal (p, r, side) = (0x%" PRpr
                    ", 0x%" PRpr ", %d) was not found on the renumbering "
                    "table\n", __func__, __FILE__, __LINE__, p, r, side);
    abort();
  }
}

/* This function assume that i does not correspond to a bad ideals or an
 * additional columns (i.e., tab[i] != RENUMBER_SPECIAL_VALUE). It will fail
 * if this assumption is not satisfied.
 */
void
renumber_get_p_r_from_index (renumber_srcptr renumber_info, p_r_values_t *p,
                             p_r_values_t * r, int *side, index_t i,
                             cado_poly pol)
{
  index_t j;
  p_r_values_t *tab = renumber_info->table;

  for (j = i; j > 0 && tab[j-1] > tab[j] && tab[j-1] != RENUMBER_SPECIAL_VALUE;)
    j--;

  *p = compute_p_from_vp (renumber_info, tab[j]);

  if (i != j)
    compute_r_side_from_p_vr (r, side, renumber_info, *p, tab[i]);
  else /* i = j */
  {
    /* If there is no rational side or if there is one but p is greater than the
       lpb on the rational side, r is the largest root of the last poly that
       has at least one root.
       If the second condition is evaluated, we knows that it means that
       renumber_info->rat >= 0.
    */
    if (renumber_info->rat == -1 ||
        *p > renumber_info->biggest_prime_below_lpb[renumber_info->rat])
    {
      *side = renumber_info->nb_polys - 1;
      while (*side >= 0 && (*p > renumber_info->biggest_prime_below_lpb[*side]
                            || !get_largest_root_mod_p (r, pol->pols[*side], *p)
                            || renumber_is_bad (NULL, NULL, renumber_info, *p,
                                                                    *r, *side)))
        (*side)--;
      ASSERT_ALWAYS (*side >= 0);
    }
    else
    {
      *r = RENUMBER_ROOT_ON_RAT_SIDE; /* root has no meaning on rat side */
      *side = renumber_info->rat;
    }
  }
}

int
renumber_get_side_from_index (renumber_srcptr renumber_info, index_t i,
                              cado_poly pol)
{
  p_r_values_t *tab = renumber_info->table;
  int side;

  if (tab[i] == RENUMBER_SPECIAL_VALUE)
  {
    if (renumber_is_additional_column (renumber_info, i))
    {
      side = 0;
      index_t index_add_col = 0;
      uint64_t b;
      for (b = renumber_info->nonmonic; b != 0; b>>=1, side++)
      {
        if (b & ((uint64_t) 1))
        {
          if (index_add_col == i)
            break;
          index_add_col++;
        }
      }
      ASSERT_ALWAYS (b != 0);
    }
    else /* i corresponds to a bad ideals. */
    {
      index_t bad = renumber_info->naddcols;
      int k;
      for (k = 0; k < renumber_info->bad_ideals.n; k++)
      {
        if (i < bad + renumber_info->bad_ideals.nb[k])
          break;
        bad += renumber_info->bad_ideals.nb[k];
      }
      ASSERT_ALWAYS (k < renumber_info->bad_ideals.n);
      side = renumber_info->bad_ideals.side[k];
    }
  }
  else
  {
    p_r_values_t p, r;
    renumber_get_p_r_from_index (renumber_info, &p, &r, &side, i, pol);
  }

  return side;
}

