/* arithmetic on polynomial with double-precision coefficients */
#include "cado.h" // IWYU pragma: keep

#ifdef HAVE_GLIBC
#include <features.h>   // __GLIBC_PREREQ
#endif

#include <cstring>
#include <cstdio> // FILE // IWYU pragma: keep
#include <cstdlib>
#include <climits>
#include <cfloat> /* for DBL_MAX */
#include <cctype> /* isspace */

#include <sstream>
#include <string>

#include <gmp.h>       // for mpz_t, mpz_clear, mpz_init, mpz_get_d, mpz_mul...
#include "mpz_poly.h"  // mpz_poly
#include "macros.h"     // for ASSERT, ASSERT_ALWAYS

#include "double_poly.h"
#include "utils_cxx.hpp"

// scan-headers: stop here

/* it's a bit nasty here. See
 * https://sourceware.org/bugzilla/show_bug.cgi?id=19439
 *
 */
#ifdef HAVE_GLIBC
#if __GLIBC_PREREQ(2, 23)
#include <cmath>
using std::isnan;
using std::isinf;
#else
/* not "the right way", but happens to work.  */
#include <math.h>
#endif
#else
/* this one is the right way */
#include <cmath>
using std::isnan;
using std::isinf;
#endif

#define STATIC_ANALYSIS_ASSERT_DATA_HEALTH(d) do {		     	   \
    ASSERT_FOR_STATIC_ANALYZER(((d)->alloc == 0) == ((d)->coeff == nullptr)); \
    ASSERT_FOR_STATIC_ANALYZER((d)->deg >= -1);                            \
    ASSERT_FOR_STATIC_ANALYZER((unsigned int) ((d)->deg+1) <= (d)->alloc);  \
} while (0)

/* By including cmath and not math.h, we get some prototypes in std:: ;
 * This happens to *also* be the case with math.h on some boxes, and that
 * may perhaps be a bug. Anyway, it seems that the cmath way is better.
 */

/* Initialize a polynomial of degree d */
void
double_poly_init (double_poly_ptr p, int d)
{
    if (d < 0) {
        p->deg = -1;
        p->alloc = 0;
        p->coeff = nullptr;
    } else {
        ASSERT_ALWAYS(d < INT_MAX);
        p->alloc = d + 1;
        p->deg = -1;
        p->coeff = (double*) malloc ((d + 1) * sizeof (double));
        FATAL_ERROR_CHECK(p->coeff == nullptr, "malloc failed");
    }
    STATIC_ANALYSIS_ASSERT_DATA_HEALTH(p);
}

/* Clear a polynomial */
void
double_poly_clear (double_poly_ptr p)
{
    free(p->coeff);
    memset(p, 0, sizeof(*p));
}

/* realloc to at least nc coefficients */
/* This never shrinks (as in mpz_poly) */
void
double_poly_realloc (double_poly_ptr p, unsigned int nc)
{
    ASSERT_ALWAYS(nc <= (unsigned int) INT_MAX);
    if (p->alloc < nc) {
        p->alloc = nc;
        checked_realloc(p->coeff, p->alloc);
    }
    STATIC_ANALYSIS_ASSERT_DATA_HEALTH(p);
}



/*
 * Set the degree of f.
 */
#if 0
void double_poly_set_degree(double_poly_ptr f, int deg)
{
    ASSERT_ALWAYS(0);
    /* do we really need this one, like this ? */
    double_poly_realloc(f, deg + 1);
    f->deg = deg;
}
#endif

/* Set r = s. */
void
double_poly_set (double_poly_ptr r, double_poly_srcptr s)
{
    if (r == s) return;
    if (s->deg == -1) { r->deg = -1; return; }
    double_poly_realloc(r, s->deg + 1);
    memcpy(r->coeff, s->coeff, (s->deg+1) * sizeof(double));
    r->deg = s->deg;
}
void double_poly_set_xi(double_poly_ptr s, int i)
{
    ASSERT_ALWAYS(i >= -1);
    double_poly_realloc(s, i+1);
    s->deg = i;
    if (i<0) return;
    if (i) memset(s->coeff, 0, sizeof(double) * i);
    s->coeff[i]=1;
}

void
double_poly_set_zero (double_poly_ptr r)
{
    r->deg = -1;
}

/* Evaluate the polynomial p at point x */
double
double_poly_eval (double_poly_srcptr p, const double x)
{
    double r;
    unsigned int k;
    const double *f = p->coeff;
    int const deg = p->deg;

    switch (deg) {
        case -1: return 0;
        case 0: return f[0];
        case 1: return f[0]+x*f[1];
        case 2: return f[0]+x*(f[1]+x*f[2]);
        case 3: return f[0]+x*(f[1]+x*(f[2]+x*f[3]));
        case 4: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*f[4])));
        case 5: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*f[5]))));
        case 6: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*f[6])))));
        case 7: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*f[7]))))));
        case 8: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*f[8])))))));
        case 9: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*(f[8]+x*f[9]))))))));
        default: for (r = f[deg], k = deg - 1; k != UINT_MAX; r = r * x + f[k--]); return r;
    }
}

/* divide f by x-r, and return the quotient (which is f(r)) */
double double_poly_div_linear(double_poly_ptr q, double_poly_srcptr p, const double r)
{
    /* supports self-assignment */
    const double *f = p->coeff;
    int const deg = p->deg;
    double u = f[deg];
    for (int k = p->deg - 1; k >= 0; k--) {
        double const c = f[k];
        q->coeff[k] = u;
        u = u * r + c;
    }
    q->deg = p->deg - 1;
    return u;
}


/* Evaluate the polynomial p at x:y */
double
double_poly_eval_homogeneous (double_poly_srcptr p, double x, double y)
{
    const double * f = p->coeff;
    int const deg = p->deg;

    switch (deg) {
        case -1: return 0;
        case 0: return f[0];
        case 1: return y*f[0]+x*f[1];
        case 2: return y*y*f[0]+x*(y*f[1]+x*f[2]);
    }

    double s = 0;
    double px = 1;
    for(int k = 0 ; k <= deg ; k++) {
        s = s * y + f[k] * px;
        px = px * x;
    }
    return s;
}

/* return e and set r such that r*2^e = x */
static long
mpz_set_d_exp (mpz_t r, double x)
{
    long e;

    if (x == 0.0)
    {
        mpz_set_ui (r, 0);
        return 0;
    }

    e = ilogb (x) - 52;
    mpz_set_d (r, ldexp (x, -e));
    return e;
}

/* Same as double_poly_eval, but uses arbitrary precision to avoid
   cancellations */
double
double_poly_eval_safe (double_poly_srcptr p, double x)
{
  mpz_t xm, vm, fm;
  long xe, ve, fe;
  const double *f = p->coeff;
  const unsigned int d = p->deg;
  unsigned int k;
  double r;

  mpz_init (xm);
  mpz_init (vm);
  mpz_init (fm);
  xe = mpz_set_d_exp (xm, x);
  ve = mpz_set_d_exp (vm, f[d]);
  for (k = d - 1; k != UINT_MAX; k--)
    {
      /* multiply by x = xm*2^xe */
      mpz_mul (vm, vm, xm);
      ve += xe;
      /* add f[k] = fm*2^fe */
      fe = mpz_set_d_exp (fm, f[k]);
      if (fe < ve)
        {
          mpz_mul_2exp (vm, vm, ve - fe);
          ve = fe;
        }
      else
        mpz_mul_2exp (fm, fm, fe - ve);
      mpz_add (vm, vm, fm);
    }
  r = mpz_get_d (vm);
  mpz_clear (xm);
  mpz_clear (vm);
  mpz_clear (fm);
  return ldexp (r, ve);
}

/* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
   refines that root by dichotomy with n iterations.
   Assumes sa is of same sign as g(a).
*/
double
double_poly_dichotomy (double_poly_srcptr p, double a, double b, double sa,
                       unsigned int n MAYBE_UNUSED)
{
  double s;

  for(;;) {
#if defined(__i386)
      /* See comment in utils/usp.c on this. We want to avoid comparison
       * involving an extended precision double ! */
      { volatile double ms = (a + b) * 0.5; s = ms; }
#else
      s = (a + b) * 0.5;
#endif
      if (s == a || s == b) return s;
      if (double_poly_eval (p, s) * sa > 0)
	a = s;
      else
	b = s;
  }
  return s;
}

/* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
   refines that root by the weighted false position method
   Assumes sa is of same sign as g(a).
*/
static double
double_poly_falseposition (double_poly_srcptr p, double a, double b, double pa)
{
  double pb;
  int side=0;
  double a0=a, b0=b, pa0=pa;

  pb = double_poly_eval(p, b);

  for(;;) {
      double s, middle;
#if defined(__i386)
      /* See above */
      { volatile double ms = (a*pb-b*pa)/(pb-pa); s = ms; }
      { volatile double ms = (a + b) * 0.5; middle = ms; }
#else
      s = (a*pb-b*pa)/(pb-pa);
      middle = (a + b) * 0.5;
#endif
      /* It may happen that because of overflow, (a*pb-b*pa)/(pb-pa)
       * reaches s==a or s==b too early. If it so happens that we're
       * doing this, while the middle cut doesn't behave this way, use
       * the middle cut instead.
       *
       * Note that almost by design, this countermeasure also cancels
       * some of the benefit of the false position method.
       */
      if (s < a || s > b || ((s == a || s == b) && !(middle == a || middle == b)))
          s = middle;
      if (s == a || s == b) return s;
      double const ps = double_poly_eval (p, s);
      if (ps * pa > 0) {
          a = s; pa = ps;
          if (side==1) pb /= 2;
          side=1;
      } else {
          b = s; pb = ps;
          if (side==-1) pa /= 2;
          side=-1;
      }
      if (isnan(b)) {
          return double_poly_dichotomy(p, a0, b0, pa0, 0);
      }
  }
}

/* Stores the derivative of f in df. f and df may be identical.
   Assumes df has been initialized with degree at least f->deg-1. */
    void
double_poly_derivative(double_poly_ptr df, double_poly_srcptr f)
{
    double_poly_realloc(df, f->deg);      /* number of coefficients ! */
    /* supports self-assignment */
    for (int n = 1; n <= f->deg; n++)
        df->coeff[n - 1] = f->coeff[n] * (double) n;
    df->deg = f->deg - 1;
}

/* Stores the product of f and g in h (h = f * g).
 * supports self-assignment
 */
void
double_poly_mul(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
    if (h == f || h == g) {
        double_poly nh;
        double_poly_init(nh, f->deg + g->deg);
        double_poly_mul(nh, f, g);
        double_poly_swap(nh, h);
        double_poly_clear(nh);
        return;
    }
    double_poly_realloc(h, f->deg + g->deg + 1);
    h->deg = f->deg + g->deg;
    memset (h->coeff, 0, sizeof(double) * (h->deg + 1));
    for (int i = 0 ; i <= f->deg ; i++) {
        for(int j = 0 ; j <= g->deg ; j++) {
            h->coeff[i+j] += f->coeff[i] * g->coeff[j];
        }
    }
}

/* Stores the sum of f and g in h (h = f + g).
 * supports self-assignment (realloc does not invalidate)
 */
/* cleans up leading coefficients too */
void
double_poly_add(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
    double_poly_realloc(h, MAX(f->deg, g->deg) + 1);
    int i = 0;
    for (; i <= f->deg && i <= g->deg; ++i)
        h->coeff[i] = f->coeff[i] + g->coeff[i];
    for (; i <= f->deg ; ++i)
        h->coeff[i] = f->coeff[i];
    for (; i <= g->deg ; ++i)
        h->coeff[i] = g->coeff[i];
    double_poly_cleandeg(h, MAX(f->deg, g->deg));
}

/* subtract */
/* cleans up leading coefficients too */
void
double_poly_sub(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
    double_poly_realloc(h, MAX(f->deg, g->deg) + 1);
    int i = 0;
    for (; i <= f->deg && i <= g->deg; ++i)
        h->coeff[i] = f->coeff[i] - g->coeff[i];
    for (; i <= f->deg ; ++i)
        h->coeff[i] = f->coeff[i];
    for (; i <= g->deg ; ++i)
        h->coeff[i] = -g->coeff[i];
    double_poly_cleandeg(h, MAX(f->deg, g->deg));
}

void double_poly_addmul(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
    cxx_double_poly fg;
    double_poly_mul(fg, f, g);
    double_poly_add(h, h, fg);
}

void double_poly_submul(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
    cxx_double_poly fg;
    double_poly_mul(fg, f, g);
    double_poly_sub(h, h, fg);
}
void double_poly_addmul_double(double_poly_ptr h, double_poly_srcptr g, double m)
{
    double_poly_realloc(h, MAX(h->deg, g->deg) + 1);
    int i = 0;
    for (; i <= h->deg && i <= g->deg; ++i)
        h->coeff[i] += m * g->coeff[i];
    for (; i <= g->deg ; ++i)
        h->coeff[i] = m * g->coeff[i];
    double_poly_cleandeg(h, MAX(h->deg, g->deg));
}

void double_poly_submul_double(double_poly_ptr h, double_poly_srcptr g, double m)
{
    double_poly_realloc(h, MAX(h->deg, g->deg) + 1);
    int i = 0;
    for (; i <= h->deg && i <= g->deg; ++i)
        h->coeff[i] -= m * g->coeff[i];
    for (; i <= g->deg ; ++i)
        h->coeff[i] = -m * g->coeff[i];
    double_poly_cleandeg(h, MAX(h->deg, g->deg));
}


void
double_poly_neg(double_poly_ptr h, double_poly_srcptr f)
{
    double_poly_realloc(h, f->deg + 1);
    for (int i = 0; i <= f->deg ; ++i)
        h->coeff[i] = -f->coeff[i];
    h->deg = f->deg;
}

/* Revert the coefficients: g(x) => f(1/x) * x^degree */
/* cleans up leading coefficients too */
void
double_poly_revert (double_poly_ptr g, double_poly_srcptr f)
{
    if (g == f) {
        const int d = f->deg;

        if (d <= 0)
            return;

        /* if d is even, nothing to do for k=d/2 */
        for (int k = 0; k <= (d - 1) / 2; k++)
        {
            double const tmp = f->coeff[k];
            f->coeff[k] = f->coeff[d - k];
            f->coeff[d - k] = tmp;
        }
    } else {
        double_poly_realloc(g, f->deg + 1);
        for(int i = 0 ; i <= f->deg ; i++) {
            g->coeff[f->deg-i] = f->coeff[i];
        }
        double_poly_cleandeg(g, f->deg);
    }
}

/* Change sign of variable: x -> -x */
static void
double_poly_neg_x (double_poly_ptr r, double_poly_srcptr s)
{
    int i;
    double_poly_realloc(r, s->deg + 1);
    r->deg = s->deg;
    for (i = 0; i < s->deg; i += 2) {
        r->coeff[i]     = s->coeff[i]; /* Even power coeff */
        r->coeff[i + 1] = -s->coeff[i + 1]; /* Odd power */
    }
    if (i == s->deg) /* iff deg is even */
        r->coeff[i] = s->coeff[i];
}

static unsigned int
recurse_roots(double_poly_srcptr poly, double *roots,
              const unsigned int sign_changes, const double s)
{
  unsigned int new_sign_changes = 0;
  if (poly->deg <= 0) {
      /* A constant polynomial (degree 0 or -\infty) has no sign changes */
  } else if (poly->deg == 1) {
      /* A polynomial of degree 1 can have at most one sign change in (0, s),
         this happens iff poly(0) = poly[0] and poly(s) have different signs */
      if (poly->coeff[0] * double_poly_eval(poly, s) < 0) {
          new_sign_changes = 1;
          roots[0] = - poly->coeff[0] / poly->coeff[1];
      }
  } else {
      /* invariant: sign_changes is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_changes-1], and roots[sign_changes] = s. */
      double a = 0.0;
      double va = poly->coeff[0]; /* value of poly at x=0 */
      for (unsigned int l = 0; l <= sign_changes; l++)
        {
          /* b is a root of dg[k+1], or s, the end of the interval */
          const double b = (l < sign_changes) ? roots[l] : s;
          const double vb = double_poly_eval (poly, b);
          if (va * vb < 0) /* root in interval [va, vb] */
            roots[new_sign_changes++] = double_poly_falseposition (poly, a, b, va);
          a = b;
          va = vb;
        }
  }

  return new_sign_changes;
}

/* Return a bound on the positive roots of p.
   Assume the leading coefficient of p is positive, then for a positive root
   r we have
   p[d]*r^d + ... + p[1]*r + p[0] = 0 thus
   p[d]*r^d <= -p[d-1]*r^(d-1) - ... - p[1]*r - p[0]
            <= max(-p[d-1],0)*r^(d-1) + ... + max(-p[1],0)*r + max(-p[0],0)
   thus q(r) <= 0 where q is the degree-d polynomial formed from p as follows:
   * q[d] = p[d]
   * q[i] = p[i] if p[i] < 0, and 0 otherwise for 0 <= i < d.
   Since q has a unique positive root, say r0, and q(r) < 0 iff r < r0,
   then positive roots of p are bounded by r0.
*/
double
double_poly_bound_roots (double_poly_srcptr p)
{
    ASSERT_ALWAYS(p->deg >= 0);
  int const d = p->deg;
  double_poly q;
  double s;

  s = p->coeff[d] > 0 ? 1.0 : -1.0;
  double_poly_init (q, d);
  for (int i = 0; i < d; i++)
    q->coeff[i] = (s * p->coeff[i] < 0) ? s * p->coeff[i] : 0.0;
  q->coeff[d] = fabs (p->coeff[d]);
  double_poly_cleandeg(q, d);
  s = 1.0;
  while (double_poly_eval (q, s) < 0)
    s = s + s;
  double_poly_clear (q);
  return s;
}

/* put in roots[0], roots[1], ..., roots[k-1] all k roots of poly in (0, s],
   and return the number k of roots */
unsigned int
double_poly_compute_roots(double *roots, double_poly_srcptr poly, double s)
{
    const int d = poly->deg;
    double_poly *dg; /* derivatives of poly */

    /* The roots of the zero polynomial are ill-defined. Bomb out */
    ASSERT_ALWAYS(d>=0);

    /* Handle constant polynomials separately */
    if (d == 0)
        return 0; /* Constant non-zero poly -> no roots */

    dg = (double_poly *) malloc (d * sizeof (double_poly));
    FATAL_ERROR_CHECK(dg == nullptr, "malloc failed");

    dg[0]->deg = poly->deg;
    dg[0]->coeff = poly->coeff;

    for (int k = 1; k < d; k++) {
        /* dg[k] is the k-th derivative, thus has degree d-k, i.e., d-k+1
           coefficients */
        double_poly_init (dg[k], d - k);
        double_poly_derivative (dg[k], dg[k - 1]);
    }

    unsigned int sign_changes = 0;
    for (int k = d; k > 0; k--)
        sign_changes = recurse_roots(dg[k - 1], roots, sign_changes, s);

    for (int k = 1; k < d; k++)
        double_poly_clear (dg[k]);
    free (dg);

    return sign_changes;
}

/* compute all roots whose absolute value is <= B */
unsigned int
double_poly_compute_all_roots_with_bound (double *roots,
                                          double_poly_srcptr poly,
                                          double B)
{
  /* Positive roots */
  double bound = double_poly_bound_roots (poly);
  if (B < bound)
    bound = B;
  unsigned int const nr_roots_pos = double_poly_compute_roots (roots, poly, bound);
  /* Negative roots */
  double_poly t; /* Copy of poly which gets sign-flipped */
  double_poly_init (t, poly->deg);
  double_poly_neg_x (t, poly);
  bound = double_poly_bound_roots (t);
  if (B < bound)
    bound = B;
  unsigned int nr_roots_neg =
    double_poly_compute_roots (roots + nr_roots_pos, t, bound);
  double_poly_clear(t);
  /* Flip sign of negative roots */
  for (unsigned int i = 0; i < nr_roots_neg; i++)
    roots[nr_roots_pos + i] *= -1.;

  /* check if zero is a root */
  if (poly->coeff[0] == 0.0)
    roots[nr_roots_pos + nr_roots_neg++] = 0.0;

  return nr_roots_pos + nr_roots_neg;
}

unsigned int
double_poly_compute_all_roots (double *roots, double_poly_srcptr poly)
{
  return double_poly_compute_all_roots_with_bound (roots, poly, DBL_MAX);
}

std::string cxx_double_poly::print_poly(std::string const& var) const
{
    std::ostringstream os;
    for(int i = 0 ; i <= x->deg ; i++) {
        if (x->coeff[i] == 0) continue;
        if (x->coeff[i] > 0 && os.str().size())
            os << "+";
        if (i == 0) {
            os << x->coeff[i];
        } else {
            if (x->coeff[i] == -1) {
                os << "-";
            } else if (x->coeff[i] != 1) {
                os << x->coeff[i] << "*";
            }
            os << var;
            if (i > 1) os << "^" << i;
        }

    }
    return os.str();
}


/* Print polynomial with floating point coefficients. Assumes f[deg] != 0
   if deg > 0. */
void 
double_poly_print (FILE *stream, double_poly_srcptr p, char *name)
{
    cxx_double_poly F;
    double_poly_set(F, p);
    std::string const s = F.print_poly(name);
    fputs(s.c_str(), stream);
}

int
double_poly_asprint (char **t, double_poly_srcptr p, char *name)
{
    cxx_double_poly F;
    double_poly_set(F, p);
    std::string const s = F.print_poly(name);
    int const n = s.size();
    *t = (char*) malloc(n + 1);
    memcpy(*t, s.c_str(), n + 1);
    return n;
}

void
double_poly_set_mpz_poly (double_poly_ptr p, mpz_poly_srcptr q)
{
    double_poly_realloc(p, q->deg + 1);
    for (int i = 0; i <= q->deg; i++)
        p->coeff[i] = mpz_get_d (mpz_poly_coeff_const(q, i));
    p->deg = q->deg;
}

/*
 * set the degree to deg, or maybe less given the possibly zero leading
 * coeffs.
 *
 * This will not zero out coefficients which have been written to beyond
 * the stored degree for the polynomial. Of course, those need to be in
 * the allocated range.
 *
 * When cleandeg is called with a degree more than the allocated range,
 * the polynomial is reallocated to the requested degree, but the degree
 * is set only to the largest possible degree given the previous
 * allocated amount.
 */
void double_poly_cleandeg(double_poly_ptr f, int deg)
{
    ASSERT_ALWAYS(deg >= -1);
    if (f->alloc < (unsigned int) (deg + 1)) {
        double_poly_realloc(f, deg + 1);
        memset(f->coeff + f->deg + 1, 0, (deg - f->deg) * sizeof(double));
        /* we're not increasing f->deg, since we know that the rest is
         * made of zeroes */
        deg = f->alloc-1;
    }
    for( ; deg >= 0 && f->coeff[deg] == 0 ; deg--);
    f->deg = deg;
}

/*
 * Compare two polynomials. Returns -1, 0, or 1 depending on whether a<b,
 * a==b, or a>b.
 *
 * Assume that the polynomials are both normalized.
 */
int double_poly_cmp(double_poly_srcptr a, double_poly_srcptr b)
{
    ASSERT_ALWAYS(0);
#ifndef NDEBUG
    if (a->deg > -1) {
        ASSERT(a->coeff[a->deg] != 0);
    }
    if (b->deg > -1) {
        ASSERT(b->coeff[b->deg] != 0);
    }
#endif // NDEBUG

    int r = (a->deg > b->deg) - (b->deg > a->deg);
    for(int d = a->deg; !r && d >= 0 ; d--) {
        double const s = a->coeff[d] - b->coeff[d];
        r = (s >= 0) - (s <= 0);
    }
    return r;
}

/*
 * f = mul * g.
 */
void double_poly_mul_double(double_poly_ptr f, double_poly_srcptr g,
    double mul)
{
    if (g->deg < 0) {
        f->deg = -1;
        return;
    }

    ASSERT(mul != 0.0);
    ASSERT(g->deg >= 0);
    ASSERT(f->deg >= g->deg);

    double_poly_set(f, g);

    for (int i = 0; i <= f->deg; i++) {
        f->coeff[i] = f->coeff[i] * mul;
    }
    double_poly_cleandeg(f, f->deg);
}

#if 0
/*
 * Return the content of f.
 */
static double double_poly_content(double_poly_srcptr f)
{
    ASSERT_ALWAYS(0);
    ASSERT(f->deg > -1);

    int64_t gcd = (int64_t)f->coeff[0];
    for (int i = 1; i <= f->deg; i++) {
        gcd = gcd_int64(gcd, (int64_t)f->coeff[i]);
    }
    return fabs((double) gcd);
}
#endif

/*
 * Return the leading coefficient of f, even with f not normalized.
 * Returns 0 for the null polynomial.
 */
static double double_poly_lc(double_poly_srcptr f)
{
    for(int deg = f->deg ; deg >= 0 ; deg--) {
        if (f->coeff[deg] != 0)
            return f->coeff[deg];
    }
    return 0;
}

/*
 * Return 0 if the polynomial does not contain inf or nan, 1 otherwise.
 */
static int double_poly_isnan_or_isinf(double_poly_srcptr f)
{
    for (int i = 0; i <= f->deg; i++) {
        if (isnan(f->coeff[i]) || isinf(f->coeff[i])) {
            return 1;
        }
    }
    return 0;
}

/*
 * Compute the pseudo division of a and b such that
 *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that deg(a) >= deg(b) and b is not the zero polynomial.
 *
 * Return 0 on failure, 1 otherwise.
 *
 * No overlap allowed.
 */
static int double_poly_pseudo_division(double_poly_ptr q, double_poly_ptr r,
    double_poly_srcptr a, double_poly_srcptr b)
{
    ASSERT(a->deg >= b->deg);
    ASSERT(b->deg != -1);

    if (q) double_poly_realloc(q, a->deg - b->deg + 1);

    int const m = a->deg;
    int const n = b->deg;
    double d = double_poly_lc(b);
    int e = m - n + 1;
    double_poly s;
    double_poly_init(s, m-n);

    if (q) double_poly_set_zero(q);

    double_poly_set(r, a);

    while (r->deg >= n) {
        double_poly_cleandeg(r, r->deg);
        double_poly_set_xi(s, r->deg - n);
        s->coeff[r->deg - n] = double_poly_lc(r);

        if (q) {
            double_poly_mul_double(q, q, d);
            double_poly_add(q, s, q);
        }
        double_poly_mul_double(r, r, d);

        double_poly_mul(s, b, s);
        int const nrdeg = r->deg - 1;
        double_poly_sub(r, r, s);
        /* We'd like to enforce this because the subtraction may miss the
         * cancellation of the leading term due to rounding.
         */
        double_poly_cleandeg(r, nrdeg);
        e--;
        //TODO: is it necessary?
        if (double_poly_isnan_or_isinf(r) || r->deg == -1) {
            double_poly_clear(s);
            return 0;
        }
    }
    double_poly_clear(s);

    ASSERT(e >= 0);

    d = pow(d, (double) e);
    if (q) double_poly_mul_double(q, q, d);
    double_poly_mul_double(r, r, d);

    return 1;
}

static int double_poly_pseudo_remainder(double_poly_ptr r,
    double_poly_srcptr a, double_poly_srcptr b)
{
    return double_poly_pseudo_division(nullptr, r, a, b);
}

//TODO: follow the modification of mpz_poly_resultant.
double double_poly_resultant(double_poly_srcptr p, double_poly_srcptr q)
{
    if (p->deg == -1 || q->deg == -1) {
        return 0;
    }

    ASSERT(p->coeff[p->deg] != 0);

    ASSERT(q->coeff[q->deg] != 0);

    double_poly a;
    double_poly b;
    double_poly r;
    double_poly_init(a, p->deg);
    double_poly_init(b, q->deg);
    double_poly_init(r, MAX(p->deg, q->deg));
    double_poly_set(a, p);
    double_poly_set(b, q);

    int s = 1;
    int d;

    /* This does not really make sense on double_poly, does it ?
     *
    double g = double_poly_content(a);
    double h = double_poly_content(b);
    double_poly_mul_double(a, a, 1 / g);
    double_poly_mul_double(b, b, 1 / h);
    double t = pow(g, (double) b->deg) * pow(h, (double) a->deg);
    */

    int pseudo_div = 1;

    double g = 1;
    double h = 1;

    if (a->deg < b->deg) {
        double_poly_swap(a, b);

        if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
            s = -1;
        }
    }

    while (b->deg > 0) {
        //TODO: verify if it is necessary.
        d = a->deg - b->deg;;

        if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
            s = -s;
        }
        pseudo_div = double_poly_pseudo_remainder(r, a, b);
        if (!pseudo_div)
            break;
        /*double_poly_degree(r);*/
        double_poly_set(a, b);

        ASSERT(d >= 0);

        double_poly_mul_double(b, r, 1 / (g * pow(h, (double) d)));
        //TODO: a is normalized, so g = a->coeff[a->deg].
        g = double_poly_lc(a);

#ifdef NDEBUG
        if (d == 0) {
            ASSERT(h == 1);
        }
#endif // NDEBUG

        h = pow(h, (double) (d - 1));
        h = pow(g, (double) d) / h;
    }

    if (pseudo_div) {
        ASSERT(a->deg > 0);

        //For now, if b->deg == -1, pseudo_div == 0.
        if (b->deg == -1) {
            ASSERT(0);
        } else {

            ASSERT(a->deg >= 0);

            h = pow(b->coeff[0], (double) a->deg) / pow(h, (double) (a->deg - 1));

            // h = (double)s * t * h;
            // see further above: got rid of t.
            h *= s;
        }
    } else {
        // fall back. Why is it needed ? Is it covered by tests or not ?
        //
        //TODO: use last version of a and b in pseudo_division.
        mpz_t res;
        mpz_init(res);
        mpz_poly f, g;
        mpz_poly_init(f, p->deg);
        mpz_poly_init(g, q->deg);
        mpz_poly_set_double_poly(f, p);
        mpz_poly_set_double_poly(g, q);

        mpz_poly_resultant(res, f, g);
        h = mpz_get_d(res);

        mpz_poly_clear(f);
        mpz_poly_clear(g);
        mpz_clear(res);
    }
    double_poly_clear(a);
    double_poly_clear(b);
    double_poly_clear(r);

    return h;
}

/* swap f and g */
void
double_poly_swap (double_poly_ptr f, double_poly_ptr g)
{
    int const i = f->deg; f->deg = g->deg; g->deg = i;
    int const j = f->alloc; f->alloc = g->alloc; g->alloc = j;
    double * t = f->coeff; f->coeff = g->coeff; g->coeff = t;
}

/* Parse a space-separated string as polynomial coefficients.
 * Exactly one space between coefficients, none at start or end of string.
 * If poly is nullptr, nothing is written, but number of coefficients is
 * returned.
 */
void double_poly_set_string(double_poly_ptr poly, const char *str)
{
    int n = 0;
    /* parse a first time to find n */
    for(const char * next = str ; next ; n++) {
        int dx;
        double coeff;
        if (sscanf(next, "%lf%n", &coeff, &dx) != 1)
            break;
        next += dx;
        ASSERT_ALWAYS(!*next || isspace(*next));
        for( ; *next && isspace(*next) ; next++);
    }
    double_poly_realloc (poly, n); /* number of coefficients ! */
    n = 0;
    for(const char * next = str ; next ; n++) {
        int dx;
        if (sscanf(next, "%lf%n", &poly->coeff[n], &dx) != 1)
            break;
        next += dx;
        for( ; *next && isspace(*next) ; next++);
    }
    double_poly_cleandeg(poly, n-1);
}

