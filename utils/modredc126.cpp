#include "cado.h" // IWYU pragma: keep
#include <cstdint>     /* AIX wants it first (it's a bug) */
#include <cstdlib>        // for abort
#include "modredc126.hpp"

/* Divide residue by 3. Returns 1 if division is possible, 0 otherwise.
   Assumes that a+3m does not overflow */

typedef ModulusREDC126 Modulus;
#include "mod_common.cpp"
#include "macros.h"

bool
ModulusREDC126::div3 (Residue &r, const Residue &a) const
{
    Residue t(*this);
    uint64_t a3 = (a.r[1] % 256 + a.r[1] / 256 +
                   a.r[0] % 256 + a.r[0] / 256) % 3;
    const uint64_t m3 = (m[0] % 256 + m[0] / 256 +
                         m[1] % 256 + m[1] / 256) % 3;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue a_backup(*this);
    set (a_backup, a);
#endif

    if (m3 == 0)
        return 0;

    set(t, a);

    /* Make t[1]:t[0] divisible by 3 */
    if (a3 != 0)
    {
        if (a3 + m3 == 3) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
        {
            u64arith_add_2_2 (&(t.r[0]), &(t.r[1]), m[0], m[1]);
        }
        else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
        {
            u64arith_add_2_2 (&(t.r[0]), &(t.r[1]), m[0], m[1]);
            u64arith_add_2_2 (&(t.r[0]), &(t.r[1]), m[0], m[1]);
        }

        /* Now t[1]:t[0] is divisible by 3 */
        ASSERT_EXPENSIVE ((t.r[0] % 3 + t.r[1] % 3) % 3 == 0);
    }

    /* a = a1 * 2^w + a0, 3|a
       Let a = a' * 3 * 2^w + a'', a'' < 3 * 2^w.
       3 | a'', a'' / 3 < 2^w
       So a / 3 = a' * w + a'' / 3
       a' = trunc(a1 / 3)
       a'' = a0 * 3^{-1} (mod 2^w)
       Hence we get the correct result with one one-word multiplication
       and one one-word truncating division by a small constant.
    */

    r.r[1] = t.r[1] / 3;
    r.r[0] = t.r[0] * UINT64_C(0xaaaaaaaaaaaaaaab); /* 1/3 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    add (t, r, r);
    add (t, t, r);
    ASSERT_EXPENSIVE (equal (a_backup, t));
#endif

    return 1;
}


/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */

bool
ModulusREDC126::div5 (Residue &r, const Residue &a) const
{
  /* inv_5[i] = -1/i (mod 5) */
  const uint64_t inv_5[5] = {0,4,2,3,1};
  const uint64_t c = UINT64_C(0xcccccccccccccccd); /* 1/5 (mod 2^64) */
  
  return divn (r, a, 5, 1, inv_5, c);
}


/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */

bool
ModulusREDC126::div7 (Residue &r, const Residue &a) const
{
  /* inv_7[i] = -1/i (mod 7) */
  const uint64_t inv_7[7] = {0,6,3,2,5,4,1};
  const uint64_t c = UINT64_C(0x6db6db6db6db6db7); /* 1/7 (mod 2^64) */
  return divn (r, a, 7, 2, inv_7, c);
}


/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */

bool
ModulusREDC126::div11 (Residue &r, const Residue &a) const
{
  /* inv_11[i] = -1/i (mod 11) */
  const uint64_t inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1}; 
  const uint64_t c = UINT64_C(0x2e8ba2e8ba2e8ba3); /* 1/11 (mod 2^64) */
  return divn (r, a, 11, 5, inv_11, c);
}


/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */

bool
ModulusREDC126::div13 (Residue &r, const Residue &a) const
{
  /* inv_13[i] = -1/i (mod 13) */
  const uint64_t inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1}; 
  const uint64_t c = UINT64_C(0x4ec4ec4ec4ec4ec5); /* 1/13 (mod 2^64) */
  return divn (r, a, 13, 3, inv_13, c);
}


void
ModulusREDC126::gcd (Integer &r, const Residue &A) const
{
    uint64_t a[2], b[2];
    int sh;

    /* Since we do REDC arithmetic, we must have m odd */
    ASSERT_EXPENSIVE (m[0] % 2 != 0);

    if (is0(A))
    {
        getmod(r);
        return;
    }

    a[0] = A.r[0];
    a[1] = A.r[1];
    b[0] = m[0];
    b[1] = m[1];

    while (a[1] != 0 || a[0] != 0)
    {
        /* Make a odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros [(unsigned char) a[0]];
            u64arith_shrd (&(a[0]), a[1], a[0], sh);
            *(int64_t *) &(a[1]) >>= sh;
        } while (sh == 8);
#else
        if (a[0] == 0) /* ctzl does not like zero input */
        {
            a[0] = a[1];
            a[1] = ((int64_t)a[1] < 0L) ? (uint64_t) (-1L) : 0;
        }
        sh = u64arith_ctz (a[0]);
        u64arith_shrd (&(a[0]), a[1], a[0], sh);
        *(int64_t *) &(a[1]) >>= sh;
#endif

        /* Try to make the low two bits of b[0] zero */
        ASSERT_EXPENSIVE (a[0] % 2 == 1);
        ASSERT_EXPENSIVE (b[0] % 2 == 1);
        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2 (&(b[0]), &(b[1]), a[0], a[1]);
        else
            u64arith_sub_2_2 (&(b[0]), &(b[1]), a[0], a[1]);

        if (b[0] == 0 && b[1] == 0)
        {
            if ((int64_t) a[1] < 0)
            {
                a[1] = -a[1];
                if (a[0] != 0)
                    a[1]--;
                a[0] = -a[0];
            }
            r = Integer(a[0], a[1]);
            return;
        }

        /* Make b odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros [(unsigned char) b[0]];
            u64arith_shrd (&(b[0]), b[1], b[0], sh);
            *(int64_t *) &(b[1]) >>= sh;
        } while (sh == 8);
#else
        if (b[0] == 0) /* ctzl does not like zero input */
        {
            b[0] = b[1];
            b[1] = ((int64_t)b[1] < 0) ? (uint64_t) (-1) : 0;
        }
        sh = u64arith_ctz (b[0]);
        u64arith_shrd (&(b[0]), b[1], b[0], sh);
        *(int64_t *) &(b[1]) >>= sh;
#endif
        ASSERT_EXPENSIVE (a[0] % 2 == 1);
        ASSERT_EXPENSIVE (b[0] % 2 == 1);

        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2 (&(a[0]), &(a[1]), b[0], b[1]);
        else
            u64arith_sub_2_2 (&(a[0]), &(a[1]), b[0], b[1]);
    }

    if ((int64_t) b[1] < 0)
    {
        b[1] = -b[1];
        if (b[0] != 0)
            b[1]--;
        b[0] = -b[0];
    }
    r = Integer(b[0], b[1]);

    return;
}


/* Simple addition chains for small multipliers, used in the powering 
   functions with small bases. r and a may overlap, t must not overlap
   with anything. */
template <int M>
static inline void
simple_mul (ModulusREDC126::Residue &r, const ModulusREDC126::Residue &a,
    ModulusREDC126::Residue &t, const ModulusREDC126 &m);
template <>
void
simple_mul<2> (ModulusREDC126::Residue &r, const ModulusREDC126::Residue &a,
    ModulusREDC126::Residue &t MAYBE_UNUSED, const ModulusREDC126 &m)
{
    m.add (r, a, a); /* r = 2*a */
}
template <>
void
simple_mul<3> (ModulusREDC126::Residue &r, const ModulusREDC126::Residue &a,
    ModulusREDC126::Residue &t, const ModulusREDC126 &m)
{
    m.add (t, a, a); /* t = 2*a */
    m.add (r, t, a); /* r = 3*a */
}
template <>
void
simple_mul<5> (ModulusREDC126::Residue &r, const ModulusREDC126::Residue &a,
    ModulusREDC126::Residue &t, const ModulusREDC126 &m)
{
    m.add (t, a, a); /* t = 2*a */
    m.add (t, t, t); /* t = 4*a */
    m.add (r, r, a); /* r = 5*a */
}
template <>
void
simple_mul<7> (ModulusREDC126::Residue &r, const ModulusREDC126::Residue &a,
    ModulusREDC126::Residue &t, const ModulusREDC126 &m)
{
    m.add (t, a, a); /* r = 2*a */
    m.add (t, t, t); /* r = 4*a */
    m.add (t, t, t); /* r = 8*a */
    m.sub (r, t, a); /* r = 7*a */
}

template <int B, typename WordType>
static inline void npow_oneWord(
    WordType mask, const WordType word, typename ModulusREDC126::Residue &t,
    typename ModulusREDC126::Residue &u, const ModulusREDC126 &m)
{
    while (mask > 0) {
        m.sqr (t, t);
        if (word & mask) {
            simple_mul<B> (t, t, u, m);
        }
        mask >>= 1;
    }
}

/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are 
   implemented. Here, e is an uint64_t */
template <int B>
static inline void
npow (ModulusREDC126::Residue &r, const uint64_t e, const ModulusREDC126 &m)
{
    uint64_t mask;
    ModulusREDC126::Residue t(m), u(m);

    if (e == 0) {
        m.set1 (r);
        return;
    }

    m.set1 (t);
    simple_mul<B> (t, t, u, m); /* t = b */

    mask = (UINT64_C(1) << 63) >> u64arith_clz (e);
    ASSERT (e & mask);
    mask >>= 1;

    npow_oneWord<B>(mask, e, t, u, m);
    m.set (r, t);
}


/* Compute r = 2^e mod m. Here, e is an uint64_t */
void ModulusREDC126::pow2 (Residue &r, const uint64_t e) const
{
  npow<2> (r, e, *this);
}


/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are 
   implemented.  Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i. e_nrwords must be
   minimal, i.e., either e_nrwords == 0 or e[e_nrwords - 1] != 0. */
template <int B>
static inline void
npow (ModulusREDC126::Residue &r, const uint64_t *e,
      const size_t e_nrwords, const ModulusREDC126 &m)
{
    ModulusREDC126::Residue t(m), u(m);
    size_t i = e_nrwords;
    uint64_t mask;

    ASSERT(i == 0 || e[i - 1] != 0);
    
    if (i == 0) {
        m.set1 (r);
        return;
    }

    m.set1 (t);
    simple_mul<B> (t, t, u, m); /* t = b */

    mask = (UINT64_C(1) << 63) >> u64arith_clz (e[i - 1]);
    mask >>= 1;

    for ( ; i > 0; i--)
    {
        npow_oneWord<B>(mask, e[i - 1], t, u, m);
        mask = UINT64_C(1) << 63;
    }
    m.set (r, t);
}

template <int B>
static inline void
npow (ModulusREDC126::Residue &r,
      const ModulusREDC126::Integer &e, const ModulusREDC126 &m) {
    if (e.size() == 2) {
        uint64_t t[2];
        e.get(t, 2);
        npow<B> (r, t, e.size(), m); /* r = b^e mod m */
    } else if (e.size() <= 1) {
        npow<B> (r, e.getWord(0), m);
    } else {
      abort();
    }
}

/* Compute r = 2^e mod m.  Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
void
ModulusREDC126::pow2 (Residue &r, const uint64_t *e, const size_t e_nrwords) const
{
  npow<2> (r, e, e_nrwords, *this);
}

void
ModulusREDC126::pow2 (Residue &r, const Integer &e) const
{
  npow<2> (r, e, *this);
}

/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
bool
ModulusREDC126::sprp (const Residue &b) const
{
  Residue r(*this), minusone(*this);
  Integer mm1;
  int po2 = 0;
  bool i = 0;

  getmod (mm1);

  if (mm1 == 1)
    return 0;

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1 -= 1;
  po2 = mm1.ctz();
  mm1 >>= po2;

  set1 (minusone);
  neg (minusone, minusone);

  /* Exponentiate */
  pow (r, b, mm1);
  i = find_minus1 (r, minusone, po2);

  return i;
}


/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise. 
   We assume m is odd. */
bool
ModulusREDC126::sprp2 () const
{
    Residue r(*this), minusone(*this);
    int i = 0, po2 = 0;
    Integer mm1;

    getmod (mm1);
    if (mm1 == 1)
        return false;

    /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
    --mm1;
    po2 = mm1.ctz();
    mm1 >>= po2;

    set1 (minusone);
    neg (minusone, minusone);

    /* Exponentiate */
    pow2 (r, mm1);
    i = find_minus1 (r, minusone, po2);

    return i;
}


bool
ModulusREDC126::isprime () const
{
    Residue b(*this), minusone(*this), r1(*this);
    Integer n, mm1;
    bool r = false;
    int po2 = 0, i;

    getmod (n);

    if (n == 1)
        return false;

    if (n.getWord(0) % 2 == 0)
        return n == 2;

    /* Set mm1 to the odd part of m-1 */
    mm1 = n - 1;
    po2 = mm1.ctz();
    mm1 >>= po2;

    set1 (minusone);
    neg (minusone, minusone);

    /* Do base 2 SPRP test */
    pow2 (r1, mm1);
    /* If n is prime and 1 or 7 (mod 8), then 2 is a square (mod n)
       and one fewer squarings must suffice. This does not strengthen the
       test but saves one squaring for composite input */
    if (n.getWord(0) % 8 == 7) {
        if (!is1 (r1))
            goto end;
    } else if (!find_minus1 (r1, minusone, po2 - ((n.getWord(0) % 8 == 7) ? 1 : 0)))
        goto end; /* Not prime */

    /* Base 3 is poor at identifying composites == 1 (mod 3), but good at
       identifying composites == 2 (mod 3). Thus we use it only for 2 (mod 3) */
    i = n.getWord(0) % 3 + n.getWord(1) % 3;
    if (i == 1 || i == 4)
    {
        npow<7> (r1, mm1, *this); /* r = 7^mm1 mod m */
        if (!find_minus1 (r1, minusone, po2))
            goto end; /* Not prime */

        set_reduced (b, 61); /* Use addition chain? */
        pow (r1, b, mm1); /* r = 61^mm1 mod m */
        if (!find_minus1 (r1, minusone, po2))
            goto end; /* Not prime */

        npow<5> (r1, mm1, *this); /* r = 5^mm1 mod m */
        if (!find_minus1 (r1, minusone, po2))
            goto end; /* Not prime */

        /* These are the base 2,5,7,61 SPSP < 10^13 and n == 1 (mod 3) */

        r = n != UINT64_C(30926647201) &&
            n != UINT64_C(45821738881) &&
            n != UINT64_C(74359744201) &&
            n != UINT64_C(90528271681) &&
            n != UINT64_C(110330267041) &&
            n != UINT64_C(373303331521) &&
            n != UINT64_C(440478111067) &&
            n != UINT64_C(1436309367751) &&
            n != UINT64_C(1437328758421) &&
            n != UINT64_C(1858903385041) &&
            n != UINT64_C(4897239482521) &&
            n != UINT64_C(5026103290981) &&
            n != UINT64_C(5219055617887) &&
            n != UINT64_C(5660137043641) &&
            n != UINT64_C(6385803726241);
    }
    else
    {
        /* Case n % 3 == 0, 2 */

        npow<3> (r1, mm1, *this); /* r = 3^mm1 mod m */
        if (!find_minus1 (r1, minusone, po2))
            goto end; /* Not prime */

        npow<5> (r1, mm1, *this); /* r = 5^mm1 mod m */
        if (!find_minus1 (r1, minusone, po2))
            goto end; /* Not prime */

        /* These are the base 2,3,5 SPSP < 10^13 and n == 2 (mod 3) */

        r = n != UINT64_C(244970876021) &&
            n != UINT64_C(405439595861) &&
            n != UINT64_C(1566655993781) &&
            n != UINT64_C(3857382025841) &&
            n != UINT64_C(4074652846961) &&
            n != UINT64_C(5783688565841);
    }

end:
#if defined(PARI)
    printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
    return r;
}



int
ModulusREDC126::jacobi (const Residue &a_par) const
{
  Integer a, m, s;
  int t = 1;

  get (a, a_par);
  getmod (m);
  
  while (a != 0) {
    while ((a & 1) == 0) { /* TODO speedup */
      a >>= 1;
      if ((m & 7) == 3 || (m & 7) == 5)
        t = -t;
    }
    s = a; /* swap a and m */
    a = m;
    m = s;
    if ((a & 3) == 3 && (m & 3) == 3)
      t = -t;
    
    /* m is odd here */
    if (a >= m)
      a %= m;
  }
  if (m != 1)
    t = 0;
  
  return t;
}


bool
ModulusREDC126::inv (Residue &r, const Residue &A) const
{
    Integer a, b, u, v;
    int t, lsh;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue tmp(*this);

    set (tmp, A);
#endif

    assertValid(A);
    ASSERT_EXPENSIVE (m[0] % 2 != 0);

    if (is0(A))
        return 0;

    getmod (b);

    /* Let A = x*2^{2w}, so we want the Montgomery representation of 1/x,
       which is 2^{2w}/x. We start by getting a = x */

    /* We simply set a = x/2^{2w} and t=0. The result before correction
       will be 2^(2w+t)/x so we have to divide by t, which may be >64,
       so we may have to do one or more full and a variable width REDC. */
    /* TODO: If b[1] > 1, we could skip one of the two REDC */
    {
        Residue x(*this);
        set(x, A);
        redc1(x, x);
        get (a, x);
    }
    /* Now a = x/2^w */
    t = -64;

    u = 1;
    v = 0; /* 0 is a valid pointer */

    /* make a odd */
    lsh = a.ctz();
    t += lsh;
    a >>= lsh;

    // Here a and b are odd, and a < b
    do {
        /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
        ASSERT_EXPENSIVE (a < b);
        ASSERT_EXPENSIVE ((a & 1) == 1);
        ASSERT_EXPENSIVE ((b & 1) == 1);
        ASSERT_EXPENSIVE ((u & 1) == 1);
        ASSERT_EXPENSIVE ((v & 1) == 0);

        do {
            b -= a;
            v += u;
            ASSERT_EXPENSIVE ((b & 1) == 0);

            lsh = b.ctz();
            t += lsh;
            b >>= lsh;
            u <<= lsh;
        } while (a < b); /* ~50% branch taken :( */

        /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
        ASSERT_EXPENSIVE ((a & 1) == 1);
        ASSERT_EXPENSIVE ((b & 1) == 1);
        ASSERT_EXPENSIVE ((u & 1) == 0);
        ASSERT_EXPENSIVE ((v & 1) == 1);

        if (a == b)
            break;
        ASSERT_EXPENSIVE (a > b);

        /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
        do {
            a -= b;
            u += v;

            ASSERT_EXPENSIVE ((a & 1) == 0);
            lsh = a.ctz();
            a >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (b < a); /* about 50% branch taken :( */
        /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
    } while (a != b);

    if (a != 1) /* Non-trivial GCD */
        return 0;

    ASSERT_ALWAYS (t >= 0);

    /* Here, the inverse of a is u/2^t mod m. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t
       with impunity. */
    for ( ; t >= 64 ; t -= 64)
        redc1 (u, u);

    if (t > 0) {
        uint64_t s[5], k;
        k = ((u.getWord(0) * invm) & ((UINT64_C(1) << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2 (&(s[0]), &(s[1]), k, m[0]);
        /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
        u64arith_add_2_2 (&(s[0]), &(s[1]), u.getWord(0), u.getWord(1));
        /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
        /* s[0] == 0 (mod 2^t) */
        ASSERT_EXPENSIVE ((s[0] & ((1UL << t) - 1)) == 0);
        s[2] = 0;
        u64arith_mul_1_1_2 (&(s[3]), &(s[4]), k, m[1]);
        u64arith_add_2_2 (&(s[1]), &(s[2]), s[3], s[4]);

        /* Now shift s[2]:s[1]:s[0] right by t */
        u64arith_shrd (&(s[0]), s[1], s[0], t);
        u64arith_shrd (&(s[1]), s[2], s[1], t);

        u = Integer(s[0], s[1]);
        t = 0;
    }
    ASSERT_FOR_STATIC_ANALYZER(t == 0); // consume t

    u.get(r.r, 2);
#ifdef WANT_ASSERT_EXPENSIVE
    mul (tmp, tmp, r);
    ASSERT_EXPENSIVE (is1 (tmp));
#endif

    return 1;
}

bool
ModulusREDC126::batchinv (Residue *r, const uint64_t *a, const size_t n,
                          const Integer *c) const
{
  Residue R(*this);

  if (n == 0)
    return 1;

  r[0] = a[0];

  /* beta' = 2^64, beta = 2^128 */
  for (size_t i = 1; i < n; i++) {
    mul_ul(r[i], r[i-1], a[i]);
    /* r[i] = beta'^{-i} \prod_{0 <= j <= i} a[j] */
  }

  /* Computes R = beta^2/r[n-1] */
  if (!inv(R, r[n - 1]))
    return 0;
  /* R = beta^2 beta'^{n-1} \prod_{0 <= j < n} a[j]^{-1} */

  if (c != NULL) {
    Residue t(*this);
    t = *c;
    mul(R, R, t);
  } else {
    redc1(R, R); /* Assume c=1 */
    redc1(R, R);
  }
  /* R = beta beta'^{n-1} c \prod_{0 <= j < n} a[j]^{-1} */

  redc1(R, R);
  /* R = beta beta'^{n-2} c \prod_{0 <= j < n} a[j]^{-1} */

  for (size_t i = n-1; i > 0; i--) {
    /* Invariant: R = beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1} */

    mul(r[i], R, r[i-1]);
    /* r[i] := R * r[i-1] / beta
            = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) * (1/beta'^{i-1} \prod_{0 <= j <= i-1} a[j]) / beta
            = c a[i]^{-1} */

    mul_ul(R, R, a[i]);
    /* R := R * a[i] / beta'
         = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) * a[i] / beta'
         = beta beta'^{i-2} c \prod_{0 <= j < i} a[j]^{-1},
       thus satisfying the invariant for i := i - 1 */
  }
  /* Here have R = beta * beta'^{-1} / a[0]. We need to convert the factor
     beta to a factor of beta', so that the beta' cancel. */
  redc1(R, R); /* R := beta * beta'^{-1} / a[0] / beta',
                  with beta = beta'^2, this is 1/a[0] */
  set(r[0], R);
  return 1;
}

#if 0

ModulusREDC126::batch_Q_to_Fp_context_t *
ModulusREDC126::batch_Q_to_Fp_init (const Integer &num, const Integer &den) const
{
  batch_Q_to_Fp_context_t *context;
  Integer ratio, remainder;

  context = (batch_Q_to_Fp_context_t *) malloc(sizeof(batch_Q_to_Fp_context_t));
  if (context == NULL)
    return NULL;

  mod_initmod_int(context->m, den);

  /* Compute ratio = floor(num / den), remainder = num % den. We assume that
    ratio fits into uint64_t, and abort if it does not. We need only the
    low word of remainder. */
  remainder = num % den;
  ratio = num - remainder;
  ratio = ratio.divexact(den);
  ASSERT_ALWAYS(ratio.size() == 1);
  ratio.get(&(context->ratio_ul), 1);
  // ASSERT_ALWAYS(remainder.size() == 1);
  remainder.get(&(context->rem_ul), 1);
  if (remainder != 0)
    context->c = den - remainder; /* c = -remainder (mod den) */

  context->den_inv = u64arith_invmod(den.get()[0]);

  return context;
}


void
modredc2ul2_batch_Q_to_Fp_clear (modredc2ul2_batch_Q_to_Fp_context_t * context)
{
  mod_clearmod(context->m);
  mod_intclear(context->c);
  free(context);
}


int
modredc2ul2_batch_Q_to_Fp (uint64_t *r,
                           const modredc2ul2_batch_Q_to_Fp_context_t *context,
                           const uint64_t k, const int neg,
                           const uint64_t *p, const size_t n)
{
  Residue *tr;
  int rc = 1;

  tr = (Residue *) malloc(n * sizeof(Residue));
  for (size_t i = 0; i < n; i++) {
    mod_init_noset0(tr[i], context->m);
  }

  if (!modredc2ul2_batchinv_ul(tr, p, n, context->c, context->m)) {
    rc = 0;
    goto clear_and_exit;
  }

  for (size_t i = 0; i < n; i++) {
    uint64_t t;
    t = ularith_post_process_inverse(mod_intget_ul(tr[i]), p[i],
                                     context->rem_ul, context->den_inv,
                                     context->ratio_ul, k);
    if (neg && t != 0)
      t = p[i] - t;
    r[i] = t;
  }

clear_and_exit:
  for (size_t i = 0; i < n; i++) {
    mod_clear(tr[i], context->m);
  }
  free(tr);
  return rc;
}
#endif
