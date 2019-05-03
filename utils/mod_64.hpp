/* Some functions for modular arithmetic with residues and modulus with up
   to 64 bits. Due to inlining, this file must be included in the caller's
   source code with #include */

/* Naming convention: all function start with mod64, for 
   MODulus 64-bits, followed by underscore, functionality of function
  (add, mul, etc), and possibly underscore and specification of what argument
  types the function takes (_u64, etc). */

#ifndef MOD_64_H
#define MOD_64_H

/**********************************************************************/
#include <cassert>
#include <climits>
#include <cstdint>
#include "macros.h"
#include "u64arith.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

/*********************************************************************/
/* Helper macros, see also u64arith.h */

/* ================= Functions that are part of the API ================= */

/* Functions for the modulus */


class Modulus64 {
    /* Integers of the same width as the modulus */
    class Integer {
        unsigned long v[1];
    public:
        Integer() : v{0} {}
        Integer(const uint64_t a) : v{a} {}
        ~Integer(){}
        
        /* The two mod*_u64s() functions import/export Integer from/to an array of 
            uint64_t. For mod_intset_u64s, the size of the array is passed as
            a parameter n. For mod_intget_u64s(), the required array size can be 
            determined via mod_intbits(); if the Integer is zero, mod_intget_u64s()
            writes 0 to the first output uint64_t. It returns the number of 
            uint64_t written. */
        // void set(uint64_t s) {v[0] = s;}
        void set(const uint64_t *s, const size_t n) {
            ASSERT_ALWAYS(n <= 1);
            if (n == 0)
                v[0] = 0;
            else
                v[0] = s[0];
        }
        uint64_t get_u64() const {
            return v[0];
        }
        /* Return the size in uint64_ts that is required in the output for get_u64s() */
        size_t size() const {return 1;}
        size_t get_u64s (uint64_t *r) const {
            r[0] = v[0];
            return 1;
        }

        bool operator==(const Integer a) const {return (v[0] == a.v[0]);}
        bool operator!=(const Integer a) const {return (v[0] != a.v[0]);}
        bool operator<(const Integer a) const {return (v[0] < a.v[0]);}
        bool operator>(const Integer a) const {return (v[0] > a.v[0]);}
        bool operator<=(const Integer a) const {return (v[0] <= a.v[0]);}
        bool operator>=(const Integer a) const {return (v[0] >= a.v[0]);}
        /* x.cmp(a) returns -1 for x<a, 0 for x==a, 1 for x>a */
        int cmp(const Integer a) const {return (*this < a) ? -1 : (*this == a) ? 0 : 1;}
        
        /* Binary operators like operator+ would be nice here but for
         * variable-size operands which are stored on the heap would
         * also require allocating temporaries, whereas most of the
         * operations here can work in-place. */
        Integer operator+=(const Integer a) {return Integer(v[0] += a.v[0]);}
        Integer operator-=(const Integer a) {return Integer(v[0] -= a.v[0]);}
        Integer operator>>=(const int i) {return Integer(v[0] >>= i);}
        Integer operator<<=(const int i) {return Integer(v[0] <<= i);}
        Integer operator/=(const Integer a) {return Integer(v[0] /= a.v[0]);}
        Integer operator%=(const Integer a) {return Integer(v[0] %= a.v[0]);}
        /* r = v/a. We require a|v. */
        Integer divexact(const Integer a) const {ASSERT_EXPENSIVE(v[0] % a.v[0] == 0); return Integer(v[0] / a.v[0]);}
        
        /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
        size_t bits() const {
            if (v[0] == 0)
                return 0;
            return 64 - u64arith_clz (v[0]);
        }
    };

private:
  uint64_t m;
public:
    Modulus64(const uint64_t s) : m(s){}
    Modulus64(const Modulus64 &s) : m(s.m){}
    Modulus64(const Integer s) {s.get_u64s(&m);}
    ~Modulus64() {}
  uint64_t getmod_u64 () const {return m;}
  void getmod_int (Integer &r) const {r = m;}

  /* Functions for residues */
  class Residue {
  private:
    friend class Modulus64;
    uint64_t r;
#if 0
    /* Allow simple assignment syntax */
    Residue &operator=(const uint64_t s){r = s;}
    Residue &operator=(const Residue s){r = s.r;}
#endif
  };

private:
    void assertValid(const Residue a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a.r < m);}
    void assertValid(const uint64_t a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a < m);}
public:
  void init (Residue &r) const {r.r = 0;}
  void init_noset0 (Residue &r MAYBE_UNUSED) const {}
  void clear (Residue &r MAYBE_UNUSED) const {}
  void set (Residue &r, const Residue s) const {assertValid(s); r = s;}
  void set_u64 (Residue &r, const uint64_t s) const {r.r = s % m;}
  void set_int (Residue &r, const Integer s) const {r.r = s.get_u64() % m;}
  /* Sets the Residue to the class represented by the integer s. Assumes that
     s is reduced (mod m), i.e. 0 <= s < m */
  void set_u64_reduced (Residue &r, const uint64_t s) const {assertValid(s); r.r = s;}
  void set_int_reduced (Residue &r, const Integer s) const {assertValid(s.get_u64()); r.r = s.get_u64();}
  void set_int64 (Residue &r, const int64_t s) {r.r = llabs(s) % m; if (s < 0) neg(r, r);}
  void set0 (Residue &r) const {r.r = 0;}
  void set1 (Residue &r) const {r.r = (m != 1);}
  void swap (Residue &a, Residue &b) const {uint64_t t = a.r; a.r = b.r; b.r = t;}
  uint64_t get_u64 (const Residue s) const {assertValid(s); return s.r;}
  void get_int (Integer &r, const Residue s) const {assertValid(s); r = s.r;}
  int equal (const Residue a, const Residue b) const {assertValid(a); assertValid(b); return (a.r == b.r);}
  int is0 (const Residue a) const {assertValid(a); return (a.r == 0);}
  int is1 (const Residue a) const {assertValid(a); return (a.r == 1);}
  void neg (Residue &r, const Residue a) const {
    assertValid(a);
    if (a.r == 0)
      r.r = a.r;
    else
      r.r = m - a.r;
  }
  void add (Residue &r, const Residue a, const Residue b) const {
      u64arith_addmod_1_1(&r.r, a.r, b.r, m);
  }
  void add1 (Residue &r, const Residue a) const {
    assertValid(a);
    r.r = a.r + 1;
    if (r.r == m)
      r.r = 0;
  }
  void add_u64 (Residue &r, const Residue a, const uint64_t b) const {
    u64arith_addmod_1_1(&r.r, a.r, b % m, m);
  }
  void sub (Residue &r, const Residue a, const Residue b) const {
    u64arith_submod_1_1(&r.r, a.r, b.r, m);
  }
  void sub_u64 (Residue &r, const Residue a, const uint64_t b) const {
    u64arith_submod_1_1(&r.r, a.r, b % m, m);
  }
  void mul (Residue &r, const Residue a, const Residue b) const {
    uint64_t t1, t2;
    assertValid(a);
    assertValid(b);
    u64arith_mul_1_1_2 (&t1, &t2, a.r, b.r);
    u64arith_divr_2_1_1 (&r.r, t1, t2, m);
  }
  void sqr (Residue &r, const Residue a) const {
    uint64_t t1, t2;
    assertValid(a);
    u64arith_mul_1_1_2 (&t1, &t2, a.r, a.r);
    u64arith_divr_2_1_1 (&r.r, t1, t2, m);
  }
  /* Computes (a * 2^wordsize) % m */
  void tomontgomery (Residue &r, const Residue a) const {
    assertValid(a);
    u64arith_divr_2_1_1 (&r.r, 0, a.r, m);
  }
  /* Computes (a / 2^wordsize) % m */
  void frommontgomery (Residue &r, const Residue a, const uint64_t invm) const {
    uint64_t tlow, thigh;
    assertValid(a);
    tlow = a.r * invm;
    u64arith_mul_1_1_2 (&tlow, &thigh, tlow, m);
    r.r = thigh + (a.r != 0 ? 1 : 0);
  }
  /* Computes (a / 2^wordsize) % m, but result can be r = m. 
     Input a must not be equal 0 */
  void redcsemi_u64_not0 (Residue &r, const uint64_t a, const uint64_t invm) const {
    uint64_t tlow, thigh;
    ASSERT (a != 0);
    tlow = a * invm; /* tlow <= 2^w-1 */
    u64arith_mul_1_1_2 (&tlow, &thigh, tlow, m);
    /* thigh:tlow <= (2^w-1) * m */
    r.r = thigh + 1; 
    /* (thigh+1):tlow <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m 
                      <= 2^w * (m + 1) - m */
    /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
         <= m */
  }
  void div2 (Residue &r, const Residue a) const {r.r = u64arith_div2mod(a.r, m);}
  int next (Residue &r) const {return (++r.r == m);}
  int finished (const Residue r) const {return (r.r == m);}

  /* prototypes of non-inline functions */
  int div3 (Residue &, const Residue) const;
  int div5 (Residue &, const Residue) const;
  int div7 (Residue &, const Residue) const;
  int div11 (Residue &, const Residue) const;
  int div13 (Residue &, const Residue) const;
  void gcd (Integer &, const Residue) const;
  void pow_u64 (Residue &, const Residue, const uint64_t) const;
  void pow2_u64 (Residue &, const uint64_t) const;
  void pow_mp (Residue &, const Residue, const uint64_t *, const int) const;
  void pow2_mp (Residue &, const uint64_t *, const int) const;
  void V_u64 (Residue &, const Residue, const uint64_t) const;
  void V_mp (Residue &, const Residue, const uint64_t *, const int) const;
  int sprp (const Residue) const;
  int sprp2 () const;
  int isprime () const;
  int inv (Residue &, const Residue) const;
  int inv_odd (Residue &, const Residue) const;
  int inv_powerof2 (Residue &, const Residue) const;
  int batchinv (Residue *, const Residue *, size_t, const Residue *) const;
  int jacobi (const Residue) const;
  void pow3_u64 (Residue &, uint64_t) const;
private:
  int find_minus1 (Residue r1, const Residue minusone, const int po2) const;
};

#endif  /* MOD_64_H */
