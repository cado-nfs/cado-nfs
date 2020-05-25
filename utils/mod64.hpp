/* A class for modular arithmetic with residues and modulus of up to 64
 * bits. */

#ifndef MOD64_HPP
#define MOD64_HPP

/**********************************************************************/
#include <cstdint>
#include <cstdlib>
#include <new>
#include "macros.h"
#include "u64arith.h"
#include "modint.hpp"
#include "mod_stdop.hpp"

class Modulus64 {
    /* Type definitions */
public:
    typedef Integer64 Integer;
    class Residue {
        friend class Modulus64;
    protected:
        uint64_t r;
    public:
        typedef Modulus64 Modulus;
        typedef Modulus::Integer Integer;
        typedef bool IsResidueType;
        Residue() = delete;
        Residue(const Modulus &m MAYBE_UNUSED) : r(0) {}
        Residue(const Modulus &m MAYBE_UNUSED, const Residue &s) : r(s.r) {}
        Residue(const Residue &&s) : r(s.r) {}
    protected:
        Residue &operator=(const Residue &s) {r = s.r; return *this;}
        Residue &operator=(const Integer &s) {r = 0; s.get(&r, 1); return *this;}
        Residue &operator=(const uint64_t s) {r = s; return *this;}
    };

    typedef ResidueStdOp<Residue> ResidueOp;

    /* Data members */
protected:
    uint64_t m;
    /* Methods used internally */
    void assertValid(const Residue &a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a.r < m);}
    void assertValid(const uint64_t a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a < m);}
    uint64_t get_u64 (const Residue &s) const {assertValid(s); return s.r;}

    /* Methods of the API */
public:
    static Integer getminmod() {return Integer(1);}
    static Integer getmaxmod() {return Integer(UINT64_MAX);}
    static void getminmod(Integer &r) {r = getminmod();}
    static void getmaxmod(Integer &r) {r = getmaxmod();}
    static bool valid(const Integer &m) {return getminmod() <= m && m <= getmaxmod();}
    
    Modulus64(const uint64_t s) : m(s){}
    Modulus64(const Modulus64 &s) : m(s.m){}
    Modulus64(const Integer &s) {s.get(&m, 1);}
    ~Modulus64() {}
    uint64_t getmod_u64 () const {return m;}
    void getmod (Integer &r) const {r = m;}

    /* Methods for residues */

    /** Allocate an array of len residues.
     *
     * Must be freed with deleteArray(), not with delete[].
     */
    Residue *newArray(const size_t len) const {
        void *t = operator new[](len * sizeof(Residue));
        if (t == NULL)
            return NULL;
        Residue *ptr = static_cast<Residue *>(t);
        for(size_t i = 0; i < len; i++) {
            new(&ptr[i]) Residue(*this);
        }
        return ptr;
    }

    void deleteArray(Residue *ptr, const size_t len) const {
        for(size_t i = len; i > 0; i++) {
            ptr[i - 1].~Residue();
        }
        operator delete[](ptr);
    }

    void set (Residue &r, const Residue &s) const {assertValid(s); r = s;}
    void set (Residue &r, const uint64_t s) const {r.r = s % m;}
    void set (Residue &r, const Integer &s) const {s.get(&r.r, 1); r.r %= m;}
    /* Sets the Residue to the class represented by the integer s. Assumes that
       s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced (Residue &r, const uint64_t s) const {assertValid(s); r.r = s;}
    void set_reduced (Residue &r, const Integer &s) const {s.get(&r.r, 1); assertValid(r);}
    void set_int64 (Residue &r, const int64_t s) const {r.r = llabs(s) % m; if (s < 0) neg(r, r);}
    void set0 (Residue &r) const {r.r = 0;}
    void set1 (Residue &r) const {r.r = (m != 1);}
    /* Exchanges the values of the two arguments */
    void swap (Residue &a, Residue &b) const {uint64_t t = a.r; a.r = b.r; b.r = t;}
    void get (Integer &r, const Residue &s) const {assertValid(s); r = Integer(s.r);}
    bool equal (const Residue &a, const Residue &b) const {assertValid(a); assertValid(b); return (a.r == b.r);}
    bool is0 (const Residue &a) const {assertValid(a); return (a.r == 0);}
    bool is1 (const Residue &a) const {assertValid(a); return (a.r == 1);}
    void neg (Residue &r, const Residue &a) const {
        assertValid(a);
        if (a.r == 0)
            r.r = a.r;
        else
            r.r = m - a.r;
    }
  void add (Residue &r, const Residue &a, const Residue &b) const {u64arith_addmod_1_1(&r.r, a.r, b.r, m);}
  void add1 (Residue &r, const Residue &a) const {
    assertValid(a);
    r.r = a.r + 1;
    if (r.r == m)
      r.r = 0;
  }
  void add (Residue &r, const Residue &a, const uint64_t b) const {
    u64arith_addmod_1_1(&r.r, a.r, b % m, m);
  }
  void sub (Residue &r, const Residue &a, const Residue &b) const {
    u64arith_submod_1_1(&r.r, a.r, b.r, m);
  }
  void sub1 (Residue &r, const Residue &a) const {
    u64arith_submod_1_1(&r.r, a.r, 1, m);
  }
  void sub (Residue &r, const Residue &a, const uint64_t b) const {
    u64arith_submod_1_1(&r.r, a.r, b % m, m);
  }
  void mul (Residue &r, const Residue &a, const Residue &b) const {
    uint64_t t1, t2;
    assertValid(a);
    assertValid(b);
    u64arith_mul_1_1_2 (&t1, &t2, a.r, b.r);
    u64arith_divr_2_1_1 (&r.r, t1, t2, m);
  }
  void sqr (Residue &r, const Residue &a) const {
    uint64_t t1, t2;
    assertValid(a);
    u64arith_mul_1_1_2 (&t1, &t2, a.r, a.r);
    u64arith_divr_2_1_1 (&r.r, t1, t2, m);
  }
  /* Computes (a * 2^wordsize) % m */
  void tomontgomery (Residue &r, const Residue &a) const {
    assertValid(a);
    u64arith_divr_2_1_1 (&r.r, 0, a.r, m);
  }
  /* Computes (a / 2^wordsize) % m */
  void frommontgomery (Residue &r, const Residue &a, const uint64_t invm) const {
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
  bool next (Residue &r) const {return (++r.r == m);}
  bool finished (const Residue &r) const {return (r.r == m);}
  bool div2 (Residue &r, const Residue &a) const {
      if (m % 2 == 0)
          return false;
      else {
          r.r = u64arith_div2mod(a.r, m);
          return true;
      }
    }

  /* Given a = V_n (x), b = V_m (x) and d = V_{n-m} (x), compute V_{m+n} (x).
   * r can be the same variable as a or b but must not be the same variable as d.
   */
  void V_dadd (Residue &r, const Residue &a, const Residue &b,
               const Residue &d) const {
    ASSERT (&r != &d);
    mul (r, a, b);
    sub (r, r, d);
  }

  /* Given a = V_n (x) and two = 2, compute V_{2n} (x).
   * r can be the same variable as a but must not be the same variable as two.
   */
  void V_dbl (Residue &r, const Residue &a, const Residue &two) const {
    ASSERT (&r != &two);
    sqr (r, a);
    sub (r, r, two);
  }

  /* prototypes of non-inline functions */
  bool div3 (Residue &, const Residue &) const;
  bool div5 (Residue &, const Residue &) const;
  bool div7 (Residue &, const Residue &) const;
  bool div11 (Residue &, const Residue &) const;
  bool div13 (Residue &, const Residue &) const;
  void gcd (Integer &, const Residue &) const;
  void pow (Residue &, const Residue &, const uint64_t) const;
  void pow (Residue &, const Residue &, const uint64_t *, const size_t) const;
  void pow (Residue &, const Residue &, const Integer &) const;
  void pow2 (Residue &, const uint64_t) const;
  void pow2 (Residue &, const uint64_t *, const size_t) const;
  void pow2 (Residue &, const Integer &) const;
  void pow3 (Residue &, uint64_t) const;
  void V (Residue &, const Residue &, const uint64_t) const;
  void V (Residue &, const Residue &, const uint64_t *, const size_t) const;
  void V (Residue &, const Residue &, const Integer &) const;
  void V (Residue &r, Residue *rp1, const Residue &b,
          const uint64_t k) const;
  bool sprp (const Residue &) const;
  bool sprp2 () const;
  bool isprime () const;
  bool inv (Residue &, const Residue &) const;
  bool inv_odd (Residue &, const Residue &) const;
  bool inv_powerof2 (Residue &, const Residue &) const;
  bool batchinv (Residue *, const Residue *, size_t, const Residue *) const;
  int jacobi (const Residue &) const;
protected:
  bool find_minus1 (Residue &r1, const Residue &minusone, const int po2) const;
};

#endif  /* MOD64_HPP */
