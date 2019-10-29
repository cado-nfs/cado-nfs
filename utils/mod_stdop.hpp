#ifndef MOD_STDOP_HPP
#define MOD_STDOP_HPP

/* Functions for residues. Slow and clean
  This class extends the Residue type of an existing Modulus type by adding
  standard C++ operators for the modular arithmetic. For this to work, each
  residue needs a reference to the modulus w.r.t. which it was declared, as
  there is no other way of passing the modulus to the arithmetic operators.
 */
template <typename T, typename = typename T::IsResidueType>
class ResidueStdOp : public T {
  typedef T Residue;
  typedef typename T::Modulus Modulus;
  typedef typename T::Integer Integer;
private:
  const Modulus &m;
  void assertValid(const ResidueStdOp &other MAYBE_UNUSED) const {ASSERT_EXPENSIVE(&m = &other.m);}
public:
  ResidueStdOp(const Modulus &m) : T(m), m(m) {}
  ResidueStdOp(const Modulus &m, const Residue &s) : T(m), m(m) {m.set(*this, s);}
  ResidueStdOp(const Modulus &m, const Integer &s) : T(m), m(m) {m.set(*this, s);}
  ResidueStdOp(const Modulus &m, const uint64_t &s) : T(m), m(m) {m.set(*this, s);}
  ResidueStdOp(const ResidueStdOp &s) : Residue(s.m), m(s.m) {*this = s;}
  ResidueStdOp(ResidueStdOp &&s) : Residue(s.m), m(s.m) {*this = s;}
  
  ResidueStdOp &operator=(const ResidueStdOp &s) {assertValid(s); m.set(*this, s); return *this;}
  ResidueStdOp &operator=(const uint64_t s) {m.set(*this, s); return *this;}  
  ResidueStdOp &operator=(const Integer &s) {m.set(*this, s); return *this;}  
  
  bool operator==(const ResidueStdOp &s) const {assertValid(s); return m.equal(*this, s);}
  bool operator!=(const ResidueStdOp &s) const {assertValid(s); return !m.equal(*this, s);}
  bool is0() const {return m.is0(*this);}
  bool is1() const {return m.is1(*this);}

  operator Integer() {Integer r; m.get(r, *this); return r;}
  
  ResidueStdOp operator-() const {ResidueStdOp r = *this; m.neg(r, r); return r;}
  /* Prefix ++ and -- */
  ResidueStdOp &operator++() {m.add1(*this, *this); return *this;}
  ResidueStdOp &operator--() {m.sub1(*this, *this); return *this;}
  /* Postfix ++ and -- */
  ResidueStdOp operator++(int) {ResidueStdOp r = *this; m.add1(*this, *this); return r;}
  ResidueStdOp operator--(int) {ResidueStdOp r = *this; m.sub1(*this, *this); return r;}

  ResidueStdOp operator+(const ResidueStdOp &s) const {assertValid(s); ResidueStdOp r(m); m.add(r, *this, s); return r;}
  ResidueStdOp operator+(const uint64_t s) const {ResidueStdOp r(m, s); m.add(r, *this, r); return r;}
  ResidueStdOp operator+(const Integer &s) const {ResidueStdOp r(m, s); m.add(r, *this, r); return r;}
  ResidueStdOp &operator+=(const ResidueStdOp &s) {assertValid(s); m.add(*this, *this, s); return *this;}
  ResidueStdOp &operator+=(const uint64_t s) {ResidueStdOp t(m, s); m.add(*this, *this, t); return *this;}
  ResidueStdOp &operator+=(const Integer &s) {ResidueStdOp t(m, s); m.add(*this, *this, t); return *this;}
  friend ResidueStdOp operator+(const uint64_t s, const ResidueStdOp &t) {return t + s;}
  friend ResidueStdOp operator+(const Integer &s, const ResidueStdOp &t) {return t + s;}
  friend Integer &operator+=(Integer &s, const ResidueStdOp &t) {ResidueStdOp r(t.m, s); t.m.add(r, r, t); s = (Integer) r; return s;}

  ResidueStdOp operator-(const ResidueStdOp &s) const {assertValid(s); ResidueStdOp r(m); m.sub(r, *this, s); return r;}
  ResidueStdOp operator-(const uint64_t s) const {ResidueStdOp r(m, s); m.sub(r, *this, r); return r;}
  ResidueStdOp operator-(const Integer &s) const {ResidueStdOp r(m, s); m.sub(r, *this, r); return r;}
  ResidueStdOp &operator-=(const ResidueStdOp &s) {assertValid(s); m.sub(*this, *this, s); return *this;}
  ResidueStdOp &operator-=(const uint64_t s) {ResidueStdOp t(m, s); m.sub(*this, *this, t); return *this;}
  ResidueStdOp &operator-=(const Integer &s) {ResidueStdOp t(m, s); m.sub(*this, *this, t); return *this;}
  friend ResidueStdOp operator-(const uint64_t s, const ResidueStdOp &t) {return ResidueStdOp(t.m, s) - t;}
  friend ResidueStdOp operator-(const Integer &s, const ResidueStdOp &t) {return ResidueStdOp(t.m, s) - t;}
  friend Integer &operator-=(Integer &s, const ResidueStdOp &t) {return s = (Integer) (ResidueStdOp(t.m, s) - t);}

  ResidueStdOp operator*(const ResidueStdOp &s) const {
      assertValid(s);
      ResidueStdOp r(m);
      if (CONSTANT_P(this == &s) && this == &s) {
          m.sqr(r, *this);
      } else {
        m.mul(r, *this, s);
      }
    return r;
  }
  ResidueStdOp operator*(const uint64_t s) const {ResidueStdOp r(m, s); m.mul(r, *this, r); return r;}
  ResidueStdOp operator*(const Integer &s) const {ResidueStdOp r(m, s); m.mul(r, *this, r); return r;}
  ResidueStdOp &operator*=(const ResidueStdOp &s) {assertValid(s); m.mul(*this, *this, s); return *this;}
  ResidueStdOp &operator*=(const uint64_t s) {ResidueStdOp t(m, s); m.mul(*this, *this, t); return *this;}
  ResidueStdOp &operator*=(const Integer &s) {ResidueStdOp t(m, s); m.mul(*this, *this, t); return *this;}
  friend ResidueStdOp operator*(const uint64_t s, const ResidueStdOp &t) {ResidueStdOp r(t.m, s); t.m.mul(r, r, t); return r;}
  friend ResidueStdOp operator*(const Integer &s, const ResidueStdOp &t) {ResidueStdOp r(t.m, s); t.m.mul(r, r, t); return r;}
  friend Integer &operator*=(Integer &s, const ResidueStdOp &t) {ResidueStdOp r(t.m, s); t.m.mul(r, r, t); s = r; return s;}
  
  ResidueStdOp operator/(const ResidueStdOp &s) const {assertValid(s); ResidueStdOp r(m); m.inv(r, s); m.mul(r, *this, r); return r;}
  ResidueStdOp &operator/=(const ResidueStdOp &s) {assertValid(s); ResidueStdOp i(m); m.inv(i, s); m.mul(*this, *this, i); return *this;}
  ResidueStdOp operator/(const uint64_t s) const {
      ResidueStdOp r(m);
      if (CONSTANT_P(s) && s == 2) {
          m.div2(r, *this);
      } else if (CONSTANT_P(s) && s == 3) {
          m.div3(r, *this);
      } else if (CONSTANT_P(s) && s == 5) {
          m.div5(r, *this);
      } else if (CONSTANT_P(s) && s == 7) {
          m.div7(r, *this);
      } else if (CONSTANT_P(s) && s == 11) {
          m.div11(r, *this);
      } else if (CONSTANT_P(s) && s == 13) {
          m.div13(r, *this);
      } else {
          ResidueStdOp r(m);
          m.set(r, s);
          m.inv(r, r);
          m.mul(r, *this, r);
      }
      return r;
  }
  ResidueStdOp &operator/=(const uint64_t s) {return *this /= ResidueStdOp(m, s);}
  friend ResidueStdOp operator/(const uint64_t s, const ResidueStdOp &t) {
      if (CONSTANT_P(s) && s == 1) { /* For 1/t, don't multiply by 1 */
          ResidueStdOp r(t.m);
          t.m.inv(r, t);
          return r;
      } else {
          return ResidueStdOp(t.m, s) / t;
      }
  }
  ResidueStdOp operator/(const Integer &s) const {return *this / ResidueStdOp(m, s);}
  ResidueStdOp &operator/=(const Integer &s) {return *this /= ResidueStdOp(m, s);}
  friend ResidueStdOp operator/(const Integer &s, const ResidueStdOp &t) {return ResidueStdOp(t.m, s) / t;}
  friend Integer &operator/=(Integer &s, const ResidueStdOp &t) {return s = (Integer) (ResidueStdOp(t.m, s) / t);}
  
  ResidueStdOp pow(const uint64_t e) {ResidueStdOp r(m); m.pow_u64(r, *this, e); return r;}
  ResidueStdOp pow(const Integer &e) {ResidueStdOp r(m); m.pow(r, *this, e); return r;}
  ResidueStdOp pow(const uint64_t *e, const size_t l) {ResidueStdOp r(m); m.pow(r, *this, e, l); return r;}
  ResidueStdOp chebyshevV(const uint64_t e) {ResidueStdOp r(m); m.V(r, *this, e); return r;}
  ResidueStdOp chebyshevV(const Integer &e) {ResidueStdOp r(m); m.V(r, *this, e); return r;}
  ResidueStdOp chebyshevV(const uint64_t *e, const size_t l) {ResidueStdOp r(m); m.V(r, *this, e, l); return r;}
  
  /* Should we call it index or gcd? I can't decide */
  Integer index() const {Integer g; m.gcd(g, *this); return g;}
  Integer gcd() const {Integer g; m.gcd(g, *this); return g;}
  int jacobi() const {return m.jacobi(*this);}
};

#endif /* MOD_STDOP_HPP */
