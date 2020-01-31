#ifndef EC_ARITH_WEIERSTRASS_H_
#define EC_ARITH_WEIERSTRASS_H_

#include <iostream>


/* Short Weierstrass elliptic curves
 *
 * Affine coordinates, with equation
 *                   y^2 = x^3 + a*x + b
 * Projective coordinates, with equation:
 *                  Y^2*Z = X^3 + a*X*Z^Z + b*Z^3
 *
 * Curve coefficient needed in computation: a
 */

template <typename MODULUS>
class ECWeierstrass
{
public:
    typedef MODULUS Modulus;
    typedef typename Modulus::Residue Residue;
    typedef typename Modulus::Integer Integer;
    const Modulus &m;
    Residue a;
    ECWeierstrass(const Modulus &m, Residue &_a) : m(m), a(m) {m.set(a, _a);}

    void print (std::ostream &os, const char *prefix) const
    {
        const char *pre = (prefix == NULL) ? "" : prefix;
        Integer A, M;
        m.get (A, a);
        m.getmod(M);
        os << pre << "Weierstrass curve: y^2 = x^3 + a*x + b" << std::endl
           << pre << "a = " << A << " over Z/mZ, m = " << M << std::endl;
    }
    friend std::ostream & operator<<(std::ostream &os, const ECWeierstrass<MODULUS> &c) {
        c.print(os, NULL);
        return os;
    }

};

template <typename MODULUS>
class ECWeierstrassAffinePoint
{
    typedef MODULUS Modulus;
    typedef ECWeierstrass<MODULUS> ECurve;
    typedef ECWeierstrassAffinePoint Point;
    typedef typename Modulus::Residue Residue;
    typedef typename Modulus::Integer Integer;
    const ECurve &curve;
    Residue x,y;
    bool finite; /* If finite is false, then the point is the point at
                    infinity and the x, y coordinates are not considered
                    for any arithmetic. */

public:    
    ECWeierstrassAffinePoint(const ECurve &c) : curve(c), x(curve.m), y(curve.m), finite(true) {}
    ECWeierstrassAffinePoint(const ECurve &c, const Residue &_x, const Residue &_y) 
    : curve(c), x(curve.m, _x), y(curve.m, _y) {finite = true;}
    ECWeierstrassAffinePoint(const ECWeierstrassAffinePoint &s)
    : curve(s.curve), x(curve.m, s.x), y(curve.m, s.y), finite(s.finite) {}
    ECWeierstrassAffinePoint(ECWeierstrassAffinePoint &&) = default;
    ~ECWeierstrassAffinePoint() {}

    Point &operator=(const Point &other) {
        ASSERT_EXPENSIVE(&curve == &other.curve);
        if (this != &other) {
            curve.m.set(x, other.x);
            curve.m.set(y, other.y);
            finite = other.finite;
        }
        return *this;
    }

    bool operator==(const Point &other) const {
        ASSERT_EXPENSIVE(&curve == &other.curve);
        if (finite != other.finite)
            return false;
        if (!finite)
            return true;
        return curve.m.equal(x, other.x) && curve.m.equal(y, other.y);
    }

    bool operator!=(const Point &other) const {
        return !(*this == other);
    }

    Point operator+(const Point &Q) const {
        Point R(curve);
        add(R, Q);
        return R;
    }

    Point & operator+=(const Point &Q) {
        add(*this, Q);
        return *this;
    }

    Point operator*(const uint64_t e) const {
        Point R(curve);
        smul(R, e);
        return R;
    }

    Point & operator*=(const uint64_t e) {
        smul(*this, e);
        return *this;
    }

    void set(const Residue &_x, const Residue &_y) {
        curve.m.set(x, _x);
        curve.m.set(y, _y);
        finite = true;
    }
    
    void set0() {
        finite = false;
    }
    
    bool is0() const {
        return !finite;
    }
    
    void swap(Point &other) {
        ASSERT_EXPENSIVE(&curve == &other.curve);
        curve.m.swap (x, other.x);
        curve.m.swap (y, other.y);
        std::swap(finite, other.finite);
    }

    void print(std::ostream &os) const {
        if (finite) {
            Integer X, Y;

            curve.m.get (X, x);
            curve.m.get (Y, y);
            os << "(" << X << " : " << Y << ")";
        } else {
            os << "(point at infinity)";
        }
    }

    void dbl (Point &R) const;
    void add (Point &R, const Point &Q) const;
    void smul (Point &R, const uint64_t e) const;
    uint64_t point_order (uint64_t known_m, uint64_t known_r, int verbose) const;

    friend std::ostream & operator<<(std::ostream &os, const ECWeierstrassAffinePoint<MODULUS> &p) {
        p.print(os);
        return os;
    }
};

template <typename MODULUS>
class ECWeierstrassProjectivePoint
{
    typedef MODULUS Modulus;
    typedef ECWeierstrass<MODULUS> ECurve;
    typedef ECWeierstrassProjectivePoint Point;
    typedef typename Modulus::Residue Residue;
    typedef typename Modulus::Integer Integer;
    const ECurve &curve;
    Residue x, y, z;

    ECWeierstrassProjectivePoint(const ECurve &c) : curve(c), x(curve.m), y(curve.m), z(curve.m) {}
    ECWeierstrassProjectivePoint(const ECurve &c, const Residue &_x, const Residue &_y, const Residue &_z) 
    : curve(c), x(curve.m), y(curve.m), z(curve.m) {
        curve.m.set(x, _x);
        curve.m.set(y, _y);
        curve.m.set(z, _z);
    }
    ECWeierstrassProjectivePoint(const ECWeierstrassProjectivePoint &s)
    : curve(s.curve), x(curve.m, s.x), y(curve.m, s.y), z(curve.m, s.z) {}
    ECWeierstrassProjectivePoint(ECWeierstrassProjectivePoint &&) = default;
    ~ECWeierstrassProjectivePoint() {}

    Point &operator=(const Point &other) {
        ASSERT_EXPENSIVE(&curve == &other.curve);
        if (this != &other) {
            curve.m.set(x, other.x);
            curve.m.set(y, other.y);
            curve.m.set(z, other.z);
        }
        return *this;
    }
    bool operator==(const Point &other) const {
        Residue t1(curve.m), t2(curve.m);
        curve.m.mul(t1, x, other.z);
        curve.m.mul(t2, z, other.x);
        if (!curve.m.equal(t1, t2))
            return false;
        curve.m.mul(t1, y, other.z);
        curve.m.mul(t2, z, other.y);
        if (!curve.m.equal(t1, t2))
            return false;
        return true;
    }

    bool operator!=(const Point &other) const {
        return !(*this == other);
    }

    Point operator+(const Point &Q) const {
        Point R(curve);
        add(R, Q);
        return R;
    }

    Point & operator+=(const Point &Q) {
        add(*this, Q);
        return *this;
    }

    Point operator*(const uint64_t e) const {
        Point R(curve);
        smul(R, e);
        return R;
    }

    Point & operator*=(const uint64_t e) {
        smul(*this, e);
        return *this;
    }

    void set(const Residue &_x, const Residue &_y, const Residue &_z) {
        curve.m.set(x, _x);
        curve.m.set(y, _y);
        curve.m.set(z, _z);
    }
    
    void swap(Point &other) {
        ASSERT_EXPENSIVE(&curve == &other.curve);
        curve.m.swap (x, other.x);
        curve.m.swap (y, other.y);
        curve.m.swap (z, other.z);
    }

    void print(std::ostream &os) const {
        Integer X, Y, Z;

        curve.m.get (X, x);
        curve.m.get (Y, y);
        curve.m.get (Z, z);
        os << "(" << X << " : " << Y << " : " << Z << ")";
    }

    /* Set P to zero (the neutral point): (0:1:0) */
    void set0 () {
        curve.m.set0 (x);
        curve.m.set1 (y);
        curve.m.set0 (z);
    }

    bool is0() const {
        return curve.m.is0(x) && curve.m.is1 (y) && curve.m.is0 (z);
    }

    void dbl (Point &R) const;
    void add (Point &R, const Point &Q) const;
    void smul (Point &R, const uint64_t e) const;
};


#if 0
/* Convert the Montgomery curve B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2 with a valid
 * point Pm into a affine Weierstrass curve y^2 = x^3 + a*x + b with a
 * valid affine point Pw.
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in Pw->x.
 *
 * The curve coefficient b of the short Weierstrass curve will not be computed.
 * a and A can be the same variable.
 * Pm and Pw can be the same variable.
 */
#define weierstrass_aff_from_montgomery MOD_APPEND_TYPE(weierstrass_aff_from_montgomery)
MAYBE_UNUSED
static int
weierstrass_aff_from_montgomery (residue_t a, ec_point_t Pw, const residue_t A,
                                 ec_point_t Pm, const modulus_t m)
{
  residue_t B, one, t, x;
  int ret;

  mod_init_noset0 (B, m);
  mod_init_noset0 (one, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (x, m);
  mod_set1 (one, m);

  ret = mod_inv (t, Pm->z, m);
  if (ret == 0)
  {
    fprintf (stderr, "%s: could not invert Z\n", __func__);
    mod_set (Pw->x, Pm->z, m);
  }
  else
  {
    mod_mul (x, Pm->x, t, m); /* x = X/Z */
    mod_add (B, x, A, m);
    mod_mul (B, B, x, m);
    mod_add (B, B, one, m);
    mod_mul (B, B, x, m); /* B = x^3 + A*x^2 + x */

    /* Now (x,1) is on the curve B*y^2 = x^3 + A*x^2 + x. */
    ret = mod_inv (Pw->y, B, m);    /* y = 1/B */
    if (ret == 0)
    {
      fprintf (stderr, "%s: could not invert B\n", __func__);
      mod_set (Pw->x, B, m);
    }
    else
    {
      mod_div3 (t, A, m);
      mod_add (Pw->x, x, t, m);
      mod_mul (Pw->x, Pw->x, Pw->y, m); /* x = (X + A/3)/B */
      mod_mul (a, t, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, Pw->y, m);
      mod_mul (a, a, Pw->y, m);         /* a = (1 - (A^2)/3)/B^2 */
    }
  }

  mod_clear (one, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (x, m);

  return ret;
}

#endif

#endif /* EC_ARITH_WEIERSTRASS_H_ */
