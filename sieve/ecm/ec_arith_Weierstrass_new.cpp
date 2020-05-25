#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include "cxx_mpz.hpp"
// IWYU pragma: no_include "modint.hpp"
#include <cstdlib>        // for abort, NULL
#include <cinttypes>       // for PRIu64
#include <cmath>           // for sqrt, ceil
#include <vector>          // for vector

#include "macros.h"        // for ASSERT
// #include "modint.hpp"      // for operator<<

#include "mod64.hpp"    // IWYU pragma: keep
#include "modredc64.hpp"        // IWYU pragma: keep
#include "modredc126.hpp"       // IWYU pragma: keep
#include "mod_mpz_new.hpp"      // IWYU pragma: keep
#include "ec_arith_Weierstrass_new.hpp"

/* Computes R=2P, with 1 inv, 4 muls (2 muls and 2 squares) and 8 add/sub.
 *
 * If result is point at infinity, return non-invertible value in R.x.
 *
 * It is permissible to let *this and R use the same memory.
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::AffinePoint::dbl (AffinePoint &R) const {
    ASSERT_EXPENSIVE(curve.is_same(R.curve));
    Residue lambda(curve.m), u(curve.m), v(curve.m);

    if (is0()) {
        R.set0();
        return;
    }
    
    curve.m.sqr (u, x);
    curve.m.add (v, u, u);
    curve.m.add (v, v, u);
    curve.m.add (v, v, curve.a); /* 3x^2 + a */
    curve.m.add (u, y, y);
    bool ret = curve.m.inv (lambda, u);    /* 1/(2*y) */
    if (ret) {
        curve.m.mul (lambda, lambda, v);
        curve.m.sqr (u, lambda);
        curve.m.sub (u, u, x);
        curve.m.sub (u, u, x);    /* x3 = u = lambda^2 - 2*x */
        curve.m.sub (v, x, u);
        curve.m.mul (v, v, lambda);
        curve.m.sub (R.y, v, y);
        curve.m.set (R.x, u);
    } else {
        curve.m.set (R.x, u);
    }
    R.finite = ret;
}

/* Computes R=P+Q, with 1 inv, 3 muls (2 muls and 1 square) and 6 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return true if it worked, false if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in R->x.
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::AffinePoint::add (AffinePoint &R, const AffinePoint &Q) const {
    ASSERT_EXPENSIVE(curve.is_same(Q.curve));
    ASSERT_EXPENSIVE(curve.is_same(R.curve));
    Residue lambda(curve.m), u(curve.m), v(curve.m);

    if (is0()) {
        R = Q;
        return;
    }

    if (Q.is0()) {
        R = *this;
        return;
    }
    
    curve.m.sub (u, Q.x, x);
    bool ret = curve.m.inv (v, u);
    if (ret) {
        curve.m.sub (u, Q.y, y);
        curve.m.mul (lambda, u, v);
        curve.m.sqr (u, lambda);
        curve.m.sub (u, u, x);
        curve.m.sub (u, u, Q.x);    /* x3 = u = lambda^2 - P.x - Q.x */
        curve.m.sub (v, x, u);
        curve.m.mul (v, v, lambda);
        curve.m.sub (R.y, v, y);
        curve.m.set (R.x, u);
        R.finite = true;
    } else if (*this == Q) {
        dbl (R);
    } else {
        curve.m.set (R.x, u);
        R.finite = false;
    }
}

/* Computes R<-eP, with double-and-add algorithm.
 *
 * If the point at infinity is due to a failed inversion, the non-invertible
 * value is returned in R.x.
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::AffinePoint::smul (AffinePoint &R, const uint64_t e) const
{
    if (e == 0) {
        R.set0();
    } else if (e == 1) {
        R = *this;
    } else {
        uint64_t i;
        AffinePoint T(*this);

        i = UINT64_C(1) << 63;
        while ((i & e) == 0)
            i >>= 1;

        i >>= 1;

        for (; i > 0; i >>= 1) {
            T.dbl (T);
            if (e & i)
                T.add (T, *this);
        }

        R = T;
    }
}

/* Return the order of the point P on the Weierstrass curve defined by the curve
 * coefficient a modulo the prime m.
 * Looks for i in Hasse interval so that i*P = O, has complexity O(m^(1/4)).
 * If the group order is known to be == r (mod m), this can be supplied in
 * the variables "known_r" and" known_m".
 * XXX For now, this function only works for 64-bit arithmetic, so only the
 * 64-bit version is defined.
 */
template <typename MODULUS>
uint64_t ECWeierstrass<MODULUS>::AffinePoint::point_order (const uint64_t known_m,
    const uint64_t known_r, const int verbose) const
{
  AffinePoint Pi(curve), Pg(curve), Q(curve);
  Residue x(curve.m), d(curve.m);
  uint64_t min, max, i, j, order, cof, p;
  uint64_t giant_step, giant_min, baby_len;
  Integer tm;

  ASSERT (known_r < known_m);

  if (is0())
      return 1;
  
  curve.m.getmod (tm);

  if (verbose >= 2)
  {
    printf ("%s: ", __func__);
    curve.print(std::cout,  NULL);
    print(std::cout);
  }

  /* XXX Here we assume m fits in a uint64_t */
  i = (uint64_t) (2. * sqrt((double) (uint64_t)tm));
  min = (uint64_t) tm - i;
  max = (uint64_t) tm + i;

  /* Giant steps visit values == r (mod m), baby steps values == 0 (mod m) */
  giant_step = ceil(sqrt(2.*(double)i / (double) known_m));
  /* Round up to multiple of m */
  giant_step = ((giant_step - 1) / known_m + 1) * known_m;

  /* We test Pi +- Pj, where Pi = (giant_min + i*giant_step), i >= 0,
     and Pj = j*P, 0 <= j <= giant_step / 2.
     To ensure we can find all values >= min, ensure
     giant_min <= min + giant_step / 2.
     We also want giant_min == r (mod m) */
  giant_min = ((min + giant_step / 2) / known_m) * known_m + known_r;
  if (giant_min > min + giant_step / 2)
    giant_min -= known_m;
  if (verbose >= 2)
    printf ("known_m = %" PRIu64 ", known_r = %" PRIu64 ", giant_step = %" PRIu64 ", "
            "giant_min = %" PRIu64 "\n", known_m, known_r, giant_step, giant_min);

  baby_len = giant_step / known_m / 2 + 1;
  std::vector<AffinePoint> baby(baby_len, Pi);

  i = known_m;
  smul (Pg, i); /* Pg = m*P for now */
  if (Pg.is0())
    goto found_inf;

  if (1 < baby_len)
    baby[1] = Pg;

  if (2 < baby_len)
    {
      Pg.dbl (Pi);
      if (Pi.is0())
        {
          i = 2 * known_m;
          goto found_inf;
        }
      baby[2] = Pi;
    }

  for (i = 3; i < baby_len; i++)
    {
      Pi.add (Pi, Pg);
      if (Pi.is0())
        {
          i *= known_m;
          goto found_inf;
        }
      baby[i] = Pi;
    }

  /* Now compute the giant steps in [giant_min, giant_max] */
  i = giant_step;
  smul (Pg, i);
  if (Pg.is0())
    goto found_inf;

  i = giant_min;
  smul (Pi, i);
  if (Pi.is0())
    goto found_inf;

  while (i <= max + giant_step - 1)
    {
      /* Compare x-coordinate with stored baby steps. This makes it
         O(sqrt(p)) complexity, strictly speaking. */
      for (j = 1; j < baby_len; j++)
        if (curve.m.equal (Pi.x, baby[j].x))
          {
            if (curve.m.equal (Pi.y, baby[j].y))
              i -= j * known_m; /* Equal, so iP = jP and (i-j)P = 0 */
            else
              {
                curve.m.neg (Pi.y, Pi.y);
                if (curve.m.equal (Pi.y, baby[j].y))
                  i += j * known_m; /* Negatives, so iP = -jP and (i+j)P = 0 */
                else
                  {
                    std::cerr << "Matching x-coordinates, but y neither equal nor negatives" << std::endl;
                    std::cerr << "giant_min = " << giant_min << ",  i = " << i << ",  j = " << j << std::endl;
                    abort();
                  }
              }
            goto found_inf;
          }

      i += giant_step;
      Pi.add (Pi, Pg);
      if (Pi.is0())
        goto found_inf;
    }

  if (i > max)
  {
      fprintf (stderr, "ec_order: Error, no match found for p = %" PRIu64 ", "
               "min = %" PRIu64 ", max = %" PRIu64 ", giant_step = %" PRIu64 
               ", giant_min = %" PRIu64 "\n",
               (uint64_t) tm, min, max, giant_step, giant_min);
      abort ();
  }

found_inf:
  /* Check that i is a multiple of the order */
  smul (Pi, i);
  if (!Pi.is0())
    {
      Integer tx1, ty1;
      curve.m.get (tx1, x);
      curve.m.get (ty1, y);
#ifndef MODMPZ_MAXBITS
      fprintf (stderr, "ec_order: Error, %" PRIu64 "*(%" PRIu64 ", %" PRIu64 ") (mod %" PRIu64 ") is "
               "not the point at infinity\n",
               (uint64_t) i, (uint64_t) tx1, (uint64_t) ty1, (uint64_t) tm);
#endif
      return 0UL;
    }

  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  /* For each prime factor of the order, reduce the exponent of
     that prime factor as far as possible */

  cof = order = i;
  for (p = 2; p * p <= cof; p += 1 + p%2)
    if (cof % p == 0)
      {
        ASSERT (order % p == 0);
        /* Remove all factors of p */
        for (order /= p, cof /= p; order % p == 0; order /= p)
          {
            ASSERT(cof % p == 0);
            cof /= p;
          }
        ASSERT (cof % p != 0);

        /* Add factors of p again one by one, stopping when we hit
           point at infinity */
        smul (Pi, order);
        while (!Pi.is0()) {
            order *= p;
            Pi.smul (Pi, p);
        }
      }
  /* Now cof is 1 or a prime */
  if (cof > 1)
    {
      ASSERT (order % cof == 0);
      smul (Pi, order / cof);
      if (Pi.is0())
        order /= cof;
    }

  /* One last check that order divides real order */
  smul (Pi, order);
  if (!Pi.is0())
    {
      Integer tx1, ty1;
      curve.m.get (tx1, x);
      curve.m.get (ty1, y);
      fprintf (stderr, "ec_order: Error, final order %" PRIu64 " is wrong\n",
               order);
      abort ();
    }

  return order;
}

/* Computes R=2P, with ? muls (? muls and ? squares) and ? add/sub.
 *
 * It is permissible to let P and Q use the same memory.
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::ProjectivePoint::dbl (ProjectivePoint &R) const
{
    ASSERT_EXPENSIVE(curve.is_same(R.curve));

    /* If input is the point at infinity or a point of order 2,  then result
     * is point at infinity. */
    if (is0() || curve.m.is0(y)) {
        R.set0();
        return;
    }

    Residue xx(curve.m), zz(curve.m), w(curve.m), u(curve.m), s(curve.m),
            ss(curve.m), sss(curve.m), r(curve.m), rr(curve.m), B(curve.m),
            t(curve.m), h(curve.m);
    /* TODO reduce number of var */

    curve.m.sqr (xx, x);
    curve.m.sqr (zz, z);

    curve.m.mul (w, curve.a, zz);
    curve.m.add (t, xx, xx);
    curve.m.add (t, t, xx);
    curve.m.add (w, w, t);

    curve.m.mul (s, y, z);
    curve.m.add (s, s, s);
    curve.m.sqr (ss, s);
    curve.m.mul (sss, ss, s);
    curve.m.mul (r, y, s);
    curve.m.sqr (rr, r);

    curve.m.add (B, x, r);
    curve.m.sqr (B, B);
    curve.m.sub (B, B, xx);
    curve.m.sub (B, B, rr);
    curve.m.sqr (h, w);
    curve.m.sub (h, h, B);
    curve.m.sub (h, h, B);

    curve.m.mul (R.x, h, s);
    curve.m.sub (R.y, B, h);
    curve.m.mul (R.y, R.y, w);
    curve.m.sub (R.y, R.y, rr);
    curve.m.sub (R.y, R.y, rr);
    curve.m.set (R.z, sss);
}

/* Computes R=P+Q, with 14 muls (12 muls and 2 squares) and 7 add/sub.
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::ProjectivePoint::add (ProjectivePoint &R, const ProjectivePoint &Q) const
{
    ASSERT_EXPENSIVE(curve.is_same(R.curve));
    ASSERT_EXPENSIVE(curve.is_same(Q.curve));
    Residue t1(curve.m), t2(curve.m), t3(curve.m), u(curve.m), uu(curve.m),
              v(curve.m), vv(curve.m), vvv(curve.m), r(curve.m), A(curve.m);
    /* TODO reduce number of var */

    if (is0()) {
        R = Q;
        return;
    }

    if (Q.is0()) {
        R = *this;
        return;
    }

    if (*this ==  Q) {
        dbl(R);
        return;
    }

    if (*this == -Q) {
        R.set0();
        return;
    }

    curve.m.mul (t1, x, Q.z);
    curve.m.mul (t2, y, Q.z);
    curve.m.mul (t3, z, Q.z);

    curve.m.mul (u, Q.y, z);
    curve.m.sub (u, u, t2);
    curve.m.sqr (uu, u);

    curve.m.mul (v, Q.x, z);
    curve.m.sub (v, v, t1);
    curve.m.sqr (vv, v);
    curve.m.mul (vvv, vv, v);

    curve.m.mul (r, vv, t1);

    curve.m.mul (A, uu, t3);
    curve.m.sub (A, A, vvv);
    curve.m.sub (A, A, r);
    curve.m.sub (A, A, r);

    curve.m.mul (R.x, v, A);

    curve.m.sub (R.y, r, A);
    curve.m.mul (R.y, R.y, u);
    curve.m.mul (t2, t2, vvv);
    curve.m.sub (R.y, R.y, t2);

    curve.m.mul (R.z, vvv, t3);
}

/* Computes P<-eP, with double-and-add algorithm.
 *    - m : modulus
 *    - a : curve coefficient
 */
template <typename MODULUS>
void ECWeierstrass<MODULUS>::ProjectivePoint::smul (ProjectivePoint &R, const uint64_t e) const
{
    if (e == 0) {
        R.set0();
    } else if (e == 1) {
        R = *this;
    } else {
        uint64_t i;
        ProjectivePoint T(*this);

        i = UINT64_C(1) << 63;
        while ((i & e) == 0)
            i >>= 1;

        i >>= 1;

        for (; i > 0; i >>= 1) {
            T.dbl (T);
            if (e & i)
                T.add (T, *this);
        }

        R = T;
    }
}

template class ECWeierstrass<Modulus64>;
template class ECWeierstrass<ModulusREDC64>;
template class ECWeierstrass<ModulusREDC126>;
template class ECWeierstrass<ModulusMPZ>;
