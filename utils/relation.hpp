#ifndef RELATION_HPP_
#define RELATION_HPP_

#ifndef __cplusplus
#error "This is C++-only"
#endif

#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "cado_poly.h"
#include "gmp_aux.h"
#include "cxx_mpz.hpp"

/* This is a "slow" relation type. We can afford having vectors and mpzs
 * and such.
 */

struct relation_ab : public std::vector<cxx_mpz> {
    typedef std::vector<cxx_mpz> super;
    operator bool() const { return !empty(); }
    relation_ab() {}
    /* Initializing with a,b means a+bx.
     * Initializing with a vector a,b,c,d means a+bx+c*x^2+d*x^3
     */
    relation_ab(uint64_t a, int64_t b) {
        assign(2,cxx_mpz());
        mpz_set_uint64((*this)[0], a);
        mpz_set_int64((*this)[1], b);
    }
    relation_ab(mpz_srcptr az, mpz_srcptr bz) {
        assign(2,cxx_mpz());
        mpz_set((*this)[0], az);
        mpz_set((*this)[1], bz);
    }
    relation_ab(super const& o): super(o) {}
    bool operator<(const relation_ab& o) const {
        return lexicographical_compare(begin(), end(), o.begin(), o.end());
    }
};

struct relation : public relation_ab {
    relation_ab& ab() { return *this; }
    relation_ab const & ab() const { return *this; }
    /* Note that for the rational side, we do not compute r !!! */
    struct pr {
        cxx_mpz p,r;
        int e;
        pr() {
            e=0;
        }
        pr(mpz_srcptr ap, mpz_srcptr ar, int ae=1) {
            mpz_set(p, ap);
            if (ar)
                mpz_set(r, ar);
            else
                mpz_set_ui(r, 0);
            e = ae;
        }
        pr(unsigned long ap, unsigned long ar, int ae=1) {
            mpz_set_ui(p, ap);
            mpz_set_ui(r, ar);
            e = ae;
        }
    };
    int rational_side;    /* index of the rational side, if any */
    int nb_polys;         /* number of polynoms, default = 2 */
    std::vector<pr> sides[NB_POLYS_MAX]; /* pr's are stored w.r.t. side */

    relation() {}
    operator bool() const { return (bool) (relation_ab) *this; }
    relation(int64_t a, uint64_t b, int rational_side = -1, int nb_polys = 2)
        : relation_ab(a,b), rational_side(rational_side), nb_polys(nb_polys)
    {}
    relation(mpz_t a, mpz_t b, int rational_side = -1, int nb_polys = 2)
        : relation_ab(a,b), rational_side(rational_side), nb_polys(nb_polys)
    {}

    void add(int side, mpz_srcptr p, mpz_srcptr r, int ae=1) {
        sides[side].push_back(pr(p, r, ae));
    }
    void add(int side, unsigned long p, unsigned long r, int ae=1) {
        sides[side].push_back(pr(p, r, ae));
    }

    /* the single-member add() functions recompute r for the algebraic
     * side, based on a,b
     */
    void addv(int side, mpz_srcptr p, int ae);
    void addv(int side, unsigned long p, int ae);
    void add(int side, mpz_srcptr p) { addv(side, p, 1); }
    void add(int side, unsigned long p) { addv(side, p, 1); }

    void print(FILE * f, const char * prefix) const;
    int parse(const char * line);


    void fixup_r(bool also_rational = false);
    void compress();
};

struct pr_cmp {
    bool operator()(relation::pr const& a, relation::pr const& b) const {
        int c = mpz_cmp(a.p, b.p);
        if (c) { return c < 0; }
        c = mpz_cmp(a.r, b.r);
        return c < 0;
    }
};

#endif	/* RELATION_HPP_ */
