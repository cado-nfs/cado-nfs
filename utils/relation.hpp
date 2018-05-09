#ifndef RELATION_H_
#define RELATION_H_

#ifndef __cplusplus
#error "This is C++-only"
#endif

#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <istream>
#include <ostream>

#include "cado_poly.h"
#include "gmp_aux.h"
#include "cxx_mpz.hpp"

struct relation_ab {
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
    cxx_mpz az;
    cxx_mpz bz;
    operator bool() const { return a || b; }
    relation_ab() {
        a=0;
        b=0;
    }
    relation_ab(int64_t a, uint64_t b) : a(a), b(b) {
        mpz_set_int64(az, a);
        mpz_set_uint64(bz, b);
    }
    relation_ab(mpz_srcptr _az, mpz_srcptr _bz) {
        mpz_set(az, _az);
        mpz_set(bz, _bz);
        a = mpz_get_int64(az);
        b = mpz_get_uint64(bz);
    }
    bool operator<(const relation_ab& o) const {
        return ((mpz_cmp(az, o.az) < 0) ||
                (mpz_cmp(az, o.az) == 0 &&
                 mpz_cmp(bz, o.bz) < 0));
    }
};

struct relation : public relation_ab {
    /* Note that for the rational side, we do not compute r !!! */
    struct pr {
        cxx_mpz p,r;
        int e = 0;
        pr() = default;
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
        pr(pr const&) = default;
        pr(pr &&) = default;
        pr& operator=(const pr&) = default;
    };
    int rational_side = -1;    /* index of the rational side, if any */
    int nb_polys = 0;         /* number of polynoms, default = 2 */
    std::vector<pr> sides[NB_POLYS_MAX]; /* pr's are stored w.r.t. side */

    relation() {}
    operator bool() const { return (bool) (relation_ab) *this; }
    relation(int64_t a, uint64_t b, int rational_side = -1, int nb_polys = 2)
        : relation_ab(a,b), rational_side(rational_side), nb_polys(nb_polys)
    {}
    relation(mpz_srcptr a, mpz_srcptr b, int rational_side = -1, int nb_polys = 2)
        : relation_ab(a,b), rational_side(rational_side), nb_polys(nb_polys)
    {}

    void add(int side, mpz_srcptr p, mpz_srcptr r) {
        sides[side].push_back(pr(p, r));
    }
    void add(int side, unsigned long p, unsigned long r) {
        sides[side].push_back(pr(p, r));
    }

    /* the single-member add() functions recompute r for the algebraic
     * side, based on a,b
     */
    void add(int side, mpz_srcptr p);
    void add(int side, unsigned long p);

    void print(FILE * f, const char * prefix) const;
    int parse(const char * line);

    void fixup_r(bool also_rational = false);
    void compress();
    friend std::istream& operator>>(std::istream&, relation&);
    friend std::ostream& operator<<(std::ostream&, relation const &);
};

extern std::istream& operator>>(std::istream&, relation&);
extern std::ostream& operator<<(std::ostream&, relation const &);

struct pr_cmp {
    bool operator()(relation::pr const& a, relation::pr const& b) const {
        int c = mpz_cmp(a.p, b.p);
        if (c) { return c < 0; }
        c = mpz_cmp(a.r, b.r);
        return c < 0;
    }
};

#endif	/* RELATION_H_ */
