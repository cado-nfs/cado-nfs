#ifndef RELATION_HPP_
#define RELATION_HPP_

#ifndef __cplusplus
#error "This is C++-only"
#endif

#include <algorithm> // max
#include <cstdint>
#include <cstdio>
#include <vector>
#include <array>
#include <tuple>
#include <istream> // std::istream // IWYU pragma: keep
#include <ostream> // std::ostream // IWYU pragma: keep
#include <gmp.h>      // mpz_srcptr

#include "cado_poly.h"
#include "gmp_aux.h"
#include "cxx_mpz.hpp"

struct relation_ab {
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
    cxx_mpz az;
    cxx_mpz bz;
    /* We need to have this, because a pair (a,b) does not uniquely
     * identify a relation: there can be several side pairs which lead to
     * a relation.
     */
    std::array<int, 2> active_sides;
    operator bool() const { return a || b; }
    relation_ab() {
        a=0;
        b=0;
        active_sides[0] = 0;
        active_sides[1] = 1;
    }
    relation_ab(int64_t a, uint64_t b) : a(a), b(b) {
        mpz_set_int64(az, a);
        mpz_set_uint64(bz, b);
        active_sides[0] = 0;
        active_sides[1] = 1;
    }
    relation_ab(mpz_srcptr _az, mpz_srcptr _bz) {
        mpz_set(az, _az);
        mpz_set(bz, _bz);
        a = mpz_get_int64(az);
        b = mpz_get_uint64(bz);
        active_sides[0] = 0;
        active_sides[1] = 1;
    }
    bool operator<(const relation_ab& o) const {
        typedef std::tuple<cxx_mpz const &, cxx_mpz const &, int, int> T;
        T me { az, bz, active_sides[0], active_sides[1] };
        T them { o.az, o.bz, o.active_sides[0], o.active_sides[1] };
        return me < them;
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
    int rational_side = -1;   /* index of the rational side, if any */
    std::array<std::vector<pr>, 2> sides; /* pr's are stored w.r.t. side */

    relation() {}
    operator bool() const { return (bool) (relation_ab) *this; }
    relation(int64_t a, uint64_t b, int rational_side = -1)
        : relation_ab(a,b)
        , rational_side(rational_side)
    {}
    relation(mpz_srcptr a, mpz_srcptr b, int rational_side = -1)
        : relation_ab(a,b)
        , rational_side(rational_side)
    {}

    void add(unsigned int side_index, mpz_srcptr p, mpz_srcptr r) {
        sides[side_index].push_back(pr(p, r));
    }
    void add(unsigned int side_index, unsigned long p, unsigned long r) {
        sides[side_index].push_back(pr(p, r));
    }

    /* the single-member add() functions recompute r for the algebraic
     * side, based on a,b
     */
    void add(unsigned int side_index, mpz_srcptr p);
    void add(unsigned int side_index, unsigned long p);

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

#endif	/* RELATION_HPP_ */
