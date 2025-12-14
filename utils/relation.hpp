#ifndef CADO_RELATION_HPP
#define CADO_RELATION_HPP

#include <cstdint>
#include <cstdio>

#include <array>
#include <istream> // std::istream // IWYU pragma: keep
#include <ostream> // std::ostream // IWYU pragma: keep
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>      // mpz_srcptr
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "gmp_aux.h"
#include "cxx_mpz.hpp"

struct relation_ab {
    int64_t a = 0;		/* only a is allowed to be negative */
    uint64_t b = 0;
    cxx_mpz az = 0;
    cxx_mpz bz = 0;
    /* We need to have this, because a pair (a,b) does not uniquely
     * identify a relation: there can be several side pairs which lead to
     * a relation.
     */
    std::array<int, 2> active_sides { 0, 1 };
    explicit operator bool() const { return a || b; }
    relation_ab() = default;
    relation_ab(int64_t a, uint64_t b)
        : a(a)
        , b(b)
        , az(a)
        , bz(b)
    {
    }
    relation_ab(mpz_srcptr _az, mpz_srcptr _bz)
        : a(mpz_get_int64(_az))
        , b(mpz_get_uint64(_bz))
        , az(_az)
        , bz(_bz)
    {
    }
    auto operator<=>(const relation_ab& o) const {
        using T = std::tuple<cxx_mpz const &, cxx_mpz const &, int, int>;
        T const me   { az,   bz,   active_sides[0],   active_sides[1] };
        T const them { o.az, o.bz, o.active_sides[0], o.active_sides[1] };
        return me <=> them;
    }
    bool operator==(const relation_ab& o) const {
        return operator<=>(o) == 0;
    }
    friend std::istream& operator>>(std::istream&, relation_ab&);
    friend std::ostream& operator<<(std::ostream&, relation_ab const &);
};

struct relation : public relation_ab {
    /* Note that for the rational side, we do not compute r !!! */
    struct pr {
        cxx_mpz p = 0,r = 0;
        int e = 0;
        pr(cxx_mpz ap, cxx_mpz ar, int ae=1)
            : p(std::move(ap))
            , r(std::move(ar))
            , e(ae)
        {
        }
        pr(unsigned long ap, unsigned long ar, int ae=1)
            : p(ap)
            , r(ar)
            , e(ae)
        {
        }
        pr() = default;
        auto operator<=>(pr const & b) const {
            if (auto c = p <=> b.p ; c != 0) return c;
            return r <=> b.r;
        }
        bool operator==(pr const & b) const {
            return operator<=>(b) == 0;
        }
    };
    int rational_side = -1;   /* index of the rational side, if any */
    std::array<std::vector<pr>, 2> sides; /* pr's are stored w.r.t. side */

    relation() = default;
    relation_ab const & ab() const { return *this; }
    operator bool() const { return bool(ab()); }
    relation(int64_t a, uint64_t b, int rational_side = -1)
        : relation_ab(a,b)
        , rational_side(rational_side)
    {}
    relation(mpz_srcptr a, mpz_srcptr b, int rational_side = -1)
        : relation_ab(a,b)
        , rational_side(rational_side)
    {}

    void add(unsigned int side_index, cxx_mpz const & p, cxx_mpz const & r) {
        sides[side_index].emplace_back(p, r);
    }
    void add(unsigned int side_index, unsigned long p, unsigned long r) {
        sides[side_index].emplace_back(p, r);
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

extern std::istream& operator>>(std::istream&, relation_ab&);
extern std::ostream& operator<<(std::ostream&, relation_ab const &);

namespace fmt {
    template <> struct formatter<relation>: ostream_formatter {};
    template <> struct formatter<relation_ab>: ostream_formatter {};
}

#endif	/* CADO_RELATION_HPP */
