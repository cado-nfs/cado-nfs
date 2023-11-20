#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <memory>
// IWYU asks for <memory> because of allocator_traits<>::value_type ; WTH ?
#include <cstdio>
#include <cstdlib>
#include <istream> // std::istream
#include <ostream> // std::ostream
#include <algorithm>
#include <sstream> // std::ostringstream // IWYU pragma: keep
#include <string>
#include <iomanip> // std::hex // IWYU pragma: keep
#include <gmp.h>

#include "macros.h" /* for ASSERT_ALWAYS */
#include "relation.hpp"
#include "relation-tools.h"
#include "misc.h"

/*
 * Convention for I/O of rels:
 *   a and b are printed in decimal
 *   primes are printed in hexadecimal.
 *
 */

int
relation::parse(const char *line)
{
    int consumed;

    if (gmp_sscanf(line, "%Zd,%Zd%n", (mpz_ptr) az, (mpz_ptr)bz, &consumed) < 2)
        return 0;
    a = mpz_get_int64(az);
    b = mpz_get_uint64(bz);

    if (line[consumed] == '@') {
        consumed++;
        int c;
        if (sscanf(line + consumed, "%u,%u%n",
                &active_sides[0],
                &active_sides[1],
                &c) < 2)
            return 0;
        consumed += c;
    } else {
        active_sides[0] = 0;
        active_sides[1] = 1;
    }
    ASSERT_ALWAYS (line[consumed] == ':');
    consumed++;

    sides[0].clear();
    sides[1].clear();

    int comma_allowed = 0;

    unsigned int side_index = 0;
    while(line[consumed] != '\0' && line[consumed] != '\n') {
        if (line[consumed] == ':') {
            side_index++;
            /* We do not support specifying sides in a relation by typing
             * in zillions of colons.  That is not well-defined, since we
             * might encounter situations where a-b*alpha is a unit in
             * one of the number fields (see #21707)
             */
            ASSERT_ALWAYS(side_index < 2);
            consumed++;
            continue;
        } else if (comma_allowed && line[consumed] == ',') {
            consumed++;
            comma_allowed = 0;
        }

        unsigned long p;
        int consumed_p;
        if (sscanf(line + consumed, "%lx%n", &p, &consumed_p) < 1)
            return 0;
        add(side_index, p);
        consumed += consumed_p;
        comma_allowed = 1;
    }
    compress();
    fixup_r();
    return 1;
}

std::istream& operator>>(std::istream& is, relation& rel)
{
    std::string s;
    if (!std::getline(is, s, '\n') || !rel.parse(s.c_str())) {
        is.setstate(std::ios_base::failbit);
        rel = relation();
    }
    return is;
}

void relation::print (FILE *file, const char *prefix) const
{
    std::ostringstream os;
    if (prefix) os << prefix;
    os << *this << '\n';
    int rc = fputs(os.str().c_str(), file);
    if (rc < 0) {
        perror("Error writing relation");
        abort();
    }
}

std::ostream& operator<<(std::ostream& os, relation const &rel)
{
    IoStreamFlagsRestorer dummy(os);
    {
        os << rel.az << ',' << rel.bz;
        os << std::hex;
        for(auto const & s : rel.sides) {
            os << ':';
            bool comma=false;
            for(auto const& v : s) {
                for(int e = v.e ; e ; e--) {
                    if (comma) os << ',';
                    os << v.p;
                    comma = true;
                }
            }
        }
    }
    return os;
}

void relation::add(unsigned int side_index, mpz_srcptr p)
{
    /* we have to compute a/b mod p. Since we're not critical here, don't
     * bother.
     */
    if ((int) active_sides[side_index] == rational_side) {
        add(side_index, p, 0);
    } else {
        pr x;
        mpz_set(x.p, p);

        mpz_set(x.r, bz);
        if (mpz_invert(x.r, x.r, x.p)) {
            mpz_mul(x.r, x.r, az);
            mpz_mod(x.r, x.r, x.p);
        } else {
            mpz_set(x.r, x.p);
        }

        sides[side_index].push_back(x);
    }
}

void relation::add(unsigned int side_index, unsigned long p)
{
    if ((int) active_sides[side_index] == rational_side) {
        add(side_index, p, 0);
    } else {
        /* use the function provided in relation-tools.c */
        add(side_index, p, relation_compute_r(a, b, p));
    }
}

void relation::fixup_r(bool also_rational)
{
    for(unsigned int side_index = 0 ; side_index < sides.size() ; side_index++) {
        int side = active_sides[side_index];
        if (!also_rational) {
            if ((int) side == rational_side)
                continue;
        }
        for(unsigned int i = 0 ; i < sides[side_index].size() ; i++) {
            if (mpz_cmp_ui(sides[side_index][i].r,0) == 0) {
                pr & x(sides[side_index][i]);

                mpz_set(x.r, bz);
                if (mpz_invert(x.r, x.r, x.p)) {
                    mpz_mul(x.r, x.r, az);
                    mpz_mod(x.r, x.r, x.p);
                } else {
                    mpz_set(x.r, x.p);
                }
            }
        }
    }
}

static inline bool operator==(relation::pr const& a, relation::pr const& b) {
    return mpz_cmp(a.p, b.p) == 0 && mpz_cmp(a.r, b.r) == 0;
}


void relation::compress()
{
    for(auto & v : sides) {
        std::sort(v.begin(), v.end(), pr_cmp());
        unsigned int j = 0;
        for(unsigned int i = 0; i < v.size() ; j++) {
            if (j < i) {
                v[j] = v[i];
            }
            for(i++ ; i < v.size() && v[i] == v[j] ; i++) {
                v[j].e += v[i].e;
            }
        }
        v.erase(v.begin() + j, v.end());
    }
}
