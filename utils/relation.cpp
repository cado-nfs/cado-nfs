#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <istream>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <ios>
#include <locale>
#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "macros.h"
#include "relation.hpp"
#include "relation-tools.h"
#include "misc.h"
#include "utils_cxx.hpp"
/*
 * Convention for I/O of rels:
 *   a and b are printed in decimal
 *   primes are printed in hexadecimal.
 *
 */

std::istream& operator>>(std::istream& is, relation_ab& rel)
{
    /* Note that this assumes that the istream is in the proper shape */
    if (is >> rel.az >> expect(",") >> rel.bz) {
        rel.a = mpz_get_int64(rel.az);
        rel.b = mpz_get_uint64(rel.bz);
    }
    return is;
}

int
relation::parse(const char *line)
{
    if (line[0] == '#')
        return 0;

    /* This is just to let the libc++ parsing happily grok the integers
     * in relations as they are, and _then_ we implement our
     * syntax verifications
     */
    struct relation_locale: std::ctype<char>
    {
        relation_locale(): std::ctype<char>(get_table()) {}

        static std::ctype_base::mask const* get_table()
        {
            static std::vector<std::ctype_base::mask>
                rc(std::ctype<char>::table_size,std::ctype_base::space);
            std::fill(&rc['0'], &rc['9'+1], std::ctype_base::digit);
            std::fill(&rc['a'], &rc['f'+1], std::ctype_base::xdigit);
            std::fill(&rc['A'], &rc['F'+1], std::ctype_base::xdigit);
            rc['-'] = std::ctype_base::punct;
            return rc.data();
        }
    };

    const std::string S(line);
    std::istringstream is(S);
    is.imbue(std::locale(std::locale(), new relation_locale()));

    if (!(is >> static_cast<relation_ab &>(*this)))
        return 0;

    if (is.peek() == '@') {
        is >> expect("@") >> active_sides[0] >> expect(",") >> active_sides[1];
    } else {
        active_sides[0] = 0;
        active_sides[1] = 1;
    }

    is >> expect(":");

    sides[0].clear();
    sides[1].clear();

    bool comma_allowed = false;
    unsigned int side_index = 0;

    is >> std::hex;

    for( ; !is.eof() ; ) {
        auto const c = is.peek();
        if (c == ':') {
            is.get();
            side_index++;
            /* We do not support specifying sides in a relation by typing
             * in zillions of colons.  That is not well-defined, since we
             * might encounter situations where a-b*alpha is a unit in
             * one of the number fields (see #21707)
             */
            ASSERT_ALWAYS(side_index < 2);
            comma_allowed = false;
            continue;
        } else if (comma_allowed && c == ',') {
            is.get();
        } else if (c == '\n') {
            break;
        }

        unsigned long p;
        is >> p;
        if (!is)
            return 0;
        add(side_index, p);
        comma_allowed = true;
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
    int const rc = fputs(os.str().c_str(), file);
    if (rc < 0) {
        perror("Error writing relation");
        abort();
    }
}

std::ostream& operator<<(std::ostream& os, relation_ab const &rel)
{
    IoStreamFlagsRestorer const dummy(os);
    os << std::dec << rel.az << ',' << rel.bz;
    if (rel.active_sides[0] != 0 || rel.active_sides[1] != 1)
        os << "@" << rel.active_sides[0] << "," << rel.active_sides[1];
    return os;
}

std::ostream& operator<<(std::ostream& os, relation const &rel)
{
    IoStreamFlagsRestorer const dummy(os);
    {
        os << static_cast<relation_ab const &>(rel);
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
    if (active_sides[side_index] == rational_side) {
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
    if (active_sides[side_index] == rational_side) {
        add(side_index, p, 0);
    } else {
        /* use the function provided in relation-tools.c */
        add(side_index, p, relation_compute_r(a, b, p));
    }
}

void relation::fixup_r(bool also_rational)
{
    for(unsigned int side_index = 0 ; side_index < sides.size() ; side_index++) {
        int const side = active_sides[side_index];
        if (!also_rational) {
            if (side == rational_side)
                continue;
        }
        for(auto & x : sides[side_index]) {
            if (mpz_cmp_ui(x.r,0) == 0) {
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

void relation::compress()
{
    for(auto & v : sides) {
        std::ranges::sort(v);
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
