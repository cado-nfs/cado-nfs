#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> /* for PRId64 */
#include <ctype.h> /* for isxdigit */
#include <string.h>
#include <errno.h>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "relation.hpp"
#include "gzip.h"
#include "timing.h"
#include "portability.h"
#include "relation-tools.h"
#include "relation.hpp"

using namespace std;

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
    int side = 0;

    if (gmp_sscanf(line, "%Zd,%Zd:%n", (mpz_ptr) az, (mpz_ptr)bz, &consumed) < 2)
        return 0;
    a = mpz_get_int64(az);
    b = mpz_get_uint64(bz);

    for(int i = 0; i < NB_POLYS_MAX; i++)
	sides[i].clear();

    if (line[consumed] == ':') {
        side++;
        consumed++;
    }

    while(line[consumed] != '\0' && line[consumed] != '\n') {
        unsigned long p;
        int consumed_p;
        if (sscanf(line + consumed, "%lx%n", &p, &consumed_p) < 1)
            return 0;
	// take care to the "::" problem in MNFS
	// printf("CONSUMED: %lu %d\n", p, consumed_p);
        add(side, p);
        consumed += consumed_p;
        if (line[consumed] == ',')
            consumed++;
        else if (line[consumed] == ':') {
            side++;
            ASSERT_ALWAYS(side < NB_POLYS_MAX);
            consumed++;
        }
    }
    nb_polys = side+1;
    compress();
    fixup_r();
    return 1;
}

std::istream& operator>>(std::istream& is, relation& rel)
{
    std::string s;
    if (!getline(is, s, '\n') || !rel.parse(s.c_str())) {
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
    os << rel.az << ',' << rel.bz;
    os << std::hex;
    for(int side = 0 ; side < rel.nb_polys ; side++) {
        os << ':';
        bool comma=false;
        for(auto const& v : rel.sides[side]) {
            for(int e = v.e ; e ; e--) {
                if (comma) os << ',';
                os << v.p;
                comma = true;
            }
        }
    }
    return os;
}

void relation::add(int side, mpz_srcptr p)
{
    /* we have to compute a/b mod p. Since we're not critical here, don't
     * bother.
     */
    if (side == rational_side) {
        add(side, p, 0);
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

        sides[side].push_back(x);
    }
}

void relation::add(int side, unsigned long p)
{
    if (side == rational_side) {
        add(side, p, 0);
    } else {
        /* use the function provided in relation-tools.c */
        add(side, p, relation_compute_r(a, b, p));
    }
}

void relation::fixup_r(bool also_rational)
{
    for(int side = 0 ; side < nb_polys ; side++) {
        if (!also_rational && side == rational_side)
            continue;
        for(unsigned int i = 0 ; i < sides[side].size() ; i++) {
            if (mpz_cmp_ui(sides[side][i].r,0) == 0) {
                pr & x(sides[side][i]);

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

inline bool operator==(relation::pr const& a, relation::pr const& b) {
    return mpz_cmp(a.p, b.p) == 0 && mpz_cmp(a.r, b.r) == 0;
}


void relation::compress()
{
    for(int side = 0 ; side < nb_polys ; side++) {
        vector<pr> & v(sides[side]);
        sort(v.begin(), v.end(), pr_cmp());
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

