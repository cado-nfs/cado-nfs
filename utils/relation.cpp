#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> /* for PRId64 */
#include <ctype.h> /* for isxdigit */
#include <string.h>
#include <errno.h>
#include <algorithm>

#include "relation.hpp"
#include "gzip.h"
#include "timing.h"
#include "portability.h"
#include "relation-tools.h"
#include "relation.hpp"

using namespace std;

/*
 * Convention for I/O of rels:
 *   "X " is prepended before each relation as an indicator for the new
 *   format.
 *   a and b are printed in hexadecimal, decimal is legacy now.
 *   primes are printed in hexadecimal.
 *
 */

int
relation::parse(const char *line)
{
    int consumed;
    int side = 0;

    if (line[0] == '-' || isdigit(line[0])) {
        /* legacy format */
        cxx_mpz a,b;
        if (gmp_sscanf(line, "%Zd,%Zd:%n", a, b, &consumed) < 2)
            return 0;
        /* With the legacy format, (a,b) encodes a-bx. Given that we
         * still have the normalization b>0, we'll put -a+bx instead.
         */
        mpz_neg(a,a);
        ab() = relation_ab(a, b);
    } else if (line[0] == 'X') {
        ab().clear();
        ab().reserve(2);
        consumed = 1;
        for(int c ; line[consumed] == ',' || line[consumed] == ' ' ; consumed += c) {
            consumed++;
            cxx_mpz x;
            if (gmp_sscanf(line + consumed, "%Zx%n", x, &c) < 1)
                return 0;
            ab().push_back(x);
        }
        ASSERT_ALWAYS(line[consumed++]==':');
    } else {
        return 0;
    }

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

void
relation::print (FILE *file, const char *prefix) const
{
    char buf[RELATION_MAX_BYTES], *p = buf;
    char * fence = buf + sizeof(buf);
    int c;

    c = strlcpy(p, prefix, fence - p);
    p += strnlen(prefix, sizeof(buf));

    c = gmp_snprintf(p, fence - p, "X");
    p += c;

    for(unsigned int i = 0 ; i < ab().size() ; i++) {
        c = gmp_snprintf(p, fence - p, "%c%Zx", i ? ',' : ' ', ab()[i]);
        p += c;
    }

    for(int side = 0 ; side < nb_polys ; side++) {
        if (p + 1 < fence) *p++ = ':';
        char * op = p;
        for(unsigned int i = 0 ; i < sides[side].size() ; i++) {
            for(int e = sides[side][i].e ; e ; e--) {
                c = gmp_snprintf(p, fence - p, "%Zx,", sides[side][i].p);
                p += c;
            }
        }
        if (p > op) p--;
    }
    if (p + 1 < fence) *p++ = '\n';
    *p = '\0';
    /* print all in one go, this gives us a chance to minimize the risk
     * of mixing up output if we don't pay attention to I/O locking too
     * much in multithreaded context. */
    size_t written = fwrite (buf, 1, p - buf, file);
    if (written != (size_t) (p - buf)) {
        perror("Error writing relation");
        abort();
    }
}

void relation::addv(int side, mpz_srcptr p, int ae)
{
    /* we have to compute a/b mod p. Since we're not critical here, don't
     * bother.
     */
    if (side == rational_side) {
        add(side, p, 0, ae);
    } else if (ab().size() == 2) {
        pr x;
        mpz_set(x.p, p);

        mpz_neg(x.r, ab()[1]);
        if (mpz_invert(x.r, x.r, x.p)) {
            mpz_mul(x.r, x.r, ab()[0]);
            mpz_mod(x.r, x.r, x.p);
        } else {
            mpz_set(x.r, x.p);
        }
        x.e = ae;

        sides[side].push_back(x);
    } else {
        abort();        /* implement me ! */
        add(side, p, 0, ae);
    }
}

void relation::addv(int side, unsigned long p, int ae)
{
    if (side == rational_side) {
        add(side, p, 0, ae);
    } else if (ab().size() == 2) {
        /* use the function provided in relation-tools.c */
        int64_t a = mpz_get_int64(ab()[0]);
        int64_t b = mpz_get_int64(ab()[1]);
        /* the function relation_compute_r takes a uint64_t so we want to
         * be sure */
        ASSERT_ALWAYS(a >= 0);
        add(side, p, relation_compute_r(a, b, p), ae);
    } else {
        abort();        /* implement me ! */
        add(side, p, 0, ae);
    }
}

void relation::fixup_r(bool also_rational)
{
    /* r does not make sense if we're not dealing with a,b in the first
     * place...
     */
    if (ab().size() != 2)
        return;
    for(int side = 0 ; side < nb_polys ; side++) {
        if (!also_rational && side == rational_side)
            continue;
        for(unsigned int i = 0 ; i < sides[side].size() ; i++) {
            if (mpz_cmp_ui(sides[side][i].r,0) == 0) {
                pr & x(sides[side][i]);

                /* compute r=-a/b mod p */
                mpz_neg(x.r, ab()[1]);
                if (mpz_invert(x.r, x.r, x.p)) {
                    mpz_mul(x.r, x.r, ab()[0]);
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
            if (v[j].e == 0)
                j--;
        }
        v.erase(v.begin() + j, v.end());
    }
}

