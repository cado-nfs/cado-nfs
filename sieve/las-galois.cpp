#include "cado.h" // IWYU pragma: keep

#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <vector>
#include <ostream>

#include <gmp.h>

#include "runtime_numeric_cast.hpp"
#include "galois_action.hpp"
#include "cxx_mpz.hpp"
#include "las-galois.hpp"
#include "relation.hpp"
#include "macros.h"
#include "verbose.h"

static void adwg(std::ostream& os, const char *comment, unsigned long *cpt,
		 relation &rel, int64_t a, int64_t b)
{
    if (b < 0) { a = -a; b = -b; }
    rel.a = a;
    rel.b = runtime_numeric_cast<uint64_t>(b);
    if (comment) os << comment;
    os << rel << '\n';
    *cpt += (*comment != '\0');
}

/* removing p^vp from the list of factors in rel. */
static void remove_galois_factors(relation &rel, int p, int vp)
{
    bool ok = false;

    for(auto & S : rel.sides) {
        for(auto & pe : S) {
	    if(mpz_cmp_ui(pe.p, p) == 0) {
                /* We should not have p several times, I guess */
                ASSERT_ALWAYS(!ok);
		ok = true;
		ASSERT_ALWAYS(pe.e >= vp);
		pe.e -= vp;
	    }
        }
    }
    /* indeed, p was present */
    ASSERT_ALWAYS(ok == 1);
}

/* adding p^vp to the list of factors in rel. */
static void add_galois_factors(relation &rel, int p, int vp)
{
    int ok[2] = {0, 0};

    for(unsigned int side = 0 ; side < rel.sides.size() ; side++) {
	for(unsigned int i = 0 ; i < rel.sides[side].size(); i++)
	    if(mpz_cmp_ui(rel.sides[side][i].p, p) == 0) {
		ok[side] = 1;
		rel.sides[side][i].e += vp;
	    }
    }
    // FIXME: are we sure this is safe?
    for(unsigned int side = 0 ; side < rel.sides.size() ; side++)
	if(ok[side] == 0)
	    /* we must add p^vp */
	    for(int i = 0; i < vp; i++)
		rel.add(side, p);
}

/* adding relations on the fly in Galois cases */
void add_relations_with_galois(const char *galois, std::ostream& os,
				      const char *comment, unsigned long *cpt,
				      relation &rel)
{
    int64_t a0, b0, a1, b1, a2, b2, a3, b3, a5, b5, aa, bb, a;
    uint64_t b;
    int d;

    a = rel.a; b = rel.b; // should be obsolete one day
    // (a0, b0) = sigma^0((a, b)) = (a, b)
    a0 = rel.a; b0 = runtime_numeric_cast<int64_t>(rel.b);

    if (strcmp(galois, "autom2.1") == 0) {
	// remember, 1/x is for plain autom
	// 1/y is for Galois filtering: x^4+1 -> DO NOT DUPLICATE RELATIONS!
	// (a-b/x) = 1/x*(-b+a*x)
	adwg(os, comment, cpt, rel, -b0, -a0);
    } else if(strcmp(galois, "autom2.2") == 0) {
	// remember, -x is for plain autom
	// -y is for Galois filtering: x^4+1 -> DO NOT DUPLICATE RELATIONS!
	// (a-(-b)*x) ~ (-a-b*x)
	adwg(os, comment, cpt, rel, -a0, b0);
    } else if(strcmp(galois, "autom3.1") == 0) {
	// x -> 1-1/x; hence 1/x*(b-(b-a)*x)
	a1 = a; b1 = (int64_t)b;
	a2 = b1; b2 = b1-a1;
	adwg(os, comment, cpt, rel, a2, b2);
	a3 = b2; b3 = b2-a2;
	adwg(os, comment, cpt, rel, a3, b3);
    } else if(strcmp(galois, "autom3.2") == 0) {
	// x -> -1-1/x; hence 1/x*(b-(-a-b)*x)
	a1 = a; b1 = (int64_t)b;
	a2 = b1; b2 = -a1-b1;
	adwg(os, comment, cpt, rel, a2, b2);
	a3 = b2; b3 = -a2-b2;
	adwg(os, comment, cpt, rel, a3, b3);
    } else if(strcmp(galois, "autom4.1") == 0) {
	// FIXME: rewrite and check
	a1 = a; b1 = (int64_t)b;
	// tricky: sig^2((a, b)) = (2b, -2a) ~ (b, -a)
	aa = b1; bb = -a1;
	if(bb < 0) { aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
	// same factorization as for (a, b)
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
	// sig((a, b)) = (-(a+b), a-b)
	aa = -(a1+b1);
	bb = a1-b1;
	int am2 = a1 & 1, bm2 = b1 & 1;
	if (am2+bm2 == 1) {
	    // (a, b) = (1, 0) or (0, 1) mod 2
	    // aa and bb are odd, aa/bb = 1 mod 2
	    // we must add "2,2" in front of f and g
	    add_galois_factors(rel, 2, 2);
	} else {
	    // (a, b) = (1, 1), aa and bb are even
	    // we must remove "2,2" in front of f and g
	    // taken from relation.cpp
	    remove_galois_factors(rel, 2, 2);
	    // remove common powers of 2
	    do {
		aa >>= 1;
		bb >>= 1;
	    } while((aa & 1) == 0 && (bb & 1) == 0);
	}
	if(bb < 0) { aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
	// sig^3((a, b)) = sig((b, -a)) = (a-b, a+b)
	aa = -aa; // FIXME: check!
	if(aa < 0) { aa = -aa; bb = -bb; }
	rel.a = bb; rel.b = (uint64_t)aa;
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
    } else if(strcmp(galois, "autom6.1") == 0) {
	// fact do not change
	adwg(os, comment, cpt, rel, a0 + b0, -a0); // (a2, b2)
	adwg(os, comment, cpt, rel, b0, -(a0+b0)); // (a4, b4)

	// fact do change
        a1 = -(2*a0+b0); b1= a0-b0;
	d = 0;
	while(((a1 % 3) == 0) && ((b1 % 3) == 0)) {
	    a1 /= 3;
	    b1 /= 3;
	    d++;
	}
	os << "# d1=" << d << "\n";
	a3 =-(2*b0+a0); b3 = 2*a0+b0;
	a5 = a0-b0;     b5 = 2*b0+a0;
	if(d == 0)
	    // we need to add 3^3
	    add_galois_factors(rel, 3, 3);
	else
	    // we need to remove 3^3
	    remove_galois_factors(rel, 3, 3);
	adwg(os, comment, cpt, rel, a1, b1); // (a1/3^d, b1/3^d)
	for(int i = 0; i < d; i++) {
	    a3 /= 3;
	    b3 /= 3;
	    a5 /= 3;
	    b5 /= 3;
	}
	adwg(os, comment, cpt, rel, a3, b3); // (a3/3^d, b3/3^d)
	adwg(os, comment, cpt, rel, a5, b5); // (a5/3^d, b5/3^d)
    }
}

void
skip_galois_roots(cxx_mpz const & q, std::vector<cxx_mpz> &roots,
                  galois_action const &gal_action)
{
    unsigned int const ord = gal_action.get_order();
    unsigned long const qq = mpz_get_ui(q);
    size_t const nroots = roots.size();

    verbose_fmt_print(0, 2, "# galois: got {} root(s) modulo q={}\n", nroots,
                                                                      q);

    /* Keep only one root among sigma-orbits. */
    unsigned long conj[ord]; // where to put conjugates

    char used[nroots]; // used roots: non-principal conjugates
    memset(used, 0, nroots);

    unsigned int n_solo_orbit = 0;
    for (size_t k = 0; k < nroots; k++) {
        if (used[k]) /* already visited during an orbit computation, skip it */
            continue;
        unsigned long rk = mpz_get_ui(roots[k]);
        unsigned long r = rk;

        /* Compute the orbit of roots[k] */
        size_t orbit_length = 0;
        do {
            conj[orbit_length] = r;
            verbose_fmt_print(0, 3, "# galois: orbit of root {} modulo q={} "
                                   "contains {}\n", roots[k], q, r);
            r = gal_action.apply(r, qq);
            orbit_length++;
        } while (r != rk && orbit_length < ord);

        verbose_fmt_print(0, 2, "# galois: orbit of root {} modulo q={} has "
                               "length {}\n", roots[k], q, orbit_length);

        /* Checks:
         *  - sigma^ord(rk) == rk
         *  - length of orbit is 1 or ord
         */
        ASSERT_ALWAYS(r == rk);
        ASSERT_ALWAYS(orbit_length == 1 || orbit_length == ord);

        n_solo_orbit += orbit_length == 1;

        /* Mark all roots of the orbit as visited */
        for (size_t l = k+1; l < nroots; l++) {
            unsigned long const rl = mpz_get_ui(roots[l]);
            for (unsigned int i = 1; i < orbit_length; i++) {
                if (rl == conj[i]) {
                    /* l is the index of an element of the orbit of rk */
                    ASSERT_ALWAYS(used[l] == 0);
                    used[l] = (char)1;
                    break;
                }
            }
        }
    }

    ASSERT_ALWAYS((nroots-n_solo_orbit) % ord == 0);

    /* All orbits are computed, we can compact roots now */
    size_t new_nroots = 0;
    for (size_t k = 0; k < nroots; k++) {
        if (used[k] == 0) { /* roots[k] is the representative of an orbit */
            if(k > new_nroots)
                mpz_swap(roots[new_nroots], roots[k]);
            new_nroots++;
        }
    }
    ASSERT_ALWAYS(new_nroots == n_solo_orbit + (nroots-n_solo_orbit)/ord);
    roots.erase(roots.begin() + new_nroots, roots.end());

    verbose_fmt_print(0, 2, "# galois: computed {} orbits for roots modulo "
                           "q={} (of which {} have length 1).\n", roots.size(),
                           q, n_solo_orbit);
}
