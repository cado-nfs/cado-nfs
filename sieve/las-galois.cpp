#include "cado.h" // IWYU pragma: keep
#include <cstdint>       // for int64_t, uint64_t
#include <cstdio>        // for fprintf, stderr
#include <cstring>       // for strcmp, memset
#include <memory>        // for allocator_traits<>::value_type
#include <ostream>
#include <vector>        // for vector
#include "cxx_mpz.hpp"
#include "galois_utils.h" // automorphism_init
#include "las-galois.hpp"
#include "macros.h"      // for ASSERT_ALWAYS
#include "mod_ul.h"     // modulusul_t
#include "relation.hpp"

static void adwg(std::ostream& os, const char *comment, unsigned long *cpt,
		 relation &rel, int64_t a, int64_t b){
    if(b < 0) { a = -a; b = -b; }
    rel.a = a; rel.b = (uint64_t)b;
    if (comment) os << comment;
    os << rel << '\n';
    *cpt += (*comment != '\0');
}

/* removing p^vp from the list of factors in rel. */
static void remove_galois_factors(relation &rel, int p, int vp){
    int ok = 0;

    for(int side = 0 ; side < 2 ; side++){
	for(unsigned int i = 0 ; i < rel.sides[side].size(); i++)
	    if(mpz_cmp_ui(rel.sides[side][i].p, p) == 0){
		ok = 1;
		ASSERT_ALWAYS(rel.sides[side][i].e >= vp);
		rel.sides[side][i].e -= vp;
	    }
    }
    /* indeed, p was present */
    ASSERT_ALWAYS(ok == 1);
}

/* adding p^vp to the list of factors in rel. */
static void add_galois_factors(relation &rel, int p, int vp){
    int ok[2] = {0, 0};

    for(int side = 0 ; side < 2 ; side++){
	for(unsigned int i = 0 ; i < rel.sides[side].size(); i++)
	    if(mpz_cmp_ui(rel.sides[side][i].p, p) == 0){
		ok[side] = 1;
		rel.sides[side][i].e += vp;
	    }
    }
    // FIXME: are we sure this is safe?
    for(int side = 0 ; side < 2 ; side++)
	if(ok[side] == 0)
	    /* we must add p^vp */
	    for(int i = 0; i < vp; i++)
		rel.add(side, p);
}

/* adding relations on the fly in Galois cases */
void add_relations_with_galois(const char *galois, std::ostream& os,
				      const char *comment, unsigned long *cpt,
				      relation &rel){
    int64_t a0, b0, a1, b1, a2, b2, a3, b3, a5, b5, aa, bb, a;
    uint64_t b;
    int d;

    a = rel.a; b = rel.b; // should be obsolete one day
    // (a0, b0) = sigma^0((a, b)) = (a, b)
    a0 = rel.a; b0 = (int64_t)rel.b;
    if(strcmp(galois, "autom2.1") == 0)
	// remember, 1/x is for plain autom
	// 1/y is for Galois filtering: x^4+1 -> DO NOT DUPLICATE RELATIONS!
	// (a-b/x) = 1/x*(-b+a*x)
	adwg(os, comment, cpt, rel, -b0, -a0);
    else if(strcmp(galois, "autom2.2") == 0)
	// remember, -x is for plain autom
	// -y is for Galois filtering: x^4+1 -> DO NOT DUPLICATE RELATIONS!
	// (a-(-b)*x) ~ (-a-b*x)
	adwg(os, comment, cpt, rel, -a0, b0);
    else if(strcmp(galois, "autom3.1") == 0){
	// x -> 1-1/x; hence 1/x*(b-(b-a)*x)
	a1 = a; b1 = (int64_t)b;
	a2 = b1; b2 = b1-a1;
	adwg(os, comment, cpt, rel, a2, b2);
	a3 = b2; b3 = b2-a2;
	adwg(os, comment, cpt, rel, a3, b3);
    }
    else if(strcmp(galois, "autom3.2") == 0){
	// x -> -1-1/x; hence 1/x*(b-(-a-b)*x)
	a1 = a; b1 = (int64_t)b;
	a2 = b1; b2 = -a1-b1;
	adwg(os, comment, cpt, rel, a2, b2);
	a3 = b2; b3 = -a2-b2;
	adwg(os, comment, cpt, rel, a3, b3);
    }
    else if(strcmp(galois, "autom4.1") == 0){
	// FIXME: rewrite and check
	a1 = a; b1 = (int64_t)b;
	// tricky: sig^2((a, b)) = (2b, -2a) ~ (b, -a)
	aa = b1; bb = -a1;
	if(bb < 0){ aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
	// same factorization as for (a, b)
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
	// sig((a, b)) = (-(a+b), a-b)
	aa = -(a1+b1);
	bb = a1-b1;
	int am2 = a1 & 1, bm2 = b1 & 1;
	if(am2+bm2 == 1){
	    // (a, b) = (1, 0) or (0, 1) mod 2
	    // aa and bb are odd, aa/bb = 1 mod 2
	    // we must add "2,2" in front of f and g
	    add_galois_factors(rel, 2, 2);
	}
	else{
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
	if(bb < 0){ aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
	// sig^3((a, b)) = sig((b, -a)) = (a-b, a+b)
	aa = -aa; // FIXME: check!
	if(aa < 0){ aa = -aa; bb = -bb; }
	rel.a = bb; rel.b = (uint64_t)aa;
        if (comment) os << comment;
        os << rel << '\n';
        *cpt += (*comment != '\0');
    }
    else if(strcmp(galois, "autom6.1") == 0){
	// fact do not change
	adwg(os, comment, cpt, rel, a0 + b0, -a0); // (a2, b2)
	adwg(os, comment, cpt, rel, b0, -(a0+b0)); // (a4, b4)

	// fact do change
        a1 = -(2*a0+b0); b1= a0-b0;
	d = 0;
	while(((a1 % 3) == 0) && ((b1 % 3) == 0)){
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
	for(int i = 0; i < d; i++){
	    a3 /= 3;
	    b3 /= 3;
	    a5 /= 3;
	    b5 /= 3;
	}
	adwg(os, comment, cpt, rel, a3, b3); // (a3/3^d, b3/3^d)
	adwg(os, comment, cpt, rel, a5, b5); // (a5/3^d, b5/3^d)
    }
}

int
skip_galois_roots(const int orig_nroots, const mpz_t q, mpz_t *roots,
		  const char *galois_autom)
{
    int imat[4];
    residueul_t mat[4];
    int nroots = orig_nroots, ord;

    if(nroots == 0)
	return 0;
    automorphism_init(&ord, imat, galois_autom);
    modulusul_t mm;
    unsigned long qq = mpz_get_ui(q);
    modul_initmod_ul(mm, qq);
    for(int i = 0; i < 4; i++){
	modul_init(mat[i], mm);
	modul_set_int64(mat[i], imat[i], mm);
    }
    if (nroots % ord) {
        fprintf(stderr, "Number of roots modulo q is not divisible by %d. Don't know how to interpret -galois.\n", ord);
        ASSERT_ALWAYS(0);
    }
    // Keep only one root among sigma-orbits.
    residueul_t r2, r3;
    modul_init(r2, mm);
    modul_init(r3, mm);
    residueul_t conj[ord]; // where to put conjugates
    for(int k = 0; k < ord; k++)
	modul_init(conj[k], mm);
    char used[nroots];     // used roots: non-principal conjugates
    memset(used, 0, nroots);
    for(int k = 0; k < nroots; k++){
	if(used[k]) continue;
	unsigned long rr0 = mpz_get_ui(roots[k]), rr;
	rr = rr0;
	// build ord-1 conjugates for roots[k]
	for(int l = 0; l < ord; l++){
	    rr = automorphism_apply(mat, rr, mm, qq);
	    modul_set_ul(conj[l], rr, mm);
	}
#if 0 // debug. 
	printf("orbit for %lu: %lu", qq, rr);
	for(int l = 0; l < ord-1; l++)
	    printf(" -> %lu", conj[l][0]);
	printf("\n");
#endif
	// check: sigma^ord(rr0) should be rr0
	ASSERT_ALWAYS(rr == rr0);
	// look at roots
	for(int l = k+1; l < nroots; l++){
	    unsigned long ss = mpz_get_ui(roots[l]);
	    modul_set_ul(r2, ss, mm);
	    for(int i = 0; i < ord-1; i++)
		if(modul_equal(r2, conj[i], mm)){
		    ASSERT_ALWAYS(used[l] == 0);
		    // l is some conjugate, we erase it
		    used[l] = (char)1;
		    break;
		}
	}
    }
    // now, compact roots
    int kk = 0;
    for(int k = 0; k < nroots; k++)
	if(used[k] == 0){
	    if(k > kk)
		mpz_set(roots[kk], roots[k]);
	    kk++;
	}
    ASSERT_ALWAYS(kk == (nroots/ord));
    nroots = kk;
    for(int k = 0; k < ord; k++)
	modul_clear(conj[k], mm);
    for(int i = 0; i < 4; i++)
	modul_clear(mat[i], mm);
    modul_clear(r2, mm);
    modul_clear(r3, mm);
    modul_clearmod(mm);
    return nroots;
}

