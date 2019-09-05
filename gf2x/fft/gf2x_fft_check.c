/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2014, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#define _XOPEN_SOURCE   600
#define _POSIX_C_SOURCE   600
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#ifdef  _POSIX_C_SOURCE
#include <sys/types.h>
#include <unistd.h>
#endif

#include "gf2x.h"

#ifdef GF2X_SOURCE_TREE
#include "fft/gf2x-fft.h"
#else
#include "gf2x/gf2x-fft.h"
#endif

#ifndef MAX
#define MAX(a,b)        ((a)<(b) ? (b) : (a))
#endif

#ifndef MIN
#define MIN(a,b)        ((a)<(b) ? (a) : (b))
#endif

/* Handy, and does not require libm */
#ifndef iceildiv
/* C99 defines division as truncating, i.e. rounding to 0. To get ceil in all
   cases, we have to distinguish by sign of arguments. Note that (afaik) the
   sign and truncation direction of the quotient with any negative operads was
   undefined in C89. */
/* iceildiv requires *unsigned* operands in spite of the name suggesting
   (signed) integer type. For negative operands, the result is wrong. */
#define iceildiv(x,y) ((x) == 0 ? 0 : ((x)-1)/(y)+1)
#endif

#ifndef ULONG_BITS
#define ULONG_BITS      ((int) (sizeof(unsigned long) * CHAR_BIT))
#endif

/* Checks the different FFT codes, and see whether they give any
 * consistent results.
 *
 * We compute f1g1+f2g2 for several sizes of f1, f2, using all the
 * available algorithms.
 *
 *
 * This tests _only_ the functionality in the libgf2x-fft that is
 * eventually shipped.
 *
 * Note that this part of the library is still clearly development work.
 *
 * When this fail, we print the polynomials computed as lists of 64-bit
 * (or 32-bit) integers, to be parsed with magma by the companion code
 * that gets printed alongside the computation if the -m arg is passed to
 * the program.
 *
 */

unsigned long next_seed;

/* When on, all computed data is printed in magma format with asserts so
 * that it is possible to trace exactly what happens throughout.
 *
 * A priori, one does this only to inspect bugs, and only for tests that
 * have been recognized to fail once.
 *
 * It is also possible to do some wide-range testing with:
 *
 * ./fft/gf2x_fft_check -m | magma -b
 *
 */
int magma = 0;

int doing_mp = 0;

void usage()
{
    fprintf(stderr, "Usage: ./gf2x_fft_check [-a <max bits of a>] [-b <max bits of b>] [-nr <nreps per size>] [-seed <random seed>] [-m] [-mul] [-mp]\n");
}

void display(const char * text, unsigned long * f, size_t nf)
{
    printf("%s:=[", text);
    for(size_t i = 0 ; i < nf ; i++) {
        if (i) printf(",");
        printf("%lu",f[i]);
    }
    printf("];\n");
}

/* returns an unsigned long with long strings of zeroes and ones */
unsigned long rand2_ulong()
{
    /* for 64-bit ulongs, we need 1 + nrounds * 6 bits as given by
     * random(), which is ok until nchunks == 5.
     *
     * Note that we don't care much about randomness in the strict sense
     * (in case anybody still wonders).
     */
    unsigned int r = random();
    unsigned long bit = r & 1UL;
    unsigned long res = -r;
    unsigned int k = 0;
    for( ; (1 << k) < ULONG_BITS ; k++);
    r >>= 1;
    for(unsigned int i = 0 ; i < 4 ; i++) {
        unsigned int chunksize = 1 + (r & (ULONG_BITS - 1));
        bit = !bit;
        res += bit;
        r >>= k;
        res <<= chunksize;
        res -= bit;
    }
    return res;
}

/* Most of this code is made of macros. It should really be done with C++
 * templates, but gf2x has no dependency on c++ for the moment, and I
 * don't feel like adding one. Anyway, bugs reside more often in the code
 * being tested (the library) than in the test code itself.
 */
#define CHOP_HEAD(P) do {                                               \
        for(size_t i = 0 ; i < nw ## P ; i++)                           \
            P[i] = rand2_ulong();                                       \
        if (n ## P % ULONG_BITS)                                        \
            P[nw ## P-1] &= (1UL << (n ## P % ULONG_BITS)) - 1; \
} while (0)

#define ALLOC1(E,x) E ## _ptr E ## _ ## x = E ## _alloc(E, 1)

#define SETUP(E)                                         \
    unsigned long * E ## _h = malloc(nwh * sizeof(unsigned long));      \
    unsigned long * E ## _hx = malloc(nwhx * sizeof(unsigned long));    \
    E ## _info_copy(E ## _copy, E);                             \
    size_t E ## _fft_sizes[3];                                  \
    E ## _info_get_alloc_sizes(E, E ## _fft_sizes);             \
    E ## _ptr E ## _temp1 = malloc(E ## _fft_sizes[1]);         \
    E ## _ptr E ## _temp2 = malloc(E ## _fft_sizes[2]);         \
    ALLOC1(E, tf1); ALLOC1(E, tg1);                             \
    ALLOC1(E, tf2); ALLOC1(E, tg2);                             \
    E ## _ptr E ## _tfs[2] = { E ## _tf1, E ## _tf2, };         \
    E ## _ptr E ## _tgs[2] = { E ## _tg1, E ## _tg2, };         \
    ALLOC1(E, th);                                              \
    ALLOC1(E, thx);                                             \
    do { } while (0)

#define ENTER(E, OP, nf, ng)                                    \
    E ## _info_t E, E ## _copy;                                 \
    E ## _ ## OP (E, nf, ng);                                   \
    do { } while (0)

#define DO_zero_th(E) E ## _zero(E, E ## _th, 1)
#define DO_zero_thx(E) E ## _zero(E, E ## _thx, 1)

#define display_full_poly_and_lift(E, P) do {                           \
        display(#P,P,nw ## P);                                          \
        printf("K" #P ":=Zseq_to_KP(" #P ");\n");                       \
        printf("L" #P ":=KP_to_LP(K" #P ");\n");                        \
} while (0)

#define display_projected_result(E, P) do {                             \
        display(#E "_" #P,E ## _ ## P,nw ## P);                         \
        printf("K" #P ":=Zseq_to_KP(" #E "_" #P ");\n");                \
} while (0)

#define display_full_transform(E, P) do {                               \
        display(#E "_t" #P,                                             \
            (unsigned long *) E ## _t ## P,                             \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long));   \
        printf("T" #P ":=Zseq_to_Lseq(" #E "_t" #P ");\n");             \
} while (0)

#define DO_dft(E, P) do {                                               \
    if (magma) {                                                        \
        display_full_poly_and_lift(E, P);                               \
    }                                                                   \
    E ## _dft(E, E ## _t ## P, P, n ## P, E ## _temp1);                 \
    if (magma) {                                                        \
        display_full_transform(E, P);                                   \
        printf("assert T" #P " eq L_evaluate(L" #P ", transform_length(T" #P "));\n");    \
    }                                                                   \
} while (0)

// &and [T" #P "[i+1] eq Evaluate(L" #P ", evalpoint(i)):i in [0..#T" #P "-1]];\n");                                             

#define DO_dft_f1(E) DO_dft(E, f1)
#define DO_dft_f2(E) DO_dft(E, f2)
#define DO_dft_g1(E) DO_dft(E, g1)
#define DO_dft_g2(E) DO_dft(E, g2)

#define DO_comp(E, k) do {                                              \
    if (magma) {                                                        \
        display_full_transform(E, h);                                   \
        printf("previous_Th:=Th;\n");                                   \
    }                                                                   \
    E ## _addcompose(E, E ## _th,                                       \
            (E ## _srcptr) E ## _tf ## k,                               \
            (E ## _srcptr) E ## _tg ## k,                               \
            E ## _temp2, E ## _temp1);                                  \
    if (magma) {                                                        \
        display_full_transform(E, h);                                   \
        printf("assert Th eq transform_add(previous_Th,"                \
                    " transform_pointwise(Tf"#k",Tg"#k"));\n");         \
    }                                                                   \
} while (0)

#define DO_comp1(E)  DO_comp(E, 1)
#define DO_comp2(E)  DO_comp(E, 2)
#define DO_comp_n(E)  do {                                              \
    E ## _addcompose_n(E, E ## _thx,                                    \
            (E ## _srcptr *) E ## _tfs, (E ## _srcptr *) E ## _tgs, 2,  \
            E ## _temp2, E ## _temp1);                                  \
    if (magma) {                                                        \
        display_full_transform(E, hx);                                  \
        printf("assert Thx eq transform_add("                           \
                "transform_pointwise(Tf1,Tg1),"                         \
                "transform_pointwise(Tf2,Tg2));\n");                    \
    }                                                                   \
} while (0)

/* Interpolation is known to be tricky because we blend it with recovery
 * of the truncation operation. In the mul case, interpolation uses the
 * fact that h is known to have at most nh bits. In the mp case, that is
 * largely unsatisfactory at the moment, since we barely truncate in the
 * Cantor case.
 */
#define DO_ift(E, P)  do {                                              \
    if (magma) {                                                        \
        printf("saved_T" #P ":=T" #P ";\n");                            \
    }                                                                   \
    E ## _ift(E, E ## _ ## P, n ## P, E ## _t ## P, E ## _temp1);       \
    if (magma) {                                                        \
        printf("/* " #E "_th has has undergone ift now */\n");          \
        display_full_transform(E, P);                                   \
        display_projected_result(E, P);                                 \
        printf("interpolated_" #P ":=Lseq_as_LP(T" #P ");\n");          \
        printf("assert L_evaluate(interpolated_"#P","                   \
                "transform_length(T"#P")) eq "                          \
                "saved_T"#P";\n");                                      \
        printf("assert interpolated_"#P " eq "                          \
                "L_interpolate(saved_T"#P ", "                          \
                "transform_length(T"#P"));\n");                         \
        printf("assert interpolated_" #P " eq "                         \
                "LP_canonical(LP_add(LP_mul(Lf1, Lg1),"                 \
                " LP_mul(Lf2, Lg2)));\n");                              \
        /* When doing MP, we tolerate modular equality here */          \
        if (doing_mp) {                                                 \
            printf("if assigned Kwrap then\n"   \
                        "\tassert (LP_to_KP(interpolated_" #P ") - "    \
                        "(Kf1*Kg1+Kf2*Kg2)) mod Kwrap eq 0;\n"          \
                    "else\n"                                            \
                        "\tassert LP_to_KP(interpolated_" #P ") eq "    \
                        "Kf1*Kg1+Kf2*Kg2;\n"                            \
                    "end if;\n");                                       \
            printf("assert LP_to_KP(interpolated_" #P ") "              \
                    "mod x^%zu div x^%zu "                              \
                    "eq Zseq_to_KP(" #E "_h);\n",                       \
                    MAX(nf,ng), MIN(nf,ng)-1);                          \
        } else {                                                        \
            printf("assert LP_to_KP(interpolated_" #P ") eq "           \
                    "Kf1*Kg1+Kf2*Kg2;\n");                              \
        }                                                               \
    }                                                                   \
} while (0)

#define DO_ift_h(E)  DO_ift(E, h)
#define DO_ift_hx(E) DO_ift(E, hx)

#define EXTRA_DISPLAY(E)  do {                                          \
    display(#E "_tf1", (unsigned long *) E ## _tf1,                     \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long));   \
    display(#E "_tg1", (unsigned long *) E ## _tg1,                       \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long)); \
    display(#E "_tf2", (unsigned long *) E ## _tf2,                       \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long)); \
    display(#E "_tg2", (unsigned long *) E ## _tg2,                       \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long)); \
} while (0)
// Either th and thx are different, meaning that addcompose and
// addcompose_n disagree -- in which case we want to print them both.
// Or they're equal, and we've proceeded to computing the ift, which
// invalidates th and thx in most cases. In the latter case, we're of
// course better off not printing any data.
#define EXTRA_DISPLAY_TH_THX(E)  do {                                           \
    display(#E "_th",  (unsigned long *) E ## _th,                      \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long)); \
    display(#E "_thx",  (unsigned long *) E ## _thx,                    \
            E ## _transform_size(E) * sizeof(E ## _elt) / sizeof(unsigned long));   \
} while (0)

#define FREE1(E, x) do {                                                \
        E ## _free(E, E ## _ ## x, 1);                                  \
    } while (0)

#define LEAVE(E) do {                                                   \
    FREE1(E, tf1); FREE1(E, tg1);                                       \
    FREE1(E, tf2); FREE1(E, tg2);                                       \
    free(E ## _temp1);                                                  \
    free(E ## _temp2);                                                  \
    FREE1(E, th);                                                       \
    FREE1(E, thx);                                                      \
    E ## _info_clear(E);                                                \
    E ## _info_clear(E ## _copy);                                       \
    free(E ## _h);                                                      \
    free(E ## _hx);                                                     \
} while (0)

#define CHECK_SELF_ADDCOMPOSE_N_CONSISTENCY(E) do {                     \
    if (memcmp(E ## _h, E ## _hx, nwh * sizeof(unsigned long)) != 0) {  \
        fprintf(stderr, #E " != " #E "x "                               \
                "(addcompose_n consistency) for "                       \
                "-%s -a %zu -b %zu -seed %u\n",                         \
                __func__ + 8,                                           \
                nf,ng,seed);                                            \
        printf("w:=%d;\n", (int) ULONG_BITS);                           \
        display("f1",f1,nwf1);                                          \
        display("g1",g1,nwg1);                                          \
        display("f2",f2,nwf2);                                          \
        display("g2",g2,nwg2);                                          \
        display(#E "_h",E ## _h,nwh);                                   \
        display(#E "_hx",E ## _hx,nwh);                         \
        EXTRA_DISPLAY(E);                                               \
        EXTRA_DISPLAY_TH_THX(E);                                        \
        fflush(stdout);                                                 \
        fflush(stderr);                                                 \
        abort();                                                        \
    }                                                                   \
} while (0)

#define CHECK_CROSS_CONSISTENCY(E, F) do {                              \
    if (memcmp(F ## _h, E ## _h, nwh * sizeof(unsigned long)) != 0) {   \
        fprintf(stderr, #F " != " #E " "                                \
                "(cross check) for "                                    \
                "-%s -a %zu -b %zu -seed %u\n",                         \
                __func__ + 8,                                           \
                nf,ng,seed);                                            \
        printf("w:=%d;\n", (int) ULONG_BITS);                           \
        display("f1",f1,nwf1);                                          \
        display("g1",g1,nwg1);                                          \
        display("f2",f2,nwf2);                                          \
        display("g2",g2,nwg2);                                          \
        display(#F "_h",F ## _h,nwh);                                   \
        display(#E "_h",E ## _h,nwh);                                   \
        EXTRA_DISPLAY(F);                                               \
        EXTRA_DISPLAY(E);                                               \
        fflush(stdout);                                                 \
        fflush(stderr);                                                 \
        abort();                                                        \
    }                                                                   \
} while (0)

/* There are four different data types at stake.
 *
 *  - list of integers (Zseq). This is what the code prints. These objects may
 *  be interpreted as one of the next three types.
 *
 *  - polynomials (KP). The function that converts from lists of integers to
 *  polynomials is Zseq_to_KP. The reverse function is KP_to_Zseq
 *
 *  - lifted polynomials (LP), to some other structure. The conversion from
 *  polynomials to lifted polynomials is done by "KP_to_LP", and the reverse
 *  operation is done by "LP_to_KP"
 *  (however pay attention that the mapping is not bijective -- several
 *  LP polynomials may correspond to the same KP polynomial).
 *
 *  - list of evaluations in the other structure (Lseq). The conversion from
 *  list of integers to list of evaluations is done by
 *  Zseq_to_Lseq.
 *
 */
/* unKS : useful ?

KS-based polynomials may be interpreted by:

unKS:=func<X|[Zseq_to_KP(X[i..i+1]):i in [1..#X by 2]]>;

I must pay attention to the fact that data undergoes FFT early on, so
that an unKS operation may be meaningful only at precise points in the
process.

*/

/* A generic recursive polyfromeval in magma goes as follows. The
 * internal functions works on scaled evaluations, and returns the
 * interpolated polynomial as well as the polyfromroots. */

/* Note that this does Lagrange interpolation, which is a bit more
 * than what the backend code might be doing. E.g. in Cantor, the
 * raw output of the non-truncated transform is bigger, and it takes some
 * work to reach the low degree polynomial that we get here.
 */
const char * interpolate_code =
"function interpolate_generic_scaled(ev, roots)\n"
"    assert #ev eq #roots;\n"
"    n:=#ev;\n"
"    if n eq 1 then\n"
"        interpolated:=Polynomial([ev[1]]);\n"
"        return interpolated, Parent(interpolated).1-roots[1];\n"
"    else\n"
"        nl := n div 2;\n"
"        il,pl := interpolate_generic_scaled(ev[1..nl], roots[1..nl]);\n"
"        ih,ph := interpolate_generic_scaled(ev[nl+1..n], roots[nl+1..n]);\n"
"        return il * ph + ih * pl, pl * ph;\n"
"    end if;\n"
"end function;\n"
"function polyfromroots_generic(roots)\n"
"    n:=#roots;\n"
"    if n eq 1 then\n"
"        return Polynomial([-roots[1], 1]);\n"
"    else\n"
"        nl := n div 2;\n"
"        pl := polyfromroots_generic(roots[1..nl]);\n"
"        ph := polyfromroots_generic(roots[nl+1..n]);\n"
"        return pl * ph;\n"
"    end if;\n"
"end function;\n"
"function interpolate_generic(ev, roots)\n"
"    assert #ev eq #roots;\n"
"    n:=#ev;\n"
"    p := polyfromroots_generic(roots);\n"
"    dp := Derivative(p);\n"
"    ev_scaled := [ ev[i] / Evaluate(dp, roots[i]) : i in [1..n]];\n"
"    interpolated, pfr := interpolate_generic_scaled(ev_scaled, roots);\n"
"    assert p eq pfr;\n"
"    return PolynomialRing(Universe(ev))!interpolated;\n"
"end function;\n"
;

void print_context_gf2x_fake_fft(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED)
{
    if (!magma) return;
    fputs(
            "\n\n/* fake fft (transform == identity) */\n"
            "clear;\n"
            , stdout);
    printf("w:=%d;\n", (int) ULONG_BITS);
    fputs(
            "KP<x>:=PolynomialRing(GF(2));\n"
            "Zseq_to_KP:=func<f|Polynomial(GF(2),&cat[Intseq(x,2,w):x in f])>;\n"
            "KP_to_Zseq:=func<P,n|Intseq(Seqint(ChangeUniverse(Eltseq(P),Integers()),2),2^w,Ceiling(n/w))>;\n"
            "L:=KP;\n"
            "ks:=-1;\n"
            "KP_to_LP:=func<x|x>;\n"
            "LP_to_KP:=func<x|x>;\n"
            "Zseq_to_Lseq:=func<x|[Zseq_to_KP(x)]>;\n"
            "Lseq_to_Zseq:=func<x|KP_to_Zseq(x[1],1+Degree(x[1]))>;\n"
            "evalpoint:=func<i|x>;\n"
            "evalpoints:=[x];\n"
            "L_evaluate:=func<P,n|[P]>;\n"
            "L_interpolate:=func<x,n|x[1]>;\n"
            "transform_length:=func<x|1>;\n"
            "transform_add:=func<a,b|[a[1]+b[1]]>;\n"
            "transform_pointwise:=func<a,b|[a[1]*b[1]]>;\n"
            "Lseq_as_LP:=func<x|Polynomial(x)>;\n"
            "LP_canonical:=func<x|x>;\n"
            "LP_add:=func<a,b|a+b>;LP_mul:=func<a,b|a*b>;\n"
            , stdout);
}

void print_context_gf2x_cantor_fft(gf2x_cantor_fft_info_srcptr p GF2X_MAYBE_UNUSED)
{
    if (!magma) return;
    fputs(
        "\n\n/* Cantor additive fft */\n"
        "clear;\n"
        "KP<x>:=PolynomialRing(GF(2));\n"
        , stdout);
#if (CANTOR_BASE_FIELD_SIZE == 128)
        fputs("L<z> := ext<GF(2) | x^128 + x^7 + x^2 + x + 1>;\n", stdout);
#elif (CANTOR_BASE_FIELD_SIZE == 64)
        fputs("L<z> := ext<GF(2) | x^64 + x^4 + x^3 + x + 1>;\n", stdout);
#else
#error "CANTOR_BASE_FIELD_SIZE must be either 0 or 64"
#endif
    printf("w:=%d;\n", (int) ULONG_BITS);
    printf("fft_k:=%d;\n",p->k);
    fputs(
        "ks:=Degree(L) div 2;\n"
        "Zseq_to_KP:=func<f|Polynomial(GF(2),&cat[Intseq(x,2,w):x in f])>;\n"
        "KP_to_Zseq:=func<P,n|Intseq(Seqint(ChangeUniverse(Eltseq(P),Integers()),2),2^w,Ceiling(n/w))>;\n"
        "LP<X> := PolynomialRing(L);\n"
        "KS:=func<X|Polynomial([KP|Polynomial(L[i..Minimum(#L,i+ks-1)]) : i in [1..#L by ks]]) where L is Eltseq(X)>;\n"
        "KP_to_LP:=func<f|Polynomial([L|Evaluate(c,L.1):c in Coefficients(KS(f))])>;\n"
        "LP_to_KP:=func<P|Evaluate(Polynomial([KP|Polynomial(L!c):c in Coefficients(P)]), KP.1^ks)>;\n"
        "beta:=L!1;\n"
        "betas:=[beta];\n"
        "for i in [1..31] do\n"
        "    ti:=Polynomial([beta,1,1]);\n"
        "    beta:=[r[1]:r in Roots(ti) | Eltseq(r[1])[1] eq 0][1];\n"
        "    Append(~betas, beta);\n"
        "end for;\n"
        "betas:=Vector(betas);\n"
        "evalpoint:=func<i|(Vector(ChangeUniverse(Intseq(i,2,32),L)),betas)>;\n"
        "evalpoints:=[evalpoint(i-1):i in [1..2^fft_k]];\n"
        "transform_length:=func<x|#x>;\n"
        "transform_add:=func<a,b|[a[i]+b[i]:i in [1..#a]]>;\n"
        "transform_pointwise:=func<a,b|[a[i]*b[i]:i in [1..#a]]>;\n"
        "Lseq_as_LP:=func<x|Polynomial(x)>;\n"
        "LP_canonical:=func<x|x>;\n"
        "LP_add:=func<a,b|a+b>;LP_mul:=func<a,b|a*b>;\n"
        "L_evaluate:=func<P,n|[Evaluate(P,x):x in evalpoints[1..n]]>;\n"
        "Zseq_to_Lseq:=func<X|[Evaluate(Zseq_to_KP(X[i..i+(Degree(L) div w)-1]),L.1):i in [1..#X by (Degree(L) div w)]]>;\n"
        "Lseq_to_Zseq:=func<X|&cat [KP_to_Zseq(Polynomial(c),Degree(L)):c in X]>;\n"
        , stdout);
    fputs(interpolate_code, stdout);
    fputs("L_interpolate:=func<P,n|interpolate_generic(P, evalpoints[1..n])>;\n", stdout);
    if (doing_mp) {
        /* Note: binomial(k, i) mod 2 is 1 if and only if i&~k==0
         * (Kummer)
         */
        printf("Lwrap:=&+[X^(2^i)*Binomial(fft_k,i):i in [0..fft_k]];\n");
        printf("Kwrap:=&+[x^(ks*2^i)*Binomial(fft_k,i):i in [0..fft_k]];\n");
        printf("LP_canonical:=func<x|x mod Lwrap>;\n");
    }
}

void print_context_gf2x_ternary_fft(gf2x_ternary_fft_info_srcptr p GF2X_MAYBE_UNUSED)
{
    if (!magma) return;
    fputs(
            "\n\n/* Schonhage ternary fft */\n"
            "clear;\n"
            , stdout);
    printf("w:=%d;\n", (int) ULONG_BITS);
    if (p->K == 0) {
        printf("/* fall back to plain */\n");
        print_context_gf2x_fake_fft(NULL);
        return;
    }
    if (p->split == 0) {
        printf("fft_M:=%zu; fft_K:=%zu; fft_split:=%d;\n",
                p->M, p->K, p->split);
        printf("fft_Np:=Ceiling(fft_M / (fft_K div 3)) * (fft_K div 3);\n");
        printf("ks:=fft_M;\n");
        printf("fft_np:=Ceiling(fft_Np/w);fft_2np:=2*fft_np;\n");
        fputs(
                "KP<x>:=PolynomialRing(GF(2));\n"
                "Zseq_to_KP:=func<f|Polynomial(GF(2),&cat[Intseq(x,2,w):x in f])>;\n"
                "KP_to_Zseq:=func<P,n|Intseq(Seqint(ChangeUniverse(Eltseq(P),Integers()),2),2^w,Ceiling(n/w))>;\n"
                "L<z> := quo<KP | x^(2*fft_Np) + x^fft_Np + 1>;\n"
                "LP<X> := PolynomialRing(L);\n"
                "KS:=func<X|Polynomial([KP|Polynomial(L[i..Minimum(#L,i+ks-1)]) : i in [1..#L by ks]]) where L is Eltseq(X)>;\n"
                "KP_to_LP:=func<f|Polynomial([L|Evaluate(c,L.1):c in Coefficients(KS(f))])>;\n"
                /* We really need Eltseq first for LP_to_KP */
                "LP_to_KP:=func<P|Evaluate(Polynomial([KP|Polynomial(Eltseq(L!c)):c in Coefficients(P)]), KP.1^ks)>;\n"
                "tritrev:=func<x|Seqint(Reverse(Intseq(x,3,Ilog(3,fft_K))),3)>;\n"
                "evalpoint:=func<i|zeta^tritrev(i)> where zeta is L.1^(fft_Np div (fft_K div 3));\n"
                "evalpoints:=[evalpoint(i-1):i in [1..fft_K]];\n"
                "transform_length:=func<x|#x>;\n"
                "transform_add:=func<a,b|[a[i]+b[i]:i in [1..#a]]>;\n"
                "transform_pointwise:=func<a,b|[a[i]*b[i]:i in [1..#a]]>;\n"
            "Lseq_as_LP:=func<x|Polynomial(x)>;\n"
            "LP_canonical:=func<x|x>;\n"
            "LP_add:=func<a,b|a+b>;LP_mul:=func<a,b|a*b>;\n"
                "L_evaluate:=func<P,n|[Evaluate(P,x):x in evalpoints[1..n]]>;\n"
                "Zseq_to_Lseq:=func<X|[Evaluate(Zseq_to_KP(X[i..i+fft_2np-1]),L.1):i in [1..#X by fft_2np]]>;\n"
                "Lseq_to_Zseq:=func<X|&cat [KP_to_Zseq(Polynomial(Eltseq(c)),2*fft_Np):c in X]>;\n"
                , stdout);
        fputs(interpolate_code, stdout);
        fputs("L_interpolate:=func<P,n|interpolate_generic(P, evalpoints[1..n])>;\n", stdout);
        if (doing_mp) {
            printf("Lwrap:=X^fft_K-1;\n");
            printf("Kwrap:=x^(ks * fft_K)-1;\n");
            printf("LP_canonical:=func<x|x mod Lwrap>;\n");
        }
    } else {
        printf("fft_M1:=%zu; fft_M2:=%zu; fft_K:=%zu; fft_split:=%d;\n",
                p->M, p->M-1, p->K, p->split);

        printf("fft_Np1:=Ceiling(fft_M1 / (fft_K div 3)) * (fft_K div 3);\n");
        printf("fft_Np2:=Ceiling(fft_M2 / (fft_K div 3)) * (fft_K div 3);\n");
        printf("ks1:=fft_M1;\n");
        printf("ks2:=fft_M2;\n");
        printf("fft_np1:=Ceiling(fft_Np1 / w);fft_2np1:=2*fft_np1;\n");
        printf("fft_np2:=Ceiling(fft_Np2 / w);fft_2np2:=2*fft_np2;\n");
        fputs(
                "KP<x>:=PolynomialRing(GF(2));\n"
                "Zseq_to_KP:=func<f|Polynomial(GF(2),&cat[Intseq(x,2,w):x in f])>;\n"
                "KP_to_Zseq:=func<P,n|Intseq(Seqint(ChangeUniverse(Eltseq(P),Integers()),2),2^w,Ceiling(n/w))>;\n"
                "L1<z> := quo<KP | x^(2*fft_Np1) + x^fft_Np1 + 1>;\n"
                "L2<z> := quo<KP | x^(2*fft_Np2) + x^fft_Np2 + 1>;\n"
                "LP1<X> := PolynomialRing(L1);\n"
                "LP2<X> := PolynomialRing(L2);\n"
                "KS1:=func<X|Polynomial([KP|Polynomial(L[i..Minimum(#L,i+ks1-1)]) : i in [1..#L by ks1]]) where L is Eltseq(X)>;\n"
                "KS2:=func<X|Polynomial([KP|Polynomial(L[i..Minimum(#L,i+ks2-1)]) : i in [1..#L by ks2]]) where L is Eltseq(X)>;\n"
                "KP_to_LP:=func<f|<Polynomial([L1|Evaluate(c,L1.1):c in Coefficients(KS1(f))]),Polynomial([L2|Evaluate(c,L2.1):c in Coefficients(KS2(f))])>>;\n"
                /* We really need Eltseq first for LP_to_KP */
                "function split_reconstruct(N,P1,M1,P2,M2)\n"
                "    assert Parent(P1) eq KP;\n"
                "    assert Parent(P2) eq KP;\n"
                "    assert M1 gt M2;\n"
                "    assert M2 gt N/2;\n"
                "    delta:=M1-M2;\n"
                "    P:=P1;\n"
                "    for i in [N-1-M1..0 by -1] do\n"
                "        c:=Coefficient(P,i+delta) + Coefficient(P2,i+delta);\n"
                "        P +:= c * x^i * (1 + x^M1);\n"
                "    end for;\n"
                "    return P;\n"
                "end function;\n"
                "LP_to_KP:=func<P|\n"
                "    split_reconstruct(\n"
                "        2 * fft_M2 * fft_K - 1,\n"
                "        Evaluate(Polynomial([KP|\n"
                "                Polynomial(Eltseq(L1!c))\n"
                "                :c in Coefficients(P[1])]), KP.1^ks1)"
                "                mod (1+x^(fft_M1*fft_K)),\n"
                "        fft_M1 * fft_K,\n"
                "        Evaluate(Polynomial([KP|\n"
                "                Polynomial(Eltseq(L2!c))\n"
                "                :c in Coefficients(P[2])]), KP.1^ks2)"
                "                mod (1+x^(fft_M2*fft_K)),\n"
                "        fft_M2 * fft_K)"
                "        >;\n"
                "tritrev:=func<x|Seqint(Reverse(Intseq(x,3,Ilog(3,fft_K))),3)>;\n"
                "evalpoint1:=func<i|zeta1^tritrev(i)> where zeta1 is L1.1^(fft_Np1 div (fft_K div 3));\n"
                "evalpoint2:=func<i|zeta2^tritrev(i)> where zeta2 is L2.1^(fft_Np2 div (fft_K div 3));\n"
                "transform_length:=func<x|#x[1]>;\n"
                "transform_add:=func<a,b|<\n"
                    "[a[1][i]+b[1][i]:i in [1..#a[1]]],\n"
                    "[a[2][i]+b[2][i]:i in [1..#a[2]]]"
                ">>;\n"
                "transform_pointwise:=func<a,b|<\n"
                    "[a[1][i]*b[1][i]:i in [1..#a[1]]],\n"
                    "[a[2][i]*b[2][i]:i in [1..#a[2]]]"
                ">>;\n"
                "Lseq_as_LP:=func<x|<Polynomial(x[1]),Polynomial(x[2])>>;\n"
                "LP_canonical:=func<x|<\n"
                    "x[1] mod (LP1.1^fft_K+1),\n"
                    "x[2] mod (LP2.1^fft_K+1)>>;\n"
                "LP_add:=func<a,b|<a[1]+b[1],a[2]+b[2]>>;\n"
                "LP_mul:=func<a,b|<a[1]*b[1],a[2]*b[2]>>;\n"
                "L_evaluate:=func<P,n|<\n"
                "[Evaluate(P[1],evalpoint1(i)) : i in [0..n-1]],\n"
                "[Evaluate(P[2],evalpoint2(i)) : i in [0..n-1]]>>;\n"
                "Zseq_to_Lseq:=func<X|<\n"
"                    [Evaluate(Zseq_to_KP(X[i..i+fft_2np1-1]),L1.1):i in [1..fft_K*fft_2np1 by fft_2np1]],\n"
"                    [Evaluate(Zseq_to_KP(X[i..i+fft_2np2-1]),L2.1):i in [fft_K*fft_2np1 + 1..#X by fft_2np2]]>"
                ">;\n"
                "Lseq_to_Zseq:=func<XX|&cat [ &cat [KP_to_Zseq(Polynomial(Eltseq(c)),Degree(Universe(X))):c in X] : X in XX]>;\n"
                , stdout);
        fputs(interpolate_code, stdout);
        fputs(
                "L_interpolate:=func<T,n|<\n"
                "interpolate_generic(T[1],[evalpoint1(i): i in [0..n-1]]),\n"
                "interpolate_generic(T[2],[evalpoint2(i): i in [0..n-1]])>>;\n"
                , stdout);
        if (doing_mp) {
            printf("Kwrap:=x^((ks1+ks2-1) * fft_K)-1;\n");
        }
    }
}

#define PRINT_CONTEXT(E) print_context_ ## E(E)

/* cpp-based static if... */
#define CROSS_gf2x_fake_fft(X) /**/
#define CROSS_gf2x_cantor_fft(X) X
#define CROSS_gf2x_ternary_fft(X) X

#define ONE_TEST(E) do {                                                \
        PRINT_CONTEXT(E);                                               \
        DO_dft_f1(E);                                                   \
        DO_dft_g1(E);                                                   \
        DO_dft_f2(E);                                                   \
        DO_dft_g2(E);                                                   \
        DO_zero_th(E);                                                  \
        DO_comp1(E);                                                    \
        DO_comp2(E);                                                    \
        DO_ift_h(E);                                                    \
        DO_zero_thx(E);                                                 \
        DO_comp_n(E);                                                   \
        DO_ift_hx(E);                                                   \
        CHECK_SELF_ADDCOMPOSE_N_CONSISTENCY(E);                         \
        CROSS_ ## E(CHECK_CROSS_CONSISTENCY(E, gf2x_fake_fft));         \
} while (0)

long randomly_pick_order_for_ternary(size_t N)
{
    size_t Kmin = 1, imin = 0;
    for( ; Kmin * Kmin * 27 < N ; Kmin *= 3, imin++) ;
    for( ; Kmin < 64 ; Kmin *= 3, imin++) ;
    size_t Kmax = Kmin, imax = imin;
    for( ; Kmax * Kmax < N * 27 ; Kmax *= 3, imax++) ;
    if (imax == imin) { imax++, Kmax *= 3; }
    int i = imin + random() % (imax - imin);
    long K = 1;
    for(int j = 0 ; j < i ; j++) K *= 3;
    if (random() & 8)
        K = -K;
    return K;
}

int docheck_mul(size_t nf, size_t ng, int nrep)
{
    // doing f1g1+f2g2
    doing_mp = 0;

    unsigned int nf1 = nf;
    unsigned int ng1 = ng;
    unsigned int nf2 = nf;
    unsigned int ng2 = ng;
    unsigned int nh = nf + ng - 1;
    unsigned int nhx = nh;

    size_t nwf1 = iceildiv(nf1, ULONG_BITS);
    size_t nwg1 = iceildiv(ng1, ULONG_BITS);
    size_t nwf2 = iceildiv(nf2, ULONG_BITS);
    size_t nwg2 = iceildiv(ng2, ULONG_BITS);

    size_t nwh = iceildiv(nh, ULONG_BITS);
    size_t nwhx = nwh;

    unsigned long * f1 = malloc(nwf1 * sizeof(unsigned long));
    unsigned long * g1 = malloc(nwg1 * sizeof(unsigned long));
    unsigned long * f2 = malloc(nwf2 * sizeof(unsigned long));
    unsigned long * g2 = malloc(nwg2 * sizeof(unsigned long));

    unsigned int seed;
    for(int i = 0 ; i < nrep ; i++) {
        seed = next_seed;
        srandom(seed);
        next_seed = random();

        long K = randomly_pick_order_for_ternary(nh);

        ENTER(gf2x_fake_fft, info_init, nf, ng);
        ENTER(gf2x_cantor_fft, info_init, nf, ng);
        ENTER(gf2x_ternary_fft, info_init, nf, ng);

        gf2x_ternary_fft_info_adjust(gf2x_ternary_fft, GF2X_FFT_ADJUST_DEPTH, K);

        SETUP(gf2x_fake_fft);
        SETUP(gf2x_cantor_fft);
        SETUP(gf2x_ternary_fft);

        CHOP_HEAD(f1);
        CHOP_HEAD(g1);
        CHOP_HEAD(f2);
        CHOP_HEAD(g2);

        ONE_TEST(gf2x_fake_fft);
        ONE_TEST(gf2x_cantor_fft);
        // ONE_TEST(gf2x_ternary_fft);
do {
    print_context_gf2x_ternary_fft(gf2x_ternary_fft);
    do {
	if (magma) {
	    do {
		display("f1", f1, nwf1);
		printf("K" "f1" ":=Zseq_to_KP(" "f1" ");\n");
		printf("L" "f1" ":=KP_to_LP(K" "f1" ");\n");
	    } while (0);
	}
	gf2x_ternary_fft_dft(gf2x_ternary_fft, gf2x_ternary_fft_tf1, f1, nf1,
			     gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "f1",
			(unsigned long *) gf2x_ternary_fft_tf1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "f1" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "f1"
		       ");\n");
	    } while (0);
	    printf("assert T" "f1" " eq L_evaluate(L" "f1"
		   ", transform_length(T" "f1" "));\n");
	}
    } while (0);
    do {
	if (magma) {
	    do {
		display("g1", g1, nwg1);
		printf("K" "g1" ":=Zseq_to_KP(" "g1" ");\n");
		printf("L" "g1" ":=KP_to_LP(K" "g1" ");\n");
	    } while (0);
	}
	gf2x_ternary_fft_dft(gf2x_ternary_fft, gf2x_ternary_fft_tg1, g1, ng1,
			     gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "g1",
			(unsigned long *) gf2x_ternary_fft_tg1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "g1" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "g1"
		       ");\n");
	    } while (0);
	    printf("assert T" "g1" " eq L_evaluate(L" "g1"
		   ", transform_length(T" "g1" "));\n");
	}
    } while (0);
    do {
	if (magma) {
	    do {
		display("f2", f2, nwf2);
		printf("K" "f2" ":=Zseq_to_KP(" "f2" ");\n");
		printf("L" "f2" ":=KP_to_LP(K" "f2" ");\n");
	    } while (0);
	}
	gf2x_ternary_fft_dft(gf2x_ternary_fft, gf2x_ternary_fft_tf2, f2, nf2,
			     gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "f2",
			(unsigned long *) gf2x_ternary_fft_tf2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "f2" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "f2"
		       ");\n");
	    } while (0);
	    printf("assert T" "f2" " eq L_evaluate(L" "f2"
		   ", transform_length(T" "f2" "));\n");
	}
    } while (0);
    do {
	if (magma) {
	    do {
		display("g2", g2, nwg2);
		printf("K" "g2" ":=Zseq_to_KP(" "g2" ");\n");
		printf("L" "g2" ":=KP_to_LP(K" "g2" ");\n");
	    } while (0);
	}
	gf2x_ternary_fft_dft(gf2x_ternary_fft, gf2x_ternary_fft_tg2, g2, ng2,
			     gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "g2",
			(unsigned long *) gf2x_ternary_fft_tg2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "g2" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "g2"
		       ");\n");
	    } while (0);
	    printf("assert T" "g2" " eq L_evaluate(L" "g2"
		   ", transform_length(T" "g2" "));\n");
	}
    } while (0);
    gf2x_ternary_fft_zero(gf2x_ternary_fft, gf2x_ternary_fft_th, 1);
    do {
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "h",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "h" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "h"
		       ");\n");
	    } while (0);
	    printf("previous_Th:=Th;\n");
	}
	gf2x_ternary_fft_addcompose(gf2x_ternary_fft, gf2x_ternary_fft_th,
				    (gf2x_ternary_fft_srcptr)
				    gf2x_ternary_fft_tf1,
				    (gf2x_ternary_fft_srcptr)
				    gf2x_ternary_fft_tg1,
				    gf2x_ternary_fft_temp2,
				    gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "h",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "h" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "h"
		       ");\n");
	    } while (0);
	    printf("assert Th eq transform_add(previous_Th,"
		   " transform_pointwise(Tf" "1" ",Tg" "1" "));\n");
	}
    } while (0);
    do {
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "h",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "h" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "h"
		       ");\n");
	    } while (0);
	    printf("previous_Th:=Th;\n");
	}
	gf2x_ternary_fft_addcompose(gf2x_ternary_fft, gf2x_ternary_fft_th,
				    (gf2x_ternary_fft_srcptr)
				    gf2x_ternary_fft_tf2,
				    (gf2x_ternary_fft_srcptr)
				    gf2x_ternary_fft_tg2,
				    gf2x_ternary_fft_temp2,
				    gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "h",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "h" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "h"
		       ");\n");
	    } while (0);
	    printf("assert Th eq transform_add(previous_Th,"
		   " transform_pointwise(Tf" "2" ",Tg" "2" "));\n");
	}
    } while (0);
    do {
	if (magma) {
	    printf("saved_T" "h" ":=T" "h" ";\n");
	}
	gf2x_ternary_fft_ift(gf2x_ternary_fft, gf2x_ternary_fft_h, nh,
			     gf2x_ternary_fft_th, gf2x_ternary_fft_temp1);
	if (magma) {
	    printf("/* " "gf2x_ternary_fft"
		   "_th has has undergone ift now */\n");
	    do {
		display("gf2x_ternary_fft" "_t" "h",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "h" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "h"
		       ");\n");
	    } while (0);
	    do {
		display("gf2x_ternary_fft" "_" "h", gf2x_ternary_fft_h, nwh);
		printf("K" "h" ":=Zseq_to_KP(" "gf2x_ternary_fft" "_" "h"
		       ");\n");
	    } while (0);
	    printf("interpolated_" "h" ":=Lseq_as_LP(T" "h" ");\n");
	    printf("assert L_evaluate(interpolated_" "h" ","
		   "transform_length(T" "h" ")) eq " "saved_T" "h" ";\n");
	    printf("assert interpolated_" "h" " eq " "L_interpolate(saved_T"
		   "h" ", " "transform_length(T" "h" "));\n");
	    printf("assert interpolated_" "h" " eq "
		   "LP_canonical(LP_add(LP_mul(Lf1, Lg1),"
		   " LP_mul(Lf2, Lg2)));\n");
	    if (doing_mp) {
		printf("if assigned Kwrap then\n"
		       "\tassert (LP_to_KP(interpolated_" "h" ") - "
		       "(Kf1*Kg1+Kf2*Kg2)) mod Kwrap eq 0;\n" "else\n"
		       "\tassert LP_to_KP(interpolated_" "h" ") eq "
		       "Kf1*Kg1+Kf2*Kg2;\n" "end if;\n");
		printf("assert LP_to_KP(interpolated_" "h" ") "
		       "mod x^%zu div x^%zu " "eq Zseq_to_KP("
		       "gf2x_ternary_fft" "_h);\n", MAX(nf, ng), MIN(nf,
								     ng) - 1);
	    } else {
		printf("assert LP_to_KP(interpolated_" "h" ") eq "
		       "Kf1*Kg1+Kf2*Kg2;\n");
	    }
	}
    } while (0);
    gf2x_ternary_fft_zero(gf2x_ternary_fft, gf2x_ternary_fft_thx, 1);
    do {
	gf2x_ternary_fft_addcompose_n(gf2x_ternary_fft, gf2x_ternary_fft_thx,
				      (gf2x_ternary_fft_srcptr *)
				      gf2x_ternary_fft_tfs,
				      (gf2x_ternary_fft_srcptr *)
				      gf2x_ternary_fft_tgs, 2,
				      gf2x_ternary_fft_temp2,
				      gf2x_ternary_fft_temp1);
	if (magma) {
	    do {
		display("gf2x_ternary_fft" "_t" "hx",
			(unsigned long *) gf2x_ternary_fft_thx,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "hx" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "hx"
		       ");\n");
	    } while (0);
	    printf("assert Thx eq transform_add("
		   "transform_pointwise(Tf1,Tg1),"
		   "transform_pointwise(Tf2,Tg2));\n");
	}
    } while (0);
    do {
	if (magma) {
	    printf("saved_T" "hx" ":=T" "hx" ";\n");
	}
	gf2x_ternary_fft_ift(gf2x_ternary_fft, gf2x_ternary_fft_hx, nhx,
			     gf2x_ternary_fft_thx, gf2x_ternary_fft_temp1);
	if (magma) {
	    printf("/* " "gf2x_ternary_fft"
		   "_th has has undergone ift now */\n");
	    do {
		display("gf2x_ternary_fft" "_t" "hx",
			(unsigned long *) gf2x_ternary_fft_thx,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		printf("T" "hx" ":=Zseq_to_Lseq(" "gf2x_ternary_fft" "_t" "hx"
		       ");\n");
	    } while (0);
	    do {
		display("gf2x_ternary_fft" "_" "hx", gf2x_ternary_fft_hx,
			nwhx);
		printf("K" "hx" ":=Zseq_to_KP(" "gf2x_ternary_fft" "_" "hx"
		       ");\n");
	    } while (0);
	    printf("interpolated_" "hx" ":=Lseq_as_LP(T" "hx" ");\n");
	    printf("assert L_evaluate(interpolated_" "hx" ","
		   "transform_length(T" "hx" ")) eq " "saved_T" "hx" ";\n");
	    printf("assert interpolated_" "hx" " eq " "L_interpolate(saved_T"
		   "hx" ", " "transform_length(T" "hx" "));\n");
	    printf("assert interpolated_" "hx" " eq "
		   "LP_canonical(LP_add(LP_mul(Lf1, Lg1),"
		   " LP_mul(Lf2, Lg2)));\n");
	    if (doing_mp) {
		printf("if assigned Kwrap then\n"
		       "\tassert (LP_to_KP(interpolated_" "hx" ") - "
		       "(Kf1*Kg1+Kf2*Kg2)) mod Kwrap eq 0;\n" "else\n"
		       "\tassert LP_to_KP(interpolated_" "hx" ") eq "
		       "Kf1*Kg1+Kf2*Kg2;\n" "end if;\n");
		printf("assert LP_to_KP(interpolated_" "hx" ") "
		       "mod x^%zu div x^%zu " "eq Zseq_to_KP("
		       "gf2x_ternary_fft" "_h);\n", MAX(nf, ng), MIN(nf,
								     ng) - 1);
	    } else {
		printf("assert LP_to_KP(interpolated_" "hx" ") eq "
		       "Kf1*Kg1+Kf2*Kg2;\n");
	    }
	}
    } while (0);
    do {
	if (memcmp
	    (gf2x_ternary_fft_h, gf2x_ternary_fft_hx,
	     nwh * sizeof(unsigned long)) != 0) {
	    fprintf(stderr,
		    "gf2x_ternary_fft" " != " "gf2x_ternary_fft" "x "
		    "(addcompose_n consistency) for "
		    "-%s -a %zu -b %zu -seed %u\n", __func__ + 8, nf, ng,
		    seed);
	    printf("w:=%d;\n", (int) ULONG_BITS);
	    display("f1", f1, nwf1);
	    display("g1", g1, nwg1);
	    display("f2", f2, nwf2);
	    display("g2", g2, nwg2);
	    display("gf2x_ternary_fft" "_h", gf2x_ternary_fft_h, nwh);
	    display("gf2x_ternary_fft" "_hx", gf2x_ternary_fft_hx, nwh);
	    do {
		display("gf2x_ternary_fft" "_tf1",
			(unsigned long *) gf2x_ternary_fft_tf1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tg1",
			(unsigned long *) gf2x_ternary_fft_tg1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tf2",
			(unsigned long *) gf2x_ternary_fft_tf2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tg2",
			(unsigned long *) gf2x_ternary_fft_tg2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
	    } while (0);
	    do {
		display("gf2x_ternary_fft" "_th",
			(unsigned long *) gf2x_ternary_fft_th,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_thx",
			(unsigned long *) gf2x_ternary_fft_thx,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
	    } while (0);
	    fflush(stdout);
	    fflush(stderr);
	    abort();
	}
    } while (0);
    do {
	if (memcmp
	    (gf2x_fake_fft_h, gf2x_ternary_fft_h,
	     nwh * sizeof(unsigned long)) != 0) {
	    fprintf(stderr,
		    "gf2x_fake_fft" " != " "gf2x_ternary_fft" " "
		    "(cross check) for " "-%s -a %zu -b %zu -seed %u\n",
		    __func__ + 8, nf, ng, seed);
	    printf("w:=%d;\n", (int) ULONG_BITS);
	    display("f1", f1, nwf1);
	    display("g1", g1, nwg1);
	    display("f2", f2, nwf2);
	    display("g2", g2, nwg2);
	    display("gf2x_fake_fft" "_h", gf2x_fake_fft_h, nwh);
	    display("gf2x_ternary_fft" "_h", gf2x_ternary_fft_h, nwh);
	    do {
		display("gf2x_fake_fft" "_tf1",
			(unsigned long *) gf2x_fake_fft_tf1,
			gf2x_fake_fft_transform_size(gf2x_fake_fft) *
			sizeof(gf2x_fake_fft_elt) / sizeof(unsigned long));
		display("gf2x_fake_fft" "_tg1",
			(unsigned long *) gf2x_fake_fft_tg1,
			gf2x_fake_fft_transform_size(gf2x_fake_fft) *
			sizeof(gf2x_fake_fft_elt) / sizeof(unsigned long));
		display("gf2x_fake_fft" "_tf2",
			(unsigned long *) gf2x_fake_fft_tf2,
			gf2x_fake_fft_transform_size(gf2x_fake_fft) *
			sizeof(gf2x_fake_fft_elt) / sizeof(unsigned long));
		display("gf2x_fake_fft" "_tg2",
			(unsigned long *) gf2x_fake_fft_tg2,
			gf2x_fake_fft_transform_size(gf2x_fake_fft) *
			sizeof(gf2x_fake_fft_elt) / sizeof(unsigned long));
	    } while (0);
	    do {
		display("gf2x_ternary_fft" "_tf1",
			(unsigned long *) gf2x_ternary_fft_tf1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tg1",
			(unsigned long *) gf2x_ternary_fft_tg1,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tf2",
			(unsigned long *) gf2x_ternary_fft_tf2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
		display("gf2x_ternary_fft" "_tg2",
			(unsigned long *) gf2x_ternary_fft_tg2,
			gf2x_ternary_fft_transform_size(gf2x_ternary_fft) *
			sizeof(gf2x_ternary_fft_elt) / sizeof(unsigned long));
	    } while (0);
	    fflush(stdout);
	    fflush(stderr);
	    abort();
	}
    } while (0);
} while (0);

        LEAVE(gf2x_fake_fft);
        LEAVE(gf2x_cantor_fft);
        LEAVE(gf2x_ternary_fft);
    }

    free(g2);
    free(f2);
    free(g1);
    free(f1);

    return 0;
}


int docheck_mp(size_t nf, size_t ng, int nrep)
{
    // doing MP(f1,g1)+MP(f2,g2)
    doing_mp = 1;

    unsigned int nf1 = nf;
    unsigned int ng1 = ng;
    unsigned int nf2 = nf;
    unsigned int ng2 = ng;
    unsigned int nh = MAX(nf, ng) - MIN(nf, ng) + 1;
    unsigned int nhx = nh;

    size_t nwf1 = iceildiv(nf1, ULONG_BITS);
    size_t nwg1 = iceildiv(ng1, ULONG_BITS);
    size_t nwf2 = iceildiv(nf2, ULONG_BITS);
    size_t nwg2 = iceildiv(ng2, ULONG_BITS);

    size_t nwh = iceildiv(nh, ULONG_BITS);
    size_t nwhx = nwh;

    unsigned long * f1 = malloc(nwf1 * sizeof(unsigned long));
    unsigned long * g1 = malloc(nwg1 * sizeof(unsigned long));
    unsigned long * f2 = malloc(nwf2 * sizeof(unsigned long));
    unsigned long * g2 = malloc(nwg2 * sizeof(unsigned long));

    unsigned int seed;
    for(int i = 0 ; i < nrep ; i++) {
        seed = next_seed;
        srandom(seed);
        next_seed = random();

        long K = randomly_pick_order_for_ternary(nh);

        ENTER(gf2x_fake_fft, info_init_mp, nf, ng);
        ENTER(gf2x_cantor_fft, info_init_mp, nf, ng);
        ENTER(gf2x_ternary_fft, info_init_mp, nf, ng);

        gf2x_ternary_fft_info_adjust(gf2x_ternary_fft, GF2X_FFT_ADJUST_DEPTH, K);

        SETUP(gf2x_fake_fft);
        SETUP(gf2x_cantor_fft);
        SETUP(gf2x_ternary_fft);


        CHOP_HEAD(f1);
        CHOP_HEAD(g1);
        CHOP_HEAD(f2);
        CHOP_HEAD(g2);

        ONE_TEST(gf2x_fake_fft);
        ONE_TEST(gf2x_cantor_fft);
        ONE_TEST(gf2x_ternary_fft);

        LEAVE(gf2x_fake_fft);
        LEAVE(gf2x_cantor_fft);
        LEAVE(gf2x_ternary_fft);
    }

    free(g2);
    free(f2);
    free(g1);
    free(f1);

    return 0;
}


int main(int argc, char *argv[])
{
    unsigned int max_nbits_a = UINT_MAX;
    unsigned int max_nbits_b = UINT_MAX;
    unsigned long seed = ULONG_MAX;
    int nrep = 10;
    int mp = 1, mul = 1;
    argc--, argv++;
    for (; argc; argc--, argv++) {
        if (strcmp(argv[0], "-a") == 0) {
            argc--, argv++;
            max_nbits_a = atol(argv[0]);
        } else if (strcmp(argv[0], "-b") == 0) {
            argc--, argv++;
            max_nbits_b = atol(argv[0]);
        } else if (strcmp(argv[0], "-seed") == 0) {
            argc--, argv++;
            seed = atol(argv[0]);
        } else if (strcmp(argv[0], "-nr") == 0) {
            argc--, argv++;
            nrep = atol(argv[0]);
        } else if (strcmp(argv[0], "-m") == 0) {
            magma = 1;
        } else if (strcmp(argv[0], "-mp") == 0) {
            mp = 1; mul = 0;
        } else if (strcmp(argv[0], "-mul") == 0) {
            mul = 1; mp = 0;
        } else {
            usage();
        }
    }
    if (max_nbits_a == UINT_MAX)
        max_nbits_a = 1000;
    if (max_nbits_b == UINT_MAX)
        max_nbits_b = max_nbits_a;
    if (seed == ULONG_MAX) {
        /* make this completely deterministic in debug mode */
#ifdef NDEBUG
        srandom(time(NULL));
#ifdef  _POSIX_C_SOURCE
        for (int i = 0; i < getpid() % 1009; i++)
            random();
#endif
#endif
        seed = random();
    }

    next_seed = seed;
    printf("// Replicate full run with %s%s-a %u -b %u -seed %lu -nr %d\n",
            (mul && !mp) ? "-mul " : "",
            (mp && !mul) ? "-mp " : "",
           max_nbits_a, max_nbits_b, seed, nrep);

    if (mul) {
        docheck_mul(max_nbits_a, max_nbits_b, nrep);

        for (unsigned int i = 0; i < 200; i++) {
            unsigned int na = 100 + (random() % max_nbits_a);
            unsigned int nb = 100 + (random() % max_nbits_b);
            docheck_mul(na, nb, nrep);
            fputc('.', stderr);
            fflush(stderr);
        }
        printf("\n");
        for (unsigned int i = 0; i < max_nbits_a;) {
            // test some fancy sizes as well.
            i += 2 * ULONG_BITS / 3 + (random() % 3);
            docheck_mul(i, i * max_nbits_b / max_nbits_a, nrep);
            fputc('.', stderr);
            fflush(stderr);
        }
        printf("\n");
    }
    if (mp) {
        docheck_mp(max_nbits_a, max_nbits_b, nrep);

        for (unsigned int i = 0; i < 200; i++) {
            unsigned int na = 100 + (random() % max_nbits_a);
            unsigned int nb = 100 + (random() % max_nbits_b);
            docheck_mp(na, nb, nrep);
            fputc('.', stderr);
            fflush(stderr);
        }
        printf("\n");
        for (unsigned int i = 0; i < max_nbits_a;) {
            // test some fancy sizes as well.
            i += 2 * ULONG_BITS / 3 + (random() % 3);
            docheck_mp(i, i * max_nbits_b / max_nbits_a, nrep);
            fputc('.', stderr);
            fflush(stderr);
        }
        printf("\n");
    }
    return 0;
}
