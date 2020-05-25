#include "cado.h" // IWYU pragma: keep

#include <ostream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <limits.h>

#include "lingen_mat_types.hpp"
#include "lingen_qcode.h"
#include "macros.h"

using namespace std;

/* if SWAP_E is defined, put E[i,j] into E[j][i] */
#define SWAP_E

/* if SWAP_PI is defined, put pi[i,j] into pi[j][i] */
#define SWAP_PI

template<int WIDTH>
class ulmat_rowmajor {
    unsigned int n;
    unsigned long (*x)[WIDTH];
    ulmat_rowmajor(ulmat_rowmajor const&) {}
public:
    ulmat_rowmajor(unsigned int m, unsigned int n) : n(n)
    {
        x = new unsigned long[m*n][WIDTH];
        memset(x, 0, m*n*sizeof(unsigned long[WIDTH]));
    }
    ~ulmat_rowmajor() { delete[] x; }
    inline unsigned long (&operator()(unsigned int i, unsigned int j))[WIDTH] {
        return x[i * n + j];
    }
    inline const unsigned long (&operator()(unsigned int i, unsigned int j)const)[WIDTH] {
        return x[i * n + j];
    }
    inline unsigned long (*row(unsigned int i))[WIDTH] { return x + i * n; }
    inline const unsigned long (*row(unsigned int i)const)[WIDTH] { return x + i * n; }
};

template<int WIDTH>
class ulmat_colmajor {
    unsigned int m;
    unsigned long (*x)[WIDTH];
    ulmat_colmajor(ulmat_colmajor const&) {}
public:
    ulmat_colmajor(unsigned int m, unsigned int n) : m(m)
    {
        x = new unsigned long[m*n][WIDTH];
        memset(x, 0, m*n*sizeof(unsigned long[WIDTH]));
    }
    ~ulmat_colmajor() { delete[] x; }
    inline unsigned long (&operator()(unsigned int i, unsigned int j))[WIDTH] {
        return x[j * m + i];
    }
    inline const unsigned long (&operator()(unsigned int i, unsigned int j)const)[WIDTH] {
        return x[j * m + i];
    }
    inline unsigned long (*column(unsigned int j))[WIDTH] { return x + j * m; }
    inline const unsigned long (*column(unsigned int j)const)[WIDTH] { return x + j * m; }
};


bool report_spontaneous_zeros(lingen_qcode_data_ptr qq, unsigned int dt, vector<int>const& is_modified)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;
    vector<pair<unsigned int, unsigned int> > magic;
    bool changed = false;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        if (is_modified[i]) {
            qq->ch[i] = 0;
            if (is_modified[i] & 1)
                changed=true;
        } else {
            qq->ch[i]++;
            magic.push_back(make_pair(i, qq->ch[i]));
        }
    }

    if (magic.empty()) {
        ASSERT_ALWAYS(changed);
        return changed;
    }

    printf("%-8u%zucols=0:", qq->t + dt, magic.size());

    // Now print this out more nicely.
    sort(magic.begin(),magic.end());
    for( ; magic.size() ; ) {
        unsigned int mi = UINT_MAX;
        for(unsigned int i = 0 ; i < magic.size() ; i++) {
            if (magic[i].second < mi) {
                mi = magic[i].second;
            }
        }
        printf(" [");
        for(unsigned int i = 0 ; i < magic.size() ; ) {
            unsigned int j;
            for(j = i; j < magic.size() ; j++) {
                if (magic[j].first-magic[i].first != j-i) break;
            }
            if (i) printf(",");
            if (magic[i].first == magic[j-1].first - 1) {
                printf("%u,%u", magic[i].first, magic[j-1].first);
            } else if (magic[i].first < magic[j-1].first) {
                printf("%u..%u", magic[i].first, magic[j-1].first);
            } else {
                printf("%u", magic[i].first);
            }
            i = j;
        }
        printf("]");
        if (mi > 1)
            printf("*%u",mi);

        vector<pair<unsigned int, unsigned int> > zz2;
        for(unsigned int i = 0 ; i < magic.size() ; i++) {
            if (magic[i].second > mi) {
                zz2.push_back(make_pair(magic[i].first,magic[i].second));
            }
        }
        magic.swap(zz2);
    }
    printf("\n");

    return changed;
}

inline void lshift1(unsigned long&x) { x <<= 1; }

template<int WIDTH> inline void lshift1(unsigned long (&x)[WIDTH])
{
    mpn_lshift(x, x, WIDTH, 1);
}

template<> inline void lshift1<1>(unsigned long (&x)[1]) { x[0] <<= 1; }
template<> inline void lshift1<2>(unsigned long (&x)[2]) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    asm(
            "addq %0,%0\n"
            "adcq %1,%1\n"
            : "=r"(x[0]), "=r"(x[1])
            : "0" (x[0]), "1"(x[1])
            :);
#else
    x[1] <<= 1;
    x[1] |= ((long)x[0]) < 0;
    x[0] <<= 1;
#endif
}

template<int WIDTH>
unsigned int lingen_qcode_do_tmpl(lingen_qcode_data_ptr qq)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;

#ifdef  SWAP_E
    ulmat_colmajor<WIDTH> E(m, b);
#else
    ulmat_rowmajor<WIDTH> E(m, b);
#endif

#ifdef  SWAP_PI
    ulmat_colmajor<WIDTH> P(b, b);
#else
    ulmat_rowmajor<WIDTH> P(b, b);
#endif

    ASSERT_ALWAYS(qq->length <= WIDTH * ULONG_BITS);
    for (unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            memcpy(E(i, j), qq->iptrs[i * b + j], WIDTH * sizeof(unsigned long));
        }
    }

    for (unsigned int i = 0; i < b; i++) {
        P(i, i)[0] = 1;
        qq->local_delta[i] = 0;
    }

    /* We need to keep track of the qq->ch array. It checks how many
     * times in a row a given column turns out to be magically zero. 
     *
     * For this, we'll use a local table (is_modified[] below), and we
     * will use it to report which columns happen to be zero, and update
     * the qq->ch array.
     */

    unsigned int e = 0;
    for ( ; e < qq->length; e++) {
	uint64_t mask = (uint64_t) 1 << (e % ULONG_BITS);
	int mpos = e / ULONG_BITS;
        vector<int> is_modified(m + n, false);
#if 0
	unsigned int md = UINT_MAX;
	for (unsigned int j = 0; j < m + n; j++)
	    if (qq->delta[j] > md)
		md = qq->delta[j];
#endif
	for (unsigned int i = 0; i < m; i++) {
	    unsigned int min_degree = UINT_MAX;
            unsigned int pivot = m + n;
	    for (unsigned int j = 0; j < m + n; j++)
		if ((E(i, j)[mpos] & mask) && (qq->delta[j] < min_degree)) {
		    min_degree = qq->delta[j];
		    pivot = j;
		}
	    ASSERT_ALWAYS(pivot < m + n);
            /* A pivot column does not count as undergoing a modification
             * the same way as the other columns. This is because for our
             * stop criterion, what we need is to count the number of
             * transvection matrices by which we multiply. However, as
             * far as detection of spontaneous zeros goes, of course we
             * must make sure that we don't take pivots into account ! */
            is_modified[pivot] |= 2;
	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k)[mpos] & mask)) continue;
                is_modified[k] |= 1;
#if defined(SWAP_E)
                mpn_xor_n ((unsigned long*)E.column(k),(unsigned long*) E.column(k),(unsigned long*) E.column(pivot), m * WIDTH);
#else
                for (unsigned int l = 0; l < m; l++)
                    for (unsigned int s = 0; s < WIDTH; s++)
                        E(l, k)[s] ^= E(l, pivot)[s];
#endif
#if defined(SWAP_PI)
                mpn_xor_n ((unsigned long*)P.column(k), (unsigned long*)P.column(k), (unsigned long*)P.column(pivot), (m + n) * WIDTH);
#else
                for (unsigned int l = 0; l < m + n; l++)
                    for (unsigned int s = 0; s < WIDTH; s++)
                        P(l, k)[s] ^= P(l, pivot)[s];
#endif
                if (qq->local_delta[pivot] > qq->local_delta[k])
                    qq->local_delta[k] = qq->local_delta[pivot];
	    }
	    for (unsigned int l = 0; l < m; l++)
                lshift1(E(l, pivot));
	    for (unsigned int l = 0; l < m + n; l++)
                lshift1(P(l, pivot));
	    qq->delta[pivot] += 1;
            qq->local_delta[pivot]++;
	}
        /* Columns which have not been changed here are those which
         * are magically zero. Count whether this happens often or
         * not.
         *
         * While we're doing this, we can also detect whether our
         * computation stays still. */
        if (!report_spontaneous_zeros(qq, e, is_modified))
            break;
    }

    for (unsigned int i = 0; i < b; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            memcpy(qq->optrs[i * b + j], P(i, j), WIDTH * sizeof(unsigned long));
        }
    }
    qq->t += e;
    return e;
}

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    if (qq->length <= ULONG_BITS) {
        return lingen_qcode_do_tmpl<1>(qq);
    } else if (qq->length <= 2 * ULONG_BITS) {
        return lingen_qcode_do_tmpl<2>(qq);
    } else if (qq->length <= 3 * ULONG_BITS) {
        return lingen_qcode_do_tmpl<3>(qq);
    } else if (qq->length <= 4 * ULONG_BITS) {
        return lingen_qcode_do_tmpl<4>(qq);
    }
    cerr << "quadratic algorithm not supported for base length " << qq->length << endl;
    abort();
    return 0;
}
