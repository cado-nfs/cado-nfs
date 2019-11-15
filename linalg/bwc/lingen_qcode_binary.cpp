#include "cado.h"

#include <algorithm>
#include <vector>
#include <type_traits>

#define LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE
#include "lingen_qcode_binary.hpp"
#include "utils.h"

#include "lingen_expected_pi_length.hpp"

void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength)
{
    qq->m = m;
    qq->b = b;
    qq->length = length;
    qq->outlength = outlength;
    qq->iptrs = (const unsigned long **) malloc(m * b * sizeof(unsigned long *));
    qq->optrs = (unsigned long **) malloc(b * b * sizeof(unsigned long *));
    memset(qq->iptrs, 0, m * b * sizeof(unsigned long *));
    memset(qq->optrs, 0, b * b * sizeof(unsigned long *));

    qq->local_delta = (unsigned int *) malloc(b * sizeof(unsigned int));
    memset(qq->local_delta, 0, b * sizeof(unsigned int));
}
void lingen_qcode_clear(lingen_qcode_data_ptr qq)
{
    free(qq->iptrs);
    free(qq->optrs);
    free(qq->local_delta);
    memset(qq, 0, sizeof(*qq));
}

#if 0
unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    /* It's not technically necessary, but the typical use case is really
     * like this, so I'm taking the liberty to shorten the code a bit
     * using this fact */
    ASSERT_ALWAYS(qq->length <= ULONG_BITS);

    /* read the complete input */
    {
        unsigned long tmask = 1;
        for(unsigned int k = 0 ; k < qq->length ; k++, tmask <<= 1) {
            unsigned long jmask = 0;
            for(unsigned int j = 0 ; j < qq->b ; j++, jmask <<= 1) {
                unsigned long * Ac = qq->A[k] + j / ULONG_BITS;
                if (!jmask) jmask = 1;
                for(unsigned int i = 0 ; i < qq->m ; i++) {
                    if (qq->iptrs[i * qq->b + j][0] & tmask) {
                        *Ac |= jmask;
                    }
                    Ac += iceildiv(qq->b, ULONG_BITS);
                }
            }
        }
    }

    /* please implement me ! */
    abort();


    /* Now, store the complete output */
    {
        unsigned long tmask = 1;
        for(unsigned int k = 0 ; k < qq->outlength ; k++, tmask <<= 1) {
            unsigned long jmask = 0;
            for(unsigned int j = 0 ; j < qq->b ; j++, jmask <<= 1) {
                unsigned long * Xc = qq->X[k] + j / ULONG_BITS;
                if (!jmask) jmask = 1;
                for(unsigned int i = 0 ; i < qq->b ; i++) {
                    if (*Xc & jmask) {
                        qq->optrs[i * qq->b + j][0] |= tmask;
                    }
                    Xc += iceildiv(qq->b, ULONG_BITS);
                }
            }
        }
    }
    return qq->length;
}
#endif

unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j)
{
    return qq->local_delta[j] + 1;
}

/* if SWAP_E is defined, put E[i,j] into E[j][i] */
#define SWAP_E

/* if SWAP_PI is defined, put pi[i,j] into pi[j][i] */
#define SWAP_PI

template<int WIDTH>
struct constant_width {
    static constexpr const size_t width = WIDTH;
    constant_width(size_t = WIDTH) {}
    constant_width(constant_width const&) = default;
    constant_width& operator=(constant_width const&) = default;
};

struct variable_width {
    size_t width;
    variable_width(size_t x = 0) : width(x) {}
    variable_width(variable_width const&) = default;
    variable_width& operator=(variable_width const&) = default;
};

template<typename width_type, typename pointer_type>
class bitarray : private width_type {
    static_assert(std::is_same<pointer_type, unsigned long *>::value ||
            std::is_same<pointer_type, const unsigned long *>::value,
            "works only for ulong* or const ulong*");
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    typename std::enable_if<std::is_same<pointer_type, unsigned long *>::value, bitarray>::type
    operator^=(bitarray const & a) {
        mpn_xor_n (x, x, a.x, width);
        return *this;
    }
    typename std::enable_if<std::is_same<pointer_type, unsigned long *>::value, void>::type
    lshift1() {
        mpn_lshift(x, x, width, 1);
    }
    typename std::enable_if<std::is_same<pointer_type, unsigned long *>::value, void>::type
    copy_from(const unsigned long * a) {
        memcpy(x, a, width * sizeof(unsigned long));
    }
    void copy_to(unsigned long * a) const {
        memcpy(a, x, width * sizeof(unsigned long));
    }
    inline bool operator[](size_t k) const {
        return x[k / ULONG_BITS] & (1UL << (k % ULONG_BITS));
    }
    typename std::enable_if<std::is_same<pointer_type, unsigned long *>::value, void>::type
    set1() { *x = 1; }
};

template<>
class bitarray<constant_width<1>, unsigned long *> : private constant_width<1> {
    typedef constant_width<1> width_type;
    typedef unsigned long * pointer_type;
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    inline bitarray operator^=(bitarray const & a) { *x ^= *a.x; return *this; }
    inline void lshift1() { *x <<= 1; }
    inline void copy_from(const unsigned long * a) { *x = *a; }
    inline void copy_to(unsigned long * a) const { *a = *x; }
    inline bool operator[](size_t k) const {
        return *x & (1UL << k);
    }
    void set1() { *x = 1; }
};

template<>
class bitarray<constant_width<1>, const unsigned long *> : private constant_width<1> {
    typedef constant_width<1> width_type;
    typedef const unsigned long * pointer_type;
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    inline void copy_to(unsigned long * a) const { *a = *x; }
    inline bool operator[](size_t k) const {
        return *x & (1UL << k);
    }
};

template<typename width_type>
class ulmat_rowmajor : public width_type {
    using width_type::width;
    unsigned int n;
    unsigned long *x;
    ulmat_rowmajor(ulmat_rowmajor const&) {}
public:
    ulmat_rowmajor(width_type w, unsigned int m, unsigned int n) : width_type(w), n(n)
    {
        x = new unsigned long[m*n*width];
        memset(x, 0, m*n*width*sizeof(unsigned long));
    }
    ~ulmat_rowmajor() { delete[] x; }
    inline bitarray<width_type, unsigned long *> operator()(unsigned int i, unsigned int j) {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (i * n + j) * width);
    }
    inline bitarray<width_type, const unsigned long *> operator()(unsigned int i, unsigned int j) const {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (i * n + j) * width);
    }
    // inline unsigned long (*row(unsigned int i))[WIDTH] { return x + i * n; }
    // inline const unsigned long (*row(unsigned int i)const)[WIDTH] { return x + i * n; }
};

template<typename width_type>
class ulmat_colmajor : public width_type {
    using width_type::width;
    unsigned int m;
    unsigned long *x;
    ulmat_colmajor(ulmat_colmajor const&) {}
public:
    ulmat_colmajor(width_type w, unsigned int m, unsigned int n) : width_type(w), m(m)
    {
        x = new unsigned long[m*n*width];
        memset(x, 0, m*n*width*sizeof(unsigned long));
    }
    ~ulmat_colmajor() { delete[] x; }
    inline bitarray<width_type, unsigned long *> operator()(unsigned int i, unsigned int j) {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (j * m + i) * width);
    }
    inline bitarray<width_type, const unsigned long *> operator()(unsigned int i, unsigned int j) const {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (j * m + i) * width);
    }
    inline bitarray<variable_width, unsigned long *> column(unsigned int j) {
        return bitarray<variable_width, unsigned long *>(m * width, x + (j * m) * width);
    }
    inline const unsigned long * column(unsigned int j) const { return (*this)(0, j); }
};

#if 0
bool report_spontaneous_zeros(lingen_qcode_data_ptr qq, unsigned int dt, std::vector<int>const& is_modified)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;
    std::vector<std::pair<unsigned int, unsigned int> > magic;
    bool changed = false;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        if (is_modified[i]) {
            qq->ch[i] = 0;
            if (is_modified[i] & 1)
                changed=true;
        } else {
            qq->ch[i]++;
            magic.push_back( { i, qq->ch[i] } );
        }
    }

    if (magic.empty()) {
        ASSERT_ALWAYS(changed);
        return changed;
    }

    printf("%-8u%zucols=0:", qq->t + dt, magic.size());

    // Now print this out more nicely.
    std::sort(magic.begin(),magic.end());
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

        std::vector<std::pair<unsigned int, unsigned int> > zz2;
        for(unsigned int i = 0 ; i < magic.size() ; i++) {
            if (magic[i].second > mi) {
                zz2.push_back( { magic[i].first,magic[i].second } );
            }
        }
        magic.swap(zz2);
    }
    printf("\n");

    return changed;
}
#endif

/* This re-implementation is inspired from the one in
 * lingen_qcode_prime.cpp */
bool generator_found(lingen_qcode_data_ptr qq, unsigned int dt, std::vector<int>const& is_modified)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;
    unsigned int newluck = 0;
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (qq->ch[j] < 0)
            continue;
        else if (!is_modified[j]) {
            newluck++, qq->ch[j]++;
        } else {
            qq->ch[j] = 0;
        }
    }

    int happy = 0;

    if (newluck) {
        /* If newluck == n, then we probably have a generator. We add an
         * extra guarantee. newluck==n, for a total of k iterations in a
         * row, means m*n*k coefficients cancelling magically. We would
         * like this to be impossible by mere chance. Thus we want n*k >
         * luck_mini, which can easily be checked */

        unsigned int luck_mini = qq->luck_mini;
        unsigned int luck_sure = 0;

        printf("t=%d, canceled columns:", qq->t + dt);
        unsigned int last = UINT_MAX;
        unsigned int afterlast = UINT_MAX;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (qq->ch[j] > 0) {
                if (j == afterlast) {
                    afterlast++;
                } else {
                    if (last != UINT_MAX) {
                        if (last == afterlast-1) {
                            printf(" %u", last);
                        } else {
                            printf(" %u-%u", last, afterlast-1);
                        }
                    }
                    last = j;
                    afterlast = j + 1;
                }
                luck_sure += ((unsigned int) qq->ch[j]) >= luck_mini;
            }
        }
        if (last != UINT_MAX) {
            if (last == afterlast-1) {
                printf(" %u", last);
            } else {
                printf(" %u-%u", last, afterlast-1);
            }
        }

        if (newluck == n && luck_sure == n) {
            if (!happy) {
                printf(", complete generator found, for sure");
            }
            happy = 1;
        }
        printf(".\n");
    }
    return happy;
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

template<typename width_type>
unsigned int lingen_qcode_do_tmpl(width_type w, lingen_qcode_data_ptr qq)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;

    int width = w.width;

#ifdef  SWAP_E
    ulmat_colmajor<width_type> E(w, m, b);
#else
    ulmat_rowmajor<width_type> E(w, m, b);
#endif

#ifdef  SWAP_PI
    ulmat_colmajor<width_type> P(w, b, b);
#else
    ulmat_rowmajor<width_type> P(w, b, b);
#endif

    ASSERT_ALWAYS(qq->length <= (unsigned long) width * ULONG_BITS);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            E(i, j).copy_from(qq->iptrs[i * b + j]);
        }
    }

    for (unsigned int i = 0; i < b; i++) {
        P(i, i).set1();
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
        std::vector<int> is_modified(m + n, false);
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
		if (E(i, j)[e] && (qq->delta[j] < min_degree)) {
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k)[e])) continue;
                is_modified[k] |= 1;
#if defined(SWAP_E)
                E.column(k) ^= E.column(pivot);
#else
                for (unsigned int l = 0; l < m; l++)
                    E(l, k) ^= E(l, pivot);
#endif
#if defined(SWAP_PI)
                P.column(k) ^= P.column(pivot);
#else
                for (unsigned int l = 0; l < m + n; l++)
                        P(l, k) ^= P(l, pivot);
#endif
                if (qq->local_delta[pivot] > qq->local_delta[k])
                    qq->local_delta[k] = qq->local_delta[pivot];
	    }
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (unsigned int l = 0; l < m; l++)
                E(l, pivot).lshift1();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (unsigned int l = 0; l < m + n; l++)
                P(l, pivot).lshift1();
	    qq->delta[pivot] += 1;
            qq->local_delta[pivot]++;
	}
        /* Columns which have not been changed here are those which
         * are magically zero. Count whether this happens often or
         * not.
         *
         * While we're doing this, we can also detect whether our
         * computation stays still. */
        if (generator_found(qq, e, is_modified))
            break;
    }

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (unsigned int i = 0; i < b; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            P(i, j).copy_to(qq->optrs[i * b + j]);
        }
    }
    qq->t += e;
    return e;
}

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    if (qq->length <= ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<1>(), qq);
    } else if (qq->length <= 2 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<2>(), qq);
    } else if (qq->length <= 3 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<3>(), qq);
    } else if (qq->length <= 4 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<4>(), qq);
    } else if (qq->length <= 5 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<5>(), qq);
    } else if (qq->length <= 6 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<6>(), qq);
    } else if (qq->length <= 7 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<7>(), qq);
    } else if (qq->length <= 8 * ULONG_BITS) {
        return lingen_qcode_do_tmpl(constant_width<8>(), qq);
    } else {
        return lingen_qcode_do_tmpl(variable_width { iceildiv(qq->length, ULONG_BITS) }, qq);
    }
    return 0;
}

int bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly & E)/*{{{*/
{
    /* There's a nasty bug. Revealed by 32-bits, but can occur on larger
     * sizes too. Let W be the word size. When E has length W + epsilon,
     * pi_deg_bound(W+epsilon) may be < W. So that we may have tmp_pi and
     * E have stride 1 and 2 respectively. However, in
     * lingen_qcode_do_tmpl, we fill the pointers in tmp_pi (optrs[]
     * there) using the stride of E. This is not proper. A kludge is to
     * make them compatible.
     *
     * (A priori we're not talking immense sizes at the recursive threshold anyway).
     */

    /* (expected_pi_length should do as well, but 1+E.get_size() is
     * firmly on the safe side !. Note that for E.get_size() == 0, we do
     * need length 2.
     */
    size_t exp_maxlen = 1 + E.get_size();

    matpoly tmp_pi(bm.d.ab, E.ncols(), E.ncols(), exp_maxlen);
    tmp_pi.zero_pad(exp_maxlen);

    bool finished = false;

    {
        lingen_qcode_data qq;
        lingen_qcode_init(qq, E.nrows(), E.ncols(), E.get_size(), tmp_pi.get_size());
        for(unsigned int i = 0 ; i < E.nrows() ; i++) {
            for(unsigned int j = 0 ; j < E.ncols() ; j++) {
                lingen_qcode_hook_input(qq, i, j, E.part(i,j));
            }
        }
        for(unsigned int i = 0 ; i < E.ncols() ; i++) {
            for(unsigned int j = 0 ; j < E.ncols() ; j++) {
                lingen_qcode_hook_output(qq, i, j, tmp_pi.part(i,j));
            }
        }
        unsigned int vdelta[E.ncols()];
        int vch[E.ncols()];
        copy(bm.delta.begin(), bm.delta.end(), vdelta);
        copy(bm.lucky.begin(), bm.lucky.end(), vch);
        lingen_qcode_hook_delta(qq, vdelta);
        lingen_qcode_hook_chance_list(qq, vch);
        qq->t = bm.t;
        qq->luck_mini = expected_pi_length(bm.d);
        lingen_qcode_do(qq);
        finished = qq->t < bm.t + qq->length;
        bm.t = qq->t;
        size_t maxlen = 0;
        for(unsigned int j = 0 ; j < tmp_pi.ncols() ; j++) {
            if (lingen_qcode_output_column_length(qq, j) > maxlen)
                maxlen = lingen_qcode_output_column_length(qq, j);
        }
        tmp_pi.truncate(maxlen);
        copy(vdelta, vdelta + E.ncols(), bm.delta.begin());
        copy(vch, vch + E.ncols(), bm.lucky.begin());
        lingen_qcode_clear(qq);
    }
    pi = std::move(tmp_pi);
    return finished;
} /* }}} */

matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E)/*{{{*/
{
    lingen_call_companion const & C = bm.companion(bm.depth(), E.get_size());
    tree_stats::sentinel dummy(bm.stats, "basecase", E.get_size(), C.total_ncalls, true);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    matpoly pi;
    bm.done = bw_lingen_basecase_raw(bm, pi, E);
    bm.stats.end_smallstep();
    E = matpoly();
    ASSERT_ALWAYS(pi.high_word_is_clear());
    return pi;
}/*}}}*/

void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate)/*{{{*/
{
    /* used by testing code */
    bmstatus bm(m,n);
    unsigned int t0 = iceildiv(m,n);
    bm.set_t0(t0);
    matpoly E(ab, m, m+n, L);
    E.zero_pad(L);
    E.fill_random(0, L, rstate);
    matpoly pi;
    bw_lingen_basecase_raw(bm, pi, E);
}/*}}}*/
