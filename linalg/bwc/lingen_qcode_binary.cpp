#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_qcode_binary.hpp"
#include "bpack.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_bw_dimensions.hpp"
#include "lingen_call_companion.hpp"
#include "lingen_matpoly_select.hpp"
#include "macros.h"
#include "timing.h"
#include "tree_stats.hpp"
#include "lingen_expected_pi_length.hpp"
#include "bblas.hpp"

static_assert(std::is_same_v<matpoly<true>::elt, unsigned long>, "wrong flags");
static_assert(std::is_same_v<matpoly<true>::ptr, unsigned long *>, "wrong flags");

/* We have two interfaces here. The first one is the one that is common
 * with qcode_prime. This goes through bw_lingen_basecase.
 *
 * The second one
 * is really legacy code, and it's only exposed for
 * the legacy lingen_binary code. (as a matter of fact, the
 * implementation of bw_lingen_basecase does build upon the legacy
 * interface presently, but that is not a reason to have it exposed).
 */
#define LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE
#ifdef LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE
/* This trampoline structure is no longer useful, really. At some point
 * we had both C and C++ 
 */
struct lingen_qcode_data_s {
    /* we have a matrix of size m*b, and one of size b*b. The second
     * dimension is not called n in order to avoid confusion with the n
     * in BW ; in fact we have b = m + n */
    unsigned int m, b;
    unsigned int t;
    unsigned long length, outlength;
    unsigned int luck_mini;

    unsigned int * pi_col_deg;

    /* where we grab our input and store our output */
    /* Note that we don't own the data corresponding to the innermost
     * level (while we do own the outermost level for iptrs and optrs).
     */
    unsigned int * delta;
    int * ch;

    const unsigned long ** iptrs;
    unsigned long ** optrs;
};


typedef struct lingen_qcode_data_s lingen_qcode_data[1];
typedef struct lingen_qcode_data_s * lingen_qcode_data_ptr;
typedef const struct lingen_qcode_data_s * lingen_qcode_data_srcptr;

static void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength);
static void lingen_qcode_clear(lingen_qcode_data_ptr qq);

static unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j);

static unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq);

static inline void lingen_qcode_hook_delta(lingen_qcode_data_ptr qq, unsigned int * delta)
{
    qq->delta = delta;
}

static inline void lingen_qcode_hook_chance_list(lingen_qcode_data_ptr qq, int * ch)
{
    qq->ch = ch;
}

static inline void lingen_qcode_hook_input(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->iptrs[i * qq->b + j] = poly;
}

static inline void lingen_qcode_hook_output(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->optrs[i * qq->b + j] = poly;
}
#endif

void lingen_qcode_init [[maybe_unused]] (lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength)
{
    qq->m = m;
    qq->b = b;
    qq->length = length;
    qq->outlength = outlength;
    qq->iptrs = (const unsigned long **) malloc(m * b * sizeof(unsigned long *));
    qq->optrs = (unsigned long **) malloc(b * b * sizeof(unsigned long *));
    memset(qq->iptrs, 0, m * b * sizeof(unsigned long *));
    memset(qq->optrs, 0, b * b * sizeof(unsigned long *));

    qq->pi_col_deg = (unsigned int *) malloc(b * sizeof(unsigned int));
    memset(qq->pi_col_deg, 0, b * sizeof(unsigned int));
}
void lingen_qcode_clear [[maybe_unused]] (lingen_qcode_data_ptr qq)
{
    free(qq->iptrs);
    free(qq->optrs);
    free(qq->pi_col_deg);
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
    return qq->pi_col_deg[j] + 1;
}

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
    static_assert(std::is_same_v<pointer_type, unsigned long *> ||
            std::is_same_v<pointer_type, const unsigned long *>,
            "works only for ulong* or const ulong*");
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    bitarray operator^=(bitarray const & a)
    requires std::is_same_v<pointer_type, unsigned long *>
    {
        mpn_xor_n (x, x, a.x, width);
        return *this;
    }
    void lshift1()
    requires std::is_same_v<pointer_type, unsigned long *>
    {
        mpn_lshift(x, x, width, 1);
    }
    void copy_from(const unsigned long * a)
    requires std::is_same_v<pointer_type, unsigned long *>
    {
        memcpy(x, a, width * sizeof(unsigned long));
    }
    void copy_to(unsigned long * a) const {
        memcpy(a, x, width * sizeof(unsigned long));
    }
    bool operator[](size_t k) const {
        return x[k / ULONG_BITS] & (1UL << (k % ULONG_BITS));
    }
    void set1()
    requires std::is_same_v<pointer_type, unsigned long *>
    { *x = 1; }
};

template<>
class bitarray<constant_width<1>, unsigned long *> : private constant_width<1> {
    using width_type = constant_width<1>;
    using pointer_type = unsigned long *;
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    bitarray operator^=(bitarray const & a) { *x ^= *a.x; return *this; }
    void lshift1() { *x <<= 1; }
    void copy_from(const unsigned long * a) { *x = *a; }
    void copy_to(unsigned long * a) const { *a = *x; }
    bool operator[](size_t k) const {
        return *x & (1UL << k);
    }
    void set1() { *x = 1; }
};

template<>
class bitarray<constant_width<1>, const unsigned long *> : private constant_width<1> {
    using width_type = constant_width<1>;
    using pointer_type = unsigned long *;
    using width_type::width;
    pointer_type x;
public:
    bitarray(width_type w, pointer_type x) : width_type(w), x(x) {}
    void copy_to(unsigned long * a) const { *a = *x; }
    bool operator[](size_t k) const {
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
    bitarray<width_type, unsigned long *> operator()(unsigned int i, unsigned int j) {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (i * n + j) * width);
    }
    bitarray<width_type, const unsigned long *> operator()(unsigned int i, unsigned int j) const {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (i * n + j) * width);
    }
    // unsigned long (*row(unsigned int i))[WIDTH] { return x + i * n; }
    // const unsigned long (*row(unsigned int i)const)[WIDTH] { return x + i * n; }
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
    bitarray<width_type, unsigned long *> operator()(unsigned int i, unsigned int j) {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (j * m + i) * width);
    }
    bitarray<width_type, const unsigned long *> operator()(unsigned int i, unsigned int j) const {
        return bitarray<width_type, unsigned long *>((width_type) *this, x + (j * m + i) * width);
    }
    bitarray<variable_width, unsigned long *> column(unsigned int j) {
        return bitarray<variable_width, unsigned long *>(m * width, x + (j * m) * width);
    }
    const unsigned long * column(unsigned int j) const { return (*this)(0, j); }
    void lshift1_column(unsigned int j) {
        unsigned long * c = x + (j * m) * width;
        mpn_lshift(c, c, width, 1);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            c[i*width] &= ~1UL;
    }
    void xor_column(unsigned int k, unsigned int pivot) {
#if 1
        column(k) ^= column(pivot);
#else
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for (unsigned int l = 0; l < m; l++)
            (*this)(l, k) ^= (*this)(l, pivot);
#endif
    }
};

/* This re-implementation is inspired from the one in
 * lingen_qcode_prime.cpp */
static bool generator_found(lingen_qcode_data_ptr qq, unsigned int dt, std::vector<bool>const& is_modified)
{
    unsigned int const m = qq->m;
    unsigned int const b = qq->b;
    unsigned int const n = b - m;
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

        unsigned int const luck_mini = qq->luck_mini;
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

#if 0
static inline void lshift1(unsigned long&x) { x <<= 1; }

template<int WIDTH>
static inline void lshift1(unsigned long (&x)[WIDTH])
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
#endif

template<typename width_type>
static unsigned int lingen_qcode_do_tmpl [[maybe_unused]] (width_type w, lingen_qcode_data_ptr qq)
{
    unsigned int const m = qq->m;
    unsigned int const b = qq->b;
    unsigned int const n = b - m;
    std::vector<unsigned int> kk(m+n, 0);
    int const width = w.width;
    ulmat_colmajor<width_type> E(w, m, b);
    ulmat_colmajor<width_type> P(w, b, b);
    ASSERT_ALWAYS(qq->length <= (unsigned long) width * ULONG_BITS);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (unsigned int i = 0; i < m; i++)
        for(unsigned int j = 0 ; j < b ; j++)
            E(i, j).copy_from(qq->iptrs[i * b + j]);

    for (unsigned int i = 0; i < b; i++) {
        P(i, i).set1();
        qq->pi_col_deg[i] = 0;
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
        std::vector<bool> is_modified(m + n, false);
	for (unsigned int i = 0; i < m; i++) {
	    unsigned int min_degree = UINT_MAX;
            unsigned int pivot = m + n;
	    for (unsigned int j = 0; j < m + n; j++)
		if (E(i, j)[e] && (qq->delta[j] < min_degree)) {
		    min_degree = qq->delta[j];
		    pivot = j;
		}
	    ASSERT_ALWAYS(pivot < m + n);
            is_modified[pivot] = true;

            kk.clear();
	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k)[e])) continue;
                is_modified[k] = true;
                qq->pi_col_deg[k] = std::max(qq->pi_col_deg[k], qq->pi_col_deg[pivot]);
                kk.push_back(k);
            }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (auto const k : kk) {
                E.xor_column(k, pivot);
                P.xor_column(k, pivot);
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
            qq->pi_col_deg[pivot]++;
	}
        /* Columns that have not been changed here are those which
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

unsigned int lingen_qcode_do [[maybe_unused]] (lingen_qcode_data_ptr qq)/*{{{*/
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
}/*}}}*/

static matpoly<true> bw_lingen_basecase_raw_old [[maybe_unused]] (bmstatus<true> & bm, matpoly<true> & E)/*{{{*/
{
    /* There's a nasty bug. Revealed by 32-bits, but can occur on larger
     * sizes too. Let W be the word size. When E has length W + epsilon,
     * pi_deg_bound(W+epsilon) may be < W. So that we may have pi and
     * E have stride 1 and 2 respectively. However, in
     * lingen_qcode_do_tmpl, we fill the pointers in pi (optrs[]
     * there) using the stride of E. This is not proper. A kludge is to
     * make them compatible.
     *
     * (A priori we're not talking immense sizes at the recursive threshold anyway).
     */

    /* (expected_pi_length should do as well, but 1+E.get_size() is
     * firmly on the safe side !. Note that for E.get_size() == 0, we do
     * need length 2.
     */
    size_t const exp_maxlen = 1 + E.get_size();

    matpoly<true> pi(&bm.d.ab, E.ncols(), E.ncols(), exp_maxlen);
    pi.zero_pad(exp_maxlen);

    bool finished = false;

    {
        lingen_qcode_data qq;
        lingen_qcode_init(qq, E.nrows(), E.ncols(), E.get_size(), pi.get_size());
        for(unsigned int i = 0 ; i < E.nrows() ; i++) {
            for(unsigned int j = 0 ; j < E.ncols() ; j++) {
                lingen_qcode_hook_input(qq, i, j, E.part(i,j));
            }
        }
        for(unsigned int i = 0 ; i < E.ncols() ; i++) {
            for(unsigned int j = 0 ; j < E.ncols() ; j++) {
                lingen_qcode_hook_output(qq, i, j, pi.part(i,j));
            }
        }
        unsigned int vdelta[E.ncols()];
        int vch[E.ncols()];
        std::ranges::copy(bm.delta, vdelta);
        std::ranges::copy(bm.lucky, vch);
        lingen_qcode_hook_delta(qq, vdelta);
        lingen_qcode_hook_chance_list(qq, vch);
        qq->t = bm.t;
        qq->luck_mini = expected_pi_length(bm.d);
        lingen_qcode_do(qq);
        finished = qq->t < bm.t + qq->length;
        bm.t = qq->t;
        unsigned int maxlen = 0;
        for(unsigned int j = 0 ; j < pi.ncols() ; j++) {
            maxlen = std::min(maxlen, lingen_qcode_output_column_length(qq, j));
        }
        pi.truncate(maxlen);
        copy(vdelta, vdelta + E.ncols(), bm.delta.begin());
        copy(vch, vch + E.ncols(), bm.lucky.begin());
        lingen_qcode_clear(qq);
    }
    bm.done = finished;
    return pi;
} /* }}} */

static bool generator_found(unsigned int t, bpack_view<uint64_t> E_t, std::vector<int> & lucky, unsigned int luck_mini)
{
    unsigned int const b = E_t.nrows();
    unsigned int const m = E_t.ncols();
    unsigned int const n = b - m;
    unsigned int const bb = E_t.nrowblocks();
    unsigned int const mb = E_t.ncolblocks();
    constexpr const unsigned int B = mat64::width;
    unsigned int newluck = 0;
    /* find spontaneous zero rows */
    for(unsigned int bi = 0 ; bi < bb ; bi++) {
        std::array<bool, B> zz;
        std::fill_n(zz.begin(), B, true);
        unsigned int nnz = 0;
        for(unsigned int bj = 0 ; bj < mb && nnz < B ; bj++) {
            for(unsigned int i = 0 ; i < B ; i++) {
                if (!zz[i]) continue;
                if (E_t.cell(bi, bj)[i]) {
                    nnz++;
                    zz[i] = false;
                }
            }
        }
        for(unsigned int i = 0 ; i < B ; i++) {
            unsigned int const ii = bi * B + i;
            if (!zz[i]) {
                lucky[ii] = 0;
            } else {
                lucky[ii]++;
                newluck++;
            }
        }
    }

    if (!newluck) return false;

    /* If newluck == n, then we probably have a generator. We add an
     * extra guarantee. newluck==n, for a total of k iterations in a
     * row, means m*n*k coefficients cancelling magically. We would
     * like this to be impossible by mere chance. Thus we want n*k >
     * luck_mini, which can easily be checked */

    unsigned int luck_sure = 0;
    int happy = 0;

    printf("t=%d, canceled columns:", t);
    unsigned int last = UINT_MAX;
    unsigned int afterlast = UINT_MAX;
    for(unsigned int j = 0 ; j < b ; j++) {
        if (lucky[j] > 0) {
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
            luck_sure += ((unsigned int) lucky[j]) >= luck_mini;
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

    return happy;
}

static matpoly<true> bw_lingen_basecase_raw_fast(bmstatus<true> & bm, matpoly<true> const & mp_E)/*{{{*/
{
    bw_dimensions<true> & d = bm.d;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const b = m + n;
    unsigned int const L = mp_E.get_size();
    /* expected_pi_length should do as well, but 1+E.get_size() is
     * firmly on the safe side !. Note that for E.get_size() == 0, we do
     * need length 2.
     */
    size_t const D = 1 + mp_E.get_size();

    matpoly<true>::arith_hard * ab = &d.ab;
    constexpr const unsigned int B = mat64::width;
    unsigned int bb = iceildiv(b, B);
    unsigned int const bX = bb * B;
    unsigned int mb = iceildiv(m, B);
    unsigned int const mX = mb * B;
    unsigned int const Lb = iceildiv(L, 64);
    unsigned int const LX = Lb * 64;
    unsigned int const Db = iceildiv(D, B);
    unsigned int const DX = Db * B;

    ASSERT_ALWAYS(Lb * sizeof(uint64_t) == mp_E.data_entry_alloc_size_in_bytes());

    mat64::vector_type E(bb*mb*LX);
    auto E_coeff=[&](unsigned int k) {
        return bpack_view<uint64_t>(&E[k*bb*mb], bb, mb);
    };

    ASSERT_ALWAYS(mX == m);
    ASSERT_ALWAYS(bX == b);

    binary_matpoly_transpose_to_polmat(E.data(),
            (unsigned long const *) mp_E.data_area(),
            mX, bX, LX);

    mat64::vector_type pi(bb * bb * DX);
    std::fill_n(std::begin(pi), pi.size(), 0);
    auto pi_coeff=[&](unsigned int k) {
        return bpack_view<uint64_t>(&pi[k*bb*bb], bb, bb);
    };
    pi_coeff(0).set(1);
    // std::vector<unsigned int> pi_row_deg(b, 0);
    auto & delta(bm.delta);
    auto & lucky(bm.lucky);
    unsigned int const luck_mini = expected_pi_length(bm.d);

    unsigned int pi_len = 1;

    unsigned int t = 0;

    for( ; t < L ; t++, bm.t++) {
        /* invariant: pi * E_orig = X^t * E */

        bpack_view<uint64_t> E_t = E_coeff(t);
        
        if (generator_found(bm.t, E_t, lucky, luck_mini))
            break;

        std::vector<unsigned int> p = E_t.ple(delta);
        for(unsigned int k = 0 ; k < pi_len ; k++)
            pi_coeff(k).propagate_row_permutations(p);
        for(unsigned int k = t + 1 ; k < L ; k++)
            E_coeff(k).propagate_row_permutations(p);
        for(unsigned int ii = 0 ; ii < p.size() ; ii++) {
            std::swap(delta[ii], delta[p[ii]]);
            std::swap(lucky[ii], lucky[p[ii]]);
        }

        bpack<uint64_t> LL(bX, mX);
        bpack<uint64_t>::extract_LU(LL.view(), E_t.view());
        LL.invert_lower_triangular();

        /* This is really the expensive part */
        for(unsigned int k = 0 ; k < pi_len ; k++)
            bpack<uint64_t>::mul_lt_ge(LL.const_view(), pi_coeff(k));
        for(unsigned int k = t + 1 ; k < L ; k++)
            bpack<uint64_t>::mul_lt_ge(LL.const_view(), E_coeff(k));
        
        /* multiply the first p.size() rows by X */
        unsigned int const bi0 = p.size() / 64;
        int const full = pi_len == DX;
        if (bi0) {
            /* we can move complete blocks */
            for(unsigned int d = pi_len - full ; d-- ; ) {
                std::copy(
                        &pi_coeff(d).cell(0,0),
                        &pi_coeff(d).cell(bi0,0),
                        &pi_coeff(d+1).cell(0,0));
            }
            std::fill_n(&pi_coeff(0).cell(0,0), bi0 * bb, 0);
            for(unsigned int d = L - 1 ; d-- > t ; ) {
                std::copy(
                        &E_coeff(d).cell(0,0),
                        &E_coeff(d).cell(bi0,0),
                        &E_coeff(d+1).cell(0,0));
            }
            std::fill_n(&E_coeff(t).cell(0,0), bi0 * mb, 0);
        }
        unsigned int const di = p.size() % 64;
        if (di) {
            /* we can move complete blocks */
            for(unsigned int bj = 0 ; bj < mb ; bj++) {
                for(unsigned int d = pi_len - full ; d-- ; ) {
                    std::copy(
                            pi_coeff(d).cell(bi0,bj).data(),
                            pi_coeff(d).cell(bi0,bj).data() + di,
                            pi_coeff(d+1).cell(bi0,bj).data());
                }
                std::fill_n(pi_coeff(0).cell(bi0,bj).data(), di, 0);
            }
            for(unsigned int bj = 0 ; bj < bb ; bj++) {
                for(unsigned int d = L - 1 ; d-- > t ; ) {
                    std::copy(
                            E_coeff(d).cell(bi0,bj).data(),
                            E_coeff(d).cell(bi0,bj).data() + di,
                            E_coeff(d+1).cell(bi0,bj).data());
                }
                std::fill_n(E_coeff(t).cell(bi0,bj).data(), di, 0);
            }
        }
        if (!full) {
            /* check if pi_len should be increased by 1 or not */
            bool incr = false;
            for(unsigned int bi = 0 ; bi < bb && !incr ; bi++) {
                for(unsigned int bj = 0 ; bj < bb && !incr ; bj++) {
                    incr = pi_coeff(pi_len).cell(bi, bj) != 0;
                }
            }
            pi_len += incr;
        }
        for(unsigned int k = 0 ; k < p.size() ; k++)
            delta[k]++;
    }


    matpoly<true> mp_pi(ab, m+n, m+n, DX);
    binary_polmat_to_matpoly_transpose(
            (unsigned long *) mp_pi.data_area(),
            pi.data(),
            bX, bX, DX);
    mp_pi.set_size(pi_len);

    bm.done = t < L;

    if (0) {
        matpoly<true> const mp_Epi = matpoly<true>::mul(mp_E, mp_pi);
        unsigned int const v = mp_Epi.valuation();
        printf("valuation check: %u\n", v);
        ASSERT_ALWAYS(v >= t);
    }
    return mp_pi;
}/*}}}*/

static matpoly<true> bw_lingen_basecase_raw(bmstatus<true> & bm, matpoly<true> & E)/*{{{*/
{
    return bw_lingen_basecase_raw_fast(bm, E);
}/*}}}*/

matpoly<true> bw_lingen_basecase(bmstatus<true> & bm, matpoly<true> & E)/*{{{*/
{
    lingen_call_companion const & C = bm.companion(bm.depth, E.get_size());
    tree_stats::sentinel const dummy(bm.stats, "basecase", E.get_size(), C.total_ncalls, true);
    bmstatus<true>::depth_sentinel ddummy(bm);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    matpoly<true> pi = bw_lingen_basecase_raw(bm, E);
    bm.stats.end_smallstep();
    E = matpoly<true>();
    ASSERT_ALWAYS(pi.high_word_is_clear());
    return pi;
}/*}}}*/

void test_basecase(matpoly<true>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate)/*{{{*/
{
    /* used by testing code */
    cxx_mpz const p=2;
    bmstatus<true> bm(m,n,p);
    unsigned int const t0 = iceildiv(m,n);
    bm.set_t0(t0);
    matpoly<true> E(ab, m, m+n, L);
    E.zero_pad(L);
    E.fill_random(0, L, rstate);
    bw_lingen_basecase_raw(bm, E);
}/*}}}*/

void test_basecase_bblas(matpoly<true>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate)/*{{{*/
{
    // constexpr const unsigned int B = mat64::width;

    /* used by testing code */
    cxx_mpz const p = 2;
    bmstatus<true> bm(m,n,p);
    unsigned int const t0 = iceildiv(m,n);
    bm.set_t0(t0);

    // ASSERT_ALWAYS(m % B == 0);
    // ASSERT_ALWAYS(n % B == 0);

    matpoly<true> mp_E(ab, m, m+n, L);
    mp_E.zero_pad(L);
    mp_E.fill_random(0, L, rstate);

    double tt;

    tt = wct_seconds();

    matpoly<true> const mp_pi = bw_lingen_basecase_raw_fast(bm, mp_E);

    tt = wct_seconds() - tt;
    printf("%.3f\n", tt);

    {
        auto const mp_Epi = matpoly<true>::mul(mp_E, mp_pi);
        printf("valuation check: %zu\n", mp_Epi.valuation());
    }
}/*}}}*/
