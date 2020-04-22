#include "cado.h"

#include <algorithm>
#include <vector>
#include <type_traits>

#define LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE
#include "lingen_qcode_binary.hpp"
#include "utils.h"

#include "lingen_expected_pi_length.hpp"
#include "gf2x.h"
#include "bblas.hpp"

/* We have two interfaces here. The first one is the one that is common
 * with qcode_prime. This goes through bw_lingen_basecase.
 *
 * The second one, which is activated if
 * LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE is #defined prior to
 * #including this file, is really legacy code, and it's only exposed for
 * the legacy lingen_binary code. (as a matter of fact, the
 * implementation of bw_lingen_basecase does build upon the legacy
 * interface presently, but that is not a reason to have it exposed).
 */
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

    unsigned int * local_delta;

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

void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength);
void lingen_qcode_clear(lingen_qcode_data_ptr qq);

unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j);

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq);

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
    inline void lshift1_column(unsigned int j) {
        unsigned long * c = x + (j * m) * width;
        mpn_lshift(c, c, width, 1);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            c[i*width] &= ~1UL;
    }
    inline void xor_column(unsigned int k, unsigned int pivot) {
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
    std::vector<unsigned int> kk(m+n, 0);
    int width = w.width;

    ulmat_colmajor<width_type> E(w, m, b);

    ulmat_colmajor<width_type> P(w, b, b);

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

#if 1   /* prepare list of k's in advance ? */
            kk.clear();
	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k)[e])) continue;
                is_modified[k] |= 1;
                if (qq->local_delta[pivot] > qq->local_delta[k])
                    qq->local_delta[k] = qq->local_delta[pivot];
                kk.push_back(k);
            }

#if 0   /* process k's by batches ? */
            unsigned int ikbig = (kk.size() / 2) * 2;
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (unsigned int ik = 0 ; ik < ikbig ; ik += 2) {
                unsigned int k0 = kk[ik];
                unsigned int k1 = kk[ik+1];
                // unsigned int k2 = kk[ik+2];
                // unsigned int k3 = kk[ik+3];
                E.xor_column(k0, pivot);
                E.xor_column(k1, pivot);
                // E.xor_column(k2, pivot);
                // E.xor_column(k3, pivot);
                P.xor_column(k0, pivot);
                P.xor_column(k1, pivot);
                // P.xor_column(k2, pivot);
                // P.xor_column(k3, pivot);
	    }
	    for (unsigned int ik = ikbig ; ik < kk.size() ; ik++) {
                unsigned int k = kk[ik];
                E.xor_column(k, pivot);
                P.xor_column(k, pivot);
            }
#else   /* process k's by batches ? */
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	    for (unsigned int ik = 0 ; ik < kk.size() ; ik++) {
                unsigned int k = kk[ik];
                E.xor_column(k, pivot);
                P.xor_column(k, pivot);
            }
#endif   /* process k's by batches ? */
#else   /* prepare list of k's in advance ? */

	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k)[e])) continue;
                is_modified[k] |= 1;
                if (qq->local_delta[pivot] > qq->local_delta[k])
                    qq->local_delta[k] = qq->local_delta[pivot];
                E.xor_column(k, pivot);
                P.xor_column(k, pivot);
            }
#endif  /* prepare list of k's in advance ? */


#if 0
            E.lshift1_column(pivot);
            P.lshift1_column(pivot);
#else
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
#endif
	    qq->delta[pivot] += 1;
            qq->local_delta[pivot]++;
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

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)/*{{{*/
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

matpoly bw_lingen_basecase_raw(bmstatus & bm, matpoly & E)/*{{{*/
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
    size_t exp_maxlen = 1 + E.get_size();

    matpoly pi(bm.d.ab, E.ncols(), E.ncols(), exp_maxlen);
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
        for(unsigned int j = 0 ; j < pi.ncols() ; j++) {
            if (lingen_qcode_output_column_length(qq, j) > maxlen)
                maxlen = lingen_qcode_output_column_length(qq, j);
        }
        pi.truncate(maxlen);
        copy(vdelta, vdelta + E.ncols(), bm.delta.begin());
        copy(vch, vch + E.ncols(), bm.lucky.begin());
        lingen_qcode_clear(qq);
    }
    bm.done = finished;
    return pi;
} /* }}} */

#if 0
typedef int (*sortfunc_t) (const void*, const void*);

static int lexcmp2(const int x[2], const int y[2])
{
    for(int i = 0 ; i < 2 ; i++) {
        int d = x[i] - y[i];
        if (d) return d;
    }
    return 0;
}
matpoly bw_lingen_basecase_raw_slow(bmstatus & bm, matpoly const & E) /*{{{*/
{
    int generator_found = 0;

    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    abdst_field ab = d.ab;
    ASSERT(E.m == m);
    ASSERT(E.n == b);

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int mi, ma;
    std::tie(mi, ma) = get_minmax_delta(bm.delta);

    unsigned int pi_room_base = expected_pi_length(d, bm.delta, E.get_size());

    matpoly pi(ab, b, b, pi_room_base);
    pi.set_size(pi_room_base);

    /* Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */

    std::vector<unsigned int> pi_lengths(b, 1);
    std::vector<unsigned int> pi_real_lengths(b, 1);
    pi.set_constant_ui(1);

    for(unsigned int i = 0 ; i < b ; i++) {
        pi_lengths[i] = 1;
        pi_lengths[i] += bm.delta[i] - mi;
    }

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    std::vector<unsigned int> pivots(m, 0);
    std::vector<int> is_pivot(b, 0);

    matpoly e(ab, m, b, 1);
    e.set_size(1);

    matpoly T(ab, b, b, 1);
    int (*ctable)[2] = new int[b][2];

    for (unsigned int t = 0; t < E.get_size() ; t++, bm.t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        std::vector<unsigned int> todo;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            /* We should never have to recompute from pi using discarded
             * columns. Discarded columns should always correspond to
             * pivots */
            ASSERT_ALWAYS(bm.lucky[j] >= 0);
            todo.push_back(j);
        }
        /* icc openmp doesn't grok todo.size() as being a constant
         * loop bound */
        unsigned int nj = todo.size();

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            size_t twords = 2 * iceildiv(t + 1, ULONG_BITS);
            unsigned long * tmp = new unsigned long[twords];

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int jl = 0 ; jl < nj ; ++jl) {
                for(unsigned int i = 0 ; i < m ; ++i) {
                    unsigned int j = todo[jl];
                    unsigned int lj = MIN(pi_real_lengths[j], t + 1);
                    /* We want to multiply coefficients at bits
                     * [ t - lj + 1 .. t ] of E
                     * and
                     * [ 0 to lj - 1 ] of pi
                     * to form coefficient t of the input.
                     */
                    unsigned int chop = (t - lj + 1) / ULONG_BITS;
                    unsigned int use_E = iceildiv(t + 1, ULONG_BITS) - chop;
                    unsigned int use_pi = iceildiv(lj, ULONG_BITS);
                    unsigned int q_read = t / ULONG_BITS - chop;
                    unsigned int r_read = t % ULONG_BITS;
                    ASSERT_ALWAYS(use_E + use_pi <= twords);
                    unsigned long bit = 0;
                    for(unsigned int k = 0 ; k < b ; ++k) {
                        gf2x_mul(tmp, E.part(i, k) + chop, use_E, pi.part(k, j), use_pi);
                        bit ^= tmp[q_read];
                    }
                    e.part(i, j)[0] = (bit >> r_read) & 1;
                }
            }
            delete[] tmp;
        }
        /* }}} */

        /* {{{ check for cancellations */

        unsigned int newluck = 0;
        for(unsigned int jl = 0 ; jl < todo.size() ; ++jl) {
            unsigned int j = todo[jl];
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                nz += (e.part(i, j)[0] & 1) == 0;
            }
            if (nz == m) {
                newluck++, bm.lucky[j]++;
            } else if (bm.lucky[j] > 0) {
                bm.lucky[j] = 0;
            }
        }


        if (newluck) {
            /* If newluck == n, then we probably have a generator. We add an
             * extra guarantee. newluck==n, for a total of k iterations in a
             * row, means m*n*k coefficients cancelling magically. We would
             * like this to be impossible by mere chance. Thus we want n*k >
             * luck_mini, which can easily be checked */

            int luck_mini = expected_pi_length(d);
            unsigned int luck_sure = 0;

            printf("t=%d, canceled columns:", bm.t);
            for(unsigned int j = 0 ; j < b ; j++) {
                if (bm.lucky[j] > 0) {
                    printf(" %u", j);
                    luck_sure += bm.lucky[j] >= luck_mini;
                }
            }

            if (newluck == n && luck_sure == n) {
                if (!generator_found) {
                    printf(", complete generator found, for sure");
                }
                generator_found = 1;
            }
            printf(".\n");
        }
        /* }}} */

        if (generator_found) break;

        /* {{{ Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        for(unsigned int j = 0; j < b; j++) {
            ctable[j][0] = bm.delta[j];
            ctable[j][1] = j;
        }
        qsort(ctable, b, 2 * sizeof(int), (sortfunc_t) & lexcmp2);
        /* }}} */

        /* {{{ Now do Gaussian elimination */

        /*
         * The matrix T is *not* used for actually storing the product of
         * the transvections, just the *list* of transvections. Then,
         * instead of applying them row-major, we apply them column-major
         * (abiding by the ordering of pivots), so that we get a better
         * opportunity to do lazy reductions.
         */

        T.set_constant_ui(1);

        is_pivot.assign(b, 0);
        unsigned int r = 0;

        std::vector<unsigned int> pivot_columns;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                unsigned long euj = e.coeff(u, j)[0] & 1;
                if (euj)
                    break;
            }
            if (u == m) continue;
            assert(r < m);
            /* }}} */
            pivots[r++] = j;
            is_pivot[j] = 1;
            pivot_columns.push_back(j);
            /* {{{ Cancel this coeff in all other columns. */
            // abelt inv;
            // abinit(ab, &inv);
            // abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                unsigned long euk = e.coeff(u, k)[0] & 1;
                if (euk == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                // abelt lambda;
                // abinit(ab, &lambda);
                // abmul(ab, lambda, inv, e.coeff(u, k, 0));
                assert(bm.delta[j] <= bm.delta[k]);
                /* {{{ Apply on both e and pi */
                // abelt tmp;
                // abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    e.part(i,k)[0] ^= e.part(i,j)[0];
                }
                if (bm.lucky[k] < 0) {
                    /* This column is already discarded, don't bother */
                    continue;
                }
                if (bm.lucky[j] < 0) {
                    /* This column is discarded. This is going to
                     * invalidate another column of pi. Not a problem,
                     * unless it's been marked as lucky previously ! */
                    ASSERT_ALWAYS(bm.lucky[k] <= 0);
                    printf("Column %u discarded from now on (through addition from column %u)\n", k, j);
                    bm.lucky[k] = -1;
                    continue;
                }
                /* We do *NOT* really update T. T is only used as
                 * storage!
                 */
                T.coeff(j,k)[0] = 1;
                // abset(ab, T.coeff(j, k, 0), lambda);
                // abclear(ab, &tmp);
                /* }}} */
                // abclear(ab, &lambda);
            }
            // abclear(ab, &inv); /* }}} */
        }
        /* }}} */

        /* {{{ apply the transformations, using the transvection
         * reordering trick */

        /* non-pivot columns are only added to and never read, so it does
         * not really matter where we put their computation, provided
         * that the columns that we do read are done at this point.
         */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            if (!is_pivot[j])
                pivot_columns.push_back(j);
        }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            // abelt_ur tmp_pi;
            // abelt_ur tmp;
            // abelt_ur_init(ab, &tmp);
            // abelt_ur_init(ab, &tmp_pi);

            for(unsigned int jl = 0 ; jl < b ; ++jl) {
                unsigned int j = pivot_columns[jl];
                /* compute column j completely. We may put this interface in
                 * matpoly, but it's really special-purposed, to the point
                 * that it really makes little sense IMO
                 *
                 * Beware: operations on the different columns are *not*
                 * independent, here ! Operations on the different degrees,
                 * on the other hand, are. As well of course as the
                 * operations on the different entries in each column.
                 */

#ifndef NDEBUG
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
                {
                    for(unsigned int kl = m ; kl < b ; kl++) {
                        unsigned int k = pivot_columns[kl];
                        unsigned long Tkj = T.part(k, j)[0];
                        ASSERT_ALWAYS(Tkj == (k == j));
                    }
                    for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                        unsigned int k = pivot_columns[kl];
                        unsigned long Tkj = T.part(k, j)[0];
                        if (!Tkj) continue;
                        ASSERT_ALWAYS(pi_lengths[k] <= pi_lengths[j]);
                        pi_real_lengths[j] = std::max(pi_real_lengths[k], pi_real_lengths[j]);
                    }
                }
#endif

                /* Icc 2019 synthetizes a pragma omp single around
                 * accesses to pi_real_lengths. I don't think it makes
                 * sense in that particular case, it's fine enough here
                 * to read the data now after the critical section above.
                 */
                unsigned int dummy = pi_real_lengths[j];

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
                for(unsigned int i = 0 ; i < b ; i++) {
                    for(unsigned int s = 0 ; s < dummy ; s += ULONG_BITS) {

                        unsigned long tmp = 0;

                        for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                            unsigned int k = pivot_columns[kl];
                            /* TODO: if column k was already a pivot on previous
                             * turn (which could happen, depending on m and n),
                             * then the corresponding entry is probably zero
                             * (exact condition needs to be written more
                             * accurately).
                             */

                            unsigned long Tkj = T.part(k, j)[0];
                            if (!Tkj) continue;
                            /* pi[i,k] has length pi_lengths[k]. Multiply
                             * that by T[k,j], which is a constant. Add
                             * to the unreduced thing. We don't have an
                             * mpfq api call for that operation.
                             */
                            tmp ^= pi.part(i, k)[s / ULONG_BITS];
                        }

                        pi.part(i, j)[s / ULONG_BITS] ^= tmp;
                    }
                }
            }
            // abelt_ur_clear(ab, &tmp);
            // abelt_ur_clear(ab, &tmp_pi);
        }
        /* }}} */

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_real_lengths[j] >= pi.capacity()) {
                if (!generator_found) {
                    pi.realloc(pi.capacity() + MAX(pi.capacity() / (m+n), 1));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm.t,
                            100 * pi.capacity() / pi_room_base);
                    for(unsigned int j = 0; j < b; j++)
                        printf(" %u", pi_real_lengths[j]);
                    printf("\n");
                } else {
                    ASSERT_ALWAYS(bm.lucky[j] <= 0);
                    if (bm.lucky[j] == 0)
                        printf("t=%u, column %u discarded from now on\n",
                                bm.t, j);
                    bm.lucky[j] = -1;
                    pi_lengths[j]++;
                    pi_real_lengths[j]++;
                    bm.delta[j]++;
                    continue;
                }
            }
            pi.multiply_column_by_x(j, pi_real_lengths[j]);
            pi_real_lengths[j]++;
            pi_lengths[j]++;
            bm.delta[j]++;
        }
        /* }}} */
    }

    unsigned int pisize = 0;
    for(unsigned int j = 0; j < b; j++) {
        if (pi_real_lengths[j] > pisize)
            pisize = pi_real_lengths[j];
    }
    /* Given the structure of the computation, there's no reason for the
     * initial estimate to go wrong.
     */
    ASSERT_ALWAYS(pisize <= pi.capacity());
    pi.set_size(pisize);

    for(unsigned int j = 0; j < b; j++) {
        for(unsigned int k = pi_real_lengths[j] ; k < pi.get_size() ; k++) {
            for(unsigned int i = 0 ; i < b ; i++) {
                ASSERT_ALWAYS(abis_zero(ab, pi.coeff(i, j, k)));
            }
        }
    }

    bm.done = generator_found;
    return pi;
}
/* }}} */
#endif

matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E)/*{{{*/
{
    lingen_call_companion const & C = bm.companion(bm.depth(), E.get_size());
    tree_stats::sentinel dummy(bm.stats, "basecase", E.get_size(), C.total_ncalls, true);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    matpoly pi = bw_lingen_basecase_raw(bm, E);
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
    bw_lingen_basecase_raw(bm, E);
}/*}}}*/

void test_basecase_bblas(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate)/*{{{*/
{
    constexpr const unsigned int B = mat64::width;
    // typedef mat64::datatype U;

    /* used by testing code */
    bmstatus bm(m,n);
    unsigned int t0 = iceildiv(m,n);
    bm.set_t0(t0);

    // ASSERT_ALWAYS(m % B == 0);
    // ASSERT_ALWAYS(n % B == 0);

    unsigned int b = m + n;
    unsigned int bb = iceildiv(b, B);
    unsigned int bX = bb * B;
    unsigned int mb = iceildiv(m, B);
    unsigned int mX = mb * B;
    unsigned int Lb = iceildiv(L, B);
    unsigned int LX = Lb * B;
    unsigned int Db = iceildiv(iceildiv(m * L, b), B);
    unsigned int DX = Db * B;


    mat64::vector_type E(bb*mb*LX);
    auto E_coeff=[&](unsigned int k) {
        return bpack_view<mat64>(&E[k*bb*mb], bb, mb);
    };

    matpoly mp_E(ab, m+n, m, L);
    mp_E.zero_pad(L);
    mp_E.fill_random(0, L, rstate);

    double tt;

    tt = wct_seconds();
    ASSERT_ALWAYS(mp_E.data_entry_size_in_bytes() == mp_E.data_entry_alloc_size_in_bytes());
    ASSERT_ALWAYS(mp_E.data_entry_size_in_words() == iceildiv(L,B));
    binary_matpoly_to_polmat(&E[0],
            (uint64_t const *) mp_E.data_area(),
            bX, mX, LX);

    mat64::vector_type pi(bb * bb * DX);
    std::fill_n(std::begin(pi), pi.size(), 0);
    auto pi_coeff=[&](unsigned int k) {
        return bpack_view<mat64>(&pi[k*bb*bb], bb, bb);
    };
    pi_coeff(0).set(1);

    std::vector<unsigned int> d(bX, 0);

    for(unsigned int t = 0 ; t < L ; t++) {
        /* invariant: pi * E_orig = X^t * E */
        /*
        {
            tt += wct_seconds();
            matpoly mp_pi(ab, m+n, m+n, DX);
            binary_polmat_to_matpoly(
                    (uint64_t *) mp_pi.data_area(),
                    &pi[0],
                    bX, bX, DX);

            unsigned int pi_len = 1 + *std::max_element(std::begin(d), std::end(d));
            mp_pi.set_size(pi_len);
            matpoly mp_Epi = matpoly::mul(mp_pi, mp_E);

            printf("t=%u val=%u\n", t, mp_Epi.valuation());
            tt -= wct_seconds();
        }
        */
        bpack_view<mat64> E_t = E_coeff(t);
        /*
        bpack<mat64> E_t_copy(E_t.nrows(), E_t.ncols());
        E_t_copy.view().set(E_t);
        */
        std::vector<unsigned int> p = E_t.ple();
        unsigned int pi_len = 1 + *std::max_element(std::begin(d), std::end(d));

        for(unsigned int k = 0 ; k < pi_len ; k++)
            pi_coeff(k).propagate_row_permutations(p);
        for(unsigned int k = t + 1 ; k < L ; k++)
            E_coeff(k).propagate_row_permutations(p);
        for(unsigned int ii = 0 ; ii < p.size() ; ii++)
            std::swap(d[ii], d[p[ii]]);

        /*
        E_t_copy.propagate_row_permutations(p);
        */

        bpack<mat64> LL(bX, mX);
        bpack<mat64>::extract_LU(LL.view(), E_t.view());
        LL.invert_lower_triangular();

        /*
        bpack<mat64> UU(bX, mX);
        UU.view().set(E_t);
        bpack<mat64>::mul_lt_ge(LL, E_t_copy);
        ASSERT_ALWAYS(E_t_copy == UU);
        */

        /* This is really the expensive part */
        for(unsigned int k = 0 ; k < pi_len ; k++)
            bpack<mat64>::mul_lt_ge(LL.const_view(), pi_coeff(k));
        for(unsigned int k = t + 1 ; k < L ; k++)
            bpack<mat64>::mul_lt_ge(LL.const_view(), E_coeff(k));
        
        /* multiply the first p.size() rows by X */
        unsigned int bi0 = p.size() / 64;
        int full = pi_len == DX;
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
        unsigned int di = p.size() % 64;
        if (di) {
            /* we can move complete blocks */
            for(unsigned int bj = 0 ; bj < mb ; bj++) {
                for(unsigned int d = pi_len - full ; d-- ; ) {
                    std::copy(
                            &pi_coeff(d).cell(bi0,bj)[0],
                            &pi_coeff(d).cell(bi0,bj)[di],
                            &pi_coeff(d+1).cell(bi0,bj)[0]);
                }
                std::fill_n(&pi_coeff(0).cell(bi0,bj)[0], di, 0);
            }
            for(unsigned int bj = 0 ; bj < bb ; bj++) {
                for(unsigned int d = L - 1 ; d-- > t ; ) {
                    std::copy(
                            &E_coeff(d).cell(bi0,bj)[0],
                            &E_coeff(d).cell(bi0,bj)[di],
                            &E_coeff(d+1).cell(bi0,bj)[0]);
                }
                std::fill_n(&E_coeff(t).cell(bi0,bj)[0], di, 0);
            }
        }
        for(unsigned int k = 0 ; k < p.size() ; k++) {
            d[k]++;
            if (d[k] == DX) d[k]--;
        }
        /*
        {
            tt += wct_seconds();
            matpoly mp_pi(ab, m+n, m+n, DX);
            binary_polmat_to_matpoly(
                    (uint64_t *) mp_pi.data_area(),
                    &pi[0],
                    bX, bX, DX);

            unsigned int pi_len = 1 + *std::max_element(std::begin(d), std::end(d));
            mp_pi.set_size(pi_len);
            matpoly mp_Epi = matpoly::mul(mp_pi, mp_E);
            mp_Epi.realloc(pi_len + LX);
            mp_Epi.zero_pad(pi_len + LX);
            mp_Epi.set_size(LX);
            mp_Epi.shrink_to_fit();

            mat64::vector_type F(bb*mb*LX);
            binary_matpoly_to_polmat(&F[0],
                    (uint64_t const *) mp_Epi.data_area(),
                    bX, mX, LX);
            for(unsigned int s = 0 ; s <= t ; s++) {
                ASSERT_ALWAYS(bpack_view<mat64>(&F[s*bb*mb], bb, mb) == 0);
            }
            printf("t=%u val=%u\n", t, mp_Epi.valuation());
            tt -= wct_seconds();
        }
        */
    }


    matpoly mp_pi(ab, m+n, m+n, DX);
    binary_polmat_to_matpoly(
            (uint64_t *) mp_pi.data_area(),
            &pi[0],
            bX, bX, DX);
    unsigned int pi_len = 1 + *std::max_element(std::begin(d), std::end(d));
    mp_pi.set_size(pi_len);

    tt = wct_seconds() - tt;
    printf("%.3f\n", tt);

    if (0) {
        matpoly mp_Epi = matpoly::mul(mp_pi, mp_E);
        printf("valuation check: %u\n", mp_Epi.valuation());
    }
    // bw_lingen_basecase_raw(bm, E);
}/*}}}*/
