#include "cado.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include "knapsack.h"



struct kns_sum {
    int64_t x;
    unsigned long v;
};

typedef int (*sortfunc_t) (const void *, const void *);

/* returns the sign of the difference a-b (which amounts to comparing a
 * and b */
static inline int int64_cmp(int64_t a, int64_t b)
{
    return (b < a) - (a < b);
}

int kns_sum_cmp(const struct kns_sum *s, const struct kns_sum *t)
{
    int d = int64_cmp(s->x, t->x);
    if (d)
	return d;
    else
	return (s->v > t->v) - (s->v < t->v);
}

int kns_sum_negcmp(const struct kns_sum *s, const struct kns_sum *t)
{
    int d = int64_cmp(-s->x, -t->x);
    if (d)
	return d;
    else
	return (s->v > t->v) - (s->v < t->v);
}

struct kns_sum *all_sums(const int64_t * t, unsigned int t_stride,
			 unsigned long k, size_t extra_alloc,
			 unsigned long offset)
{
    struct kns_sum *r =
	malloc((extra_alloc + (1UL << k)) * sizeof(struct kns_sum));
    int64_t *t2 = malloc(k * sizeof(int64_t));
    int64_t x = offset;
    for (unsigned long i = 0; i < k; i++, t += t_stride) {
	x -= *t;
	t2[i] = *t * 2;
    }
    r[0].x = x;
    unsigned long v = 0;
    r[0].v = 0;
    for (unsigned long i = 1; i < (1UL << k); i++) {
	unsigned int s = __builtin_ffs(i) - 1;
	unsigned long sh = 1UL << s;
	unsigned int down = v & sh;
	v ^= sh;
	x += down ? -t2[s] : t2[s];
	r[i].x = x;
	r[i].v = v;
    }
    free(t2);

    return r;
}

#define ALLOC_EXTRA     256

int knapsack_solve(knapsack_object_ptr ks)
{
    const int64_t * tab = ks->tab;
    unsigned int stride = ks->stride;
    unsigned int offset = ks->offset;
    unsigned int nelems = ks->nelems;
    int64_t bound = ks->bound;
    knapsack_object_callback_t cb = ks->cb;
    void * cb_arg = ks->cb_arg;

    /* purpose: find a combination of the 64-bit integers in tab[], with
     * coefficients -1 or +1, which is within the interval
     * [-bound..bound].
     */

    unsigned int k1 = nelems / 2;
    unsigned int k2 = nelems - k1;
    struct kns_sum *s1 = all_sums(tab + offset, stride, k1, 0, bound);
    struct kns_sum *s2 = all_sums(tab + k1 * stride + offset, stride, k2, 0, 0);

    unsigned int n1 = (1UL << k1);
    unsigned int n2 = (1UL << k2);

    qsort(s1, n1, sizeof(struct kns_sum), (sortfunc_t) & kns_sum_cmp);
    qsort(s2, n2, sizeof(struct kns_sum), (sortfunc_t) & kns_sum_negcmp);

    int64_t ebound = 2 * bound;
    int res = 0;

    /* elements in u1 are sorted in increasing order, in [-B/2,B/2[ */
    /* elements in u2 are sorted in DEcreasing order, or more
     * accurately, in increasing order of their opposite. It is therefore
     * best to consider u2 as containing the OPPOSITES off the y values,
     * even though in reality we don't have to take the negation. These
     * OPPOSITES, for comparison purposes, are thus considered in
     * [-B/2,B/2[
     *
     * We search for x1 - (-x2) in [0,epsilon[ mod B. There are three
     * possible cases:
     *  
     *  x1 - (-x2) in [-B, -B + epsilon[.
     *          this implies x1 in [(-x2)-B, (-x2)-B+epsilon[.
     *          this implies x1 in [-B/2, -B/2+epsilon-1[, provided -x2
     *          has its maximal value B/2-1
     *  x1 - (-x2) in [0, 0 + epsilon[.
     *          this is the main case.
     *  x1 - (-x2) in [B, B + epsilon[.
     *          this implies x1 in [(-x2)+B, (-x2)+B+epsilon[.
     *          Since -x2 >= B/2, this case clearly impossible.
     */

    /* first treat case 1, which is exceptional */
    for (unsigned int u2 = n2; u2--;) {
	/* x2, at the end of the array, is a large negative number. -x2,
	 * (which is not the value present in s2), is at its peak, thus a
	 * large positive number. It must still be very large for the
	 * interval we're going to consider to be non-empty.
	 */
	/*      do we have (-x2)-B+epsilon >= -B/2 ? */
	/* iow, do we have (-x2) >= B/2 - epsilon ? */
	/* iow, do we have (-x2) >= -(-B/2 + epsilon) ? */
	int64_t against = INT64_MIN + ebound;
	if (int64_cmp(-against, -s2[u2].x) > 0)
	    break;
	for (unsigned int u1 = 0; u1 < n1; u1++) {
	    // continue as long as we have:
	    // x1 < (-x2)-B+epsilon
	    // note that we know that (-x2)-B+epsilon >= -B/2, thus
	    // computing (-x2)+epsilon will wrap around to a negative number
	    // in the range [-B/2, -B/2+epsilon[
	    int64_t cut = -s2[u2].x + ebound;	// wraps around.
	    if (int64_cmp(s1[u1].x, cut) >= 0)
		break;

	    unsigned long v = s2[u2].v << k1 | s1[u1].v;
	    int64_t x = s1[u1].x + s2[u2].x;
	    res += cb(cb_arg, v, x - bound);
	}
    }

    unsigned int u1 = 0;

    for (; u1 < n1; u1++) {
	if (s1[u1].x + s2[0].x >= 0) {
	    break;
	}
    }

    for (unsigned int u2 = 0; u2 < n2; u2++) {

	// printf("u2=%d\n",u2);
	/* strategy: first expand the available interval, later restrict
	 * it */

	// compared to the previous turn, x2 has decreased.
	int64_t last_x = s1[u1].x + s2[u2].x;
	for (; u1 < n1; u1++) {
	    int64_t x = s1[u1].x + s2[u2].x;
	    // printf("x=%"PRId64"\n", x);
	    /* if there is _no_ element in s2 such that x1+x2 >= 0, then
	     * we'll notice this since x will drop _below_ last_x. In
	     * such a case, we know that we have a gap of size mor than
	     * B/2 in s2. This gap will be exposed one again for the next
	     * turn.
	     */
	    if (x >= 0 || x < last_x)
		break;
	    // printf("u1=%d\n",u1+1);
	}

	/* Now u1 is the first index such that x1 + x2 >= 0 */
	// assert(u1 == n1 || s2[u2].x + s1[u1].x >= 0);
	// assert(s2[u2].x + s1[u1+pd].x >= ebound);

	for (unsigned int h = 0; u1 + h < n1; h++) {
	    int64_t x = s1[u1 + h].x + s2[u2].x;
	    if (x >= ebound || x < 0)
		break;
	    unsigned long v = s2[u2].v << k1 | s1[u1 + h].v;
	    // fprintf(stderr, "u1=%u u2=%u h=%u: %"PRId64"\n", u1,u2,h,x);
	    assert(x >= 0);
	    assert(x < ebound);
	    res += cb(cb_arg, v, x - bound);
	}
    }
    return res;
}

void knapsack_object_init(knapsack_object_ptr ptr)
{
    memset(ptr, 0, sizeof(knapsack_object));
    ptr->stride = 1;
}

void knapsack_object_clear(knapsack_object_ptr ptr)
{
    memset(ptr, 0, sizeof(knapsack_object));
}

#ifdef  DEMO
int print_solution(knapsack_object_ptr ks, unsigned long v, int64_t x)
{
    char * signs = malloc(ks->nelems + 1);
    memset(signs, '\0', ks->nelems + 1);
    for (unsigned int s = 0; s < ks->nelems; s++)
        signs[s] = (v & (1UL << s)) ? '+' : '-';
    printf("%lx (%s) %ld\n", v, signs, x);
    free(signs);
    return 1;
}

#if 1
#define NELEMS  24
const int64_t tab[NELEMS] = {
    4162221379340189282L, -6647306186101789600L, -1421709520630629187L,
    7978579249304052465L, -5210946216003197439L, 1071743655218434855L,
    2848467872950476511L, -1370801619961069543L, -480116246366186631L,
    7246359761961066352L, -3820114215392891062L, 1960455329265570848L,
    2169371464082239491L, 6918027011352649575L, -687610025789514251L,
    2178270899400006382L, -2751086820472252564L, 326442006929102621L,
    -7009969887660261263L, -5003156455825490387L, 76450565619814227L,
    7595450547102556048L, 4069562109599364928L, -510019920501521658L,
};

const int winning[NELEMS] =
    { 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0,
1 };
    // 0x9b12e5
#elif 1
// passes in 25 seconds on a U9400 1.4GHz.
#define NELEMS 48
const int64_t tab[NELEMS] = {
    -1538040154316158950L, 8258967369930019632L, -4004688999714422921L,
    1159889400949648397L, 5249521295396238696L, 3650114626103672295L,
    -7931626704593073451L, 7674135400418350404L, -2222438450190335011L,
    7860030247650033749L, 4262500207389102065L, 7160554140656125087L,
    3467957067211623500L, -2390204226671589277L, 7688679274618040461L,
    5203982469041323703L, -6554011846537959286L, -1746110400544620451L,
    -5555145176637642721L, -1550913742863465293L, 8257438977466066771L,
    -438015496418932509L, -2243959839000145910L, -3640817532474311269L,
    -1990598364962709281L, 8659906920932625770L, 664803566864969265L,
    2489632956508442703L, 2742566267756791913L, -4145814854589824541L,
    -4289798920807923686L, 4733442344418006680L, 429022769516086920L,
    2715560624548339171L, 6738808048152006654L, 4539026061699952726L,
    -4489030772365149543L, 8876121350766866461L, 1804278128331619321L,
    -2057342976905341813L, -6373799295888301404L, -8752958562853688943L,
    -7543544823992563848L, -6571543819875671449L, 1291491683018881388L,
    3459301977037432953L, 7219978704923241933L, 7852402516694331902L
};
#else
#define NELEMS  5
const int64_t tab[NELEMS] = {
    3697953231077617412, -9135238050079885539, -187993369997840440,
    -7852759334484681542, 2603467885480253854
};
#endif

int main()
{
    knapsack_object ks;
    knapsack_object_init(ks);
    ks->tab = tab;
    ks->nelems = NELEMS;
    ks->bound = 6;
    ks->cb = (knapsack_object_callback_t) print_solution;
    ks->cb_arg = ks;
    knapsack_solve(ks);
    knapsack_object_clear(ks);
    return 0;
}
#endif
