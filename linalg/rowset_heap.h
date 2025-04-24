#ifndef CADO_ROWSET_HEAP_H
#define CADO_ROWSET_HEAP_H

/*  Data type for buckets used in the priority queue stuff (w/ tests) */
/* We'll arrange our slices in a priority heap, so that we'll always
 * easily access the lightest one.
 *
 * The number of rows in a bucket does not really matter, except that we
 * absolutely want to make sure that it never exceeds the bucket size.
 * Hence if a bucket becomes full, then it's always considered heavier
 * than anything.
 */

struct bucket {
    unsigned long s;    /* total weight */
    unsigned long room; /* nb of rows that can still be stored */
    int i;              /* bucket index */
};

#ifdef __cplusplus
extern "C" {
#endif

extern int heap_index_compare(const struct bucket * a, const struct bucket * b);
extern int heap_compare(const struct bucket * a, const struct bucket * b);
extern void push_heap(struct bucket * __first, struct bucket * __last);
extern void pop_heap(struct bucket * __first, struct bucket * __last);
extern void make_heap(struct bucket * __first, struct bucket * __last);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_ROWSET_HEAP_H */
