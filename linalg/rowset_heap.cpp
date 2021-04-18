#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#include <queue>          // for make_heap, pop_heap, push_heap
#include "rowset_heap.h"

int heap_index_compare(const struct bucket * a, const struct bucket * b)
{
    return (int32_t) a->i - (int32_t) b->i;
}

int heap_compare(const struct bucket * a, const struct bucket * b)
{
    /* Returns whether b is a better candidate to be on top of the heap
     * IOW: Whether b is lighter than a
     * */
    if (b->room == 0) return 0;
    if (a->room == 0) return 1;
    return b->s < a->s;
}

struct hcmp {
    bool operator() (const bucket& a, const bucket& b) const {
        return heap_compare(&a, &b);
    }
};

void push_heap(struct bucket * __first, struct bucket * __last)
{
    std::push_heap(__first,  __last, hcmp());
}

void pop_heap(struct bucket * __first, struct bucket * __last)
{
    std::pop_heap(__first,  __last, hcmp());
}

void make_heap(struct bucket * __first, struct bucket * __last)
{
    std::make_heap(__first,  __last, hcmp());
}

/* }}} */
