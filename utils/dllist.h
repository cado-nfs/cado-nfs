#ifndef CADO_UTILS_DLLLIST_H
#define CADO_UTILS_DLLLIST_H

/*
 * our doubly linked lists are circular. As in the linux kernel, we only
 * want an object to carry a struct dllist_head *somewhere*
 *
 * all locking must be handled externally.
 *
 * references:
 * https://www.oreilly.com/library/view/linux-device-drivers/0596000081/ch10s05.html
 * https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git/tree/include/linux/list.h
 */

#include <stddef.h>     // size_t
#include <stdlib.h>
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

#define xxxDLLIST_DEBUG

struct dllist_head {
#ifdef DLLIST_DEBUG
    int kind;
#endif
 struct dllist_head *prev, *next;
};

#define dllist_entry(ptr__, type__, member__) \
        ((type__ *)((char *)(ptr__)-offsetof(type__,member__)))

#define dllist_for_each(L__, ptr__) \
    for(struct dllist_head * ptr__ = (L__)->next ; ptr__ != (L__) ; ptr__ = ptr__->next)

#define dllist_for_each_safe(L__, ptr__) \
    for(struct dllist_head * ptr__ = (L__)->next, * nptr__ ; nptr__ = ptr__->next, ptr__ != (L__) ; ptr__ = nptr__)

#define list_for_each_entry_safe(pos, n, head, member)			\
	for (pos = list_first_entry(head, typeof(*pos), member),	\
		n = list_next_entry(pos, member);			\
	     !list_entry_is_head(pos, head, member); 			\
	     pos = n, n = list_next_entry(n, member))


#define dllist_find(L__, node__, type__, member__, condition__) do {	\
        node__ = NULL;							\
        dllist_for_each((L__), ptr__) {					\
            node__ = dllist_entry(ptr__, type__, member__);     	\
            if (condition__) break;					\
            node__ = NULL;						\
        }								\
    } while (0)


#define dllist_find_pos(L__, pos, node__, type__, member__, condition__) do {	\
        node__ = NULL;							\
        pos = 0;   							\
        dllist_for_each((L__), ptr__) {					\
            node__ = dllist_entry(ptr__, type__, member__);     	\
            if (condition__) break;					\
            node__ = NULL;						\
            pos++;      						\
        }								\
    } while (0)

static inline int
dllist_is_consistent(struct dllist_head * L MAYBE_UNUSED)
{
#ifdef DLLIST_DEBUG
    if (L->kind != 0) return 0;
    dllist_for_each(L, ptr) {
        if (ptr->kind != 1) return 0;
    }
#endif
    return 1;
}

static inline void
dllist_init_head(struct dllist_head * L) {
#ifdef DLLIST_DEBUG
    L->kind = 0;
#endif
    L->prev = L->next = L;
}

static inline void
dllist_init_node(struct dllist_head * L) {
#ifdef DLLIST_DEBUG
    L->kind = 1;
#endif
    L->prev = L->next = L;
}

static inline int
dllist_is_empty(struct dllist_head * L) {
  return L->prev == L;
}

static inline int
dllist_is_singleton (struct dllist_head * L)
{
    /* no need to check both prev and next */
    return L != L->prev && L->prev == L->next;
}

/* attach new element at head of list.
 */
static inline void
dllist_push_front(struct dllist_head * L, struct dllist_head * new_head)
{
#ifdef DLLIST_DEBUG
  new_head->kind = 1;
#endif
  L->next->prev = new_head;
  new_head->next = L->next;
  new_head->prev = L;
  L->next = new_head;
}

/* attach new element at tail of list.
 */
static inline void
dllist_push_back (struct dllist_head * L, struct dllist_head * new_tail)
{
#ifdef DLLIST_DEBUG
  new_head->kind = 1;
#endif
    L->prev->next = new_tail;
    new_tail->prev = L->prev;
    new_tail->next = L;
    L->prev = new_tail;
}

/* detach element in the list that is pointed to by L.
 *
 * Note that this must not be called from a standalone list head, because
 * it might make the whole list loose!
 */
static inline void
dllist_pop (struct dllist_head * L)
{
    struct dllist_head * o_next = L->next;
    struct dllist_head * o_prev = L->prev;
    o_next->prev = o_prev;
    o_prev->next = o_next;
    L->prev = L->next = L;
}


static inline size_t
dllist_length (struct dllist_head * L) {
    size_t len = 0;
    dllist_for_each(L, ptr) {
        len++;
    }
    return len;
}

/* Get the n-th node, or NULL if the list is not large enough.
 *
 * time O(max(n, len(L)))
 *
 */
static inline struct dllist_head *
dllist_get_nth (struct dllist_head * L, size_t n) {
    size_t len = 0;
    dllist_for_each(L, ptr) {
        if (len++ == n)
            return ptr;
    }
    return NULL;
}

/* Get the first node. The list must not be empty.
 */
static inline struct dllist_head *
dllist_get_first_node (struct dllist_head * L) {
    struct dllist_head * x = L->next;
    return x == L ? NULL : x;
}

/* move the entire contents of L0 to the back of L. L0 is made empty.
 * Note that we expect L0 to be a bare list head. */
    static inline void
dllist_bulk_move_back(struct dllist_head * L, struct dllist_head * L0)
{
    if (!dllist_is_empty(L0)) {
        L->prev->next = L0->next;
        L0->next->prev = L->prev;
        L->prev = L0->prev;
        L0->prev->next = L;
        L0->prev = L0->next = L0;
    }
}

/* move the entire contents of L0 to the head of L. L0 is made empty */
static inline void
dllist_bulk_move_head(struct dllist_head * L, struct dllist_head * L0)
{
    if (!dllist_is_empty(L0)) {
        L->next->prev = L0->prev;
        L0->prev->next = L->next;
        L->next = L0->next;
        L0->next->prev = L;
        L0->prev = L0->next = L0;
    }
}

#ifdef __cplusplus
}
#endif

#endif
