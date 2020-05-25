#ifndef CADO_UTILS_DLLLIST_H
#define CADO_UTILS_DLLLIST_H

/* Need size_t */
#include <stddef.h>
#include <stdlib.h>
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

struct dllist_s {
 void *data;
 struct dllist_s *prev, *next;
};
typedef struct dllist_s dllist[1];
typedef struct dllist_s *dllist_ptr;

static inline void
dll_init(dllist node) {
  node->data = node->prev = node->next = NULL;
}

/* Invariants: 
   A node has prev==NULL iff it is the head node.
   A node has next==NULL iff it is the last node. 
   The head node has data==NULL. */

static inline int
dll_is_empty(dllist head) {
  return (head->next == NULL);
}

static inline void
dll_insert(dllist node, void *data) {
  dllist_ptr new_node = (dllist_ptr) malloc(sizeof(dllist));
  ASSERT_ALWAYS(new_node != NULL);
  new_node->prev = node;
  new_node->next = node->next;
  new_node->data = data;
  node->next = new_node;
  if (new_node->next != NULL)
    new_node->next->prev = new_node;
}

static inline void
dll_append (dllist head, void *data) {
 while (head->next != NULL)
   head = head->next;
 dll_insert(head, data);
}

static inline void
dll_delete (dllist node) {
  /* Can't delete head node */
  ASSERT_ALWAYS(node->prev != NULL);
  node->prev->next = node->next;
  if (node->next != NULL)
    node->next->prev = node->prev;
  free(node);
}

static inline size_t
dll_length (dllist head) {
  size_t len = 0;
 while (head->next != NULL) {
   head = head->next;
   len ++;
  }
  return len;
}

/* Get the n-th node. Requires n < dll_length().
   For n=0, returns the first node after head. */
static inline dllist_ptr
dll_get_nth (dllist head, size_t n) {
  dllist_ptr node = head->next;
  for (node = head->next; node != NULL && n != 0; n--, node = node->next);
  ASSERT_ALWAYS(node != NULL);
  return node;
}

static inline dllist_ptr
dll_find (dllist head, void *data) {
  dllist_ptr next = head->next;

  while (next != NULL) {
    dllist_ptr node = next;
    if (node->data == data)
      break;
    next = node->next;
  }

  return next;
}

#ifdef __cplusplus
}
#endif

#endif
