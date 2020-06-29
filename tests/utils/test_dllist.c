#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include "tests_common.h"
#include "dllist.h"
#include "portability.h" // lrand48 // IWYU pragma: keep

void
test_dllist(size_t len)
{
  dllist head;
  dllist_ptr node;
  
  dll_init(head);
  for (size_t i = 0; i < len; i++) {
    size_t cur_len = dll_length(head);
    if (cur_len != i) {
      fprintf(stderr, "i = %zu, but cur_len = %zu\n", i, cur_len);
      exit(EXIT_FAILURE);
    }
    for (size_t j = 0; j < i; j++) {
      node = dll_find(head, (void *) j);
      /* Must be found */
      if (node == NULL || node->data != (void *) j) {
        fprintf(stderr, "i = %zu, j = %zu: dll_find() did not find node\n", i, j);
        exit(EXIT_FAILURE);
      }

      node = dll_get_nth(head, j);
      if (node == NULL || node->data != (void *)j) {
        fprintf(stderr, "i = %zu, j = %zu: dll_get_nth() found wrong node\n", i, j);
        exit(EXIT_FAILURE);
      }
    }

    node = dll_find(head, (void *) i);
    /* Must not be found */
    if (node != NULL) {
      fprintf(stderr, "i = %zu: dll_find() incorrectly found a node\n", i);
      exit(EXIT_FAILURE);
    }
    dll_append(head, (void *) i);
  }

  /* Insert a node at head */
  dll_insert(head, (void *) len);
  node = dll_find(head, (void *) len);
  if (node == NULL || node->data != (void *) len) {
    fprintf(stderr, "len = %zu: dll_find() did not find node after head\n", len);
    exit(EXIT_FAILURE);
  }
  
  node = dll_get_nth(head, 0);
  if (node == NULL || node->data != (void *) len) {
    fprintf(stderr, "len = %zu: dll_get_nth() did not find node after head\n", len);
    exit(EXIT_FAILURE);
  }

  /* Delete nodes again, in random order */
  for (size_t i = 0; i < len + 1; i++) {
    size_t index = lrand48() % (len + 1 - i);
    node = dll_get_nth(head, index);
    dll_delete(node);
  }
  
  if (!dll_is_empty(head)) {
    fprintf(stderr, "len = %zu: dll_is_empty() returned false, but list should be empty\n", len);
    exit(EXIT_FAILURE);
  }
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 100, i;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  for (i = 0; i < iter; i++)
    test_dllist(i);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
