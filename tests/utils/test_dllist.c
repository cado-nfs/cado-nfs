#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include "tests_common.h"
#include "macros.h"
#include "dllist.h"
#include "portability.h" // IWYU pragma: keep

struct junk {
    size_t a;
    struct dllist_head head;
};

void
test_dllist(size_t len MAYBE_UNUSED)
{
    struct dllist_head all;

    dllist_init_head(&all);

    for (size_t i = 0; i < len; i++) {
        size_t cur_len = dllist_length(&all);
        if (cur_len != i) {
            fprintf(stderr, "i = %zu, but cur_len = %zu\n", i, cur_len);
            exit(EXIT_FAILURE);
        }
        for (size_t j = 0; j < i; j++) {
            struct junk * node = NULL;
            size_t pos = 0;
            dllist_find_pos(&all, pos, node, struct junk, head, node->a == j);

            /* Must be found */
            if (node == NULL) {
                fprintf(stderr, "i = %zu, j = %zu: dllist_find() did not find node\n", i, j);
                exit(EXIT_FAILURE);
            }

            struct dllist_head * ptr = dllist_get_nth(&all, pos);
            ASSERT_ALWAYS(ptr != NULL);
            struct junk * node2 = dllist_entry(ptr, struct junk, head);

            if (node2 == NULL || node2->a != node->a) {
                fprintf(stderr, "i = %zu, j = %zu: dllist_get_nth() found wrong node\n", i, j);
                exit(EXIT_FAILURE);
            }

            printf("Found i=%zu at position %zu\n", j, pos);
        }

        struct junk * node = NULL;
        dllist_find(&all, node, struct junk, head, node->a == i);
        /* Must not be found */
        if (node != NULL) {
            fprintf(stderr, "i = %zu: dllist_find() incorrectly found a node\n", i);
            exit(EXIT_FAILURE);
        }
        struct junk * foo = malloc(sizeof(struct junk));
        foo->a = i;
        if ((i + rand()) % 2 == 0) {
            dllist_push_front(&all, &foo->head);
        } else {
            dllist_push_back(&all, &foo->head);
        }
    }

    /* Insert a node at head */
    {
        struct junk * foo = malloc(sizeof(struct junk));
        foo->a = len;
        dllist_push_front(&all, &foo->head);
    }
    {
        struct junk * node = NULL;
        dllist_find(&all, node, struct junk, head, node->a == len);
        if (!node || node->a != len) {
            fprintf(stderr, "len = %zu: dllist_find() did not find node after head\n", len);
            exit(EXIT_FAILURE);
        }
    }

    {
        struct dllist_head * ptr = dllist_get_nth(&all, 0);
        ASSERT_ALWAYS(ptr != NULL);
        struct junk * node = dllist_entry(ptr, struct junk, head);
        if (!node || node->a != len) {
            fprintf(stderr, "len = %zu: dllist_get_nth() did not find node after head\n", len);
            exit(EXIT_FAILURE);
        }
    }

    /* Delete nodes again, in random order */
    for (size_t i = 0; i <= len; i++) {
        size_t index = gmp_urandomm_ui(state, len + 1 - i);
        struct dllist_head * ptr = dllist_get_nth(&all, index);
        ASSERT_ALWAYS(ptr != NULL);
        dllist_pop(ptr);
        struct junk * node = dllist_entry(ptr, struct junk, head);
        free(node);
    }

    if (!dllist_is_empty(&all)) {
        fprintf(stderr, "len = %zu: dllist_is_empty() returned false, but list should be empty\n", len);
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
