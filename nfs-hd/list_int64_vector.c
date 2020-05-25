#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <limits.h>
#include <inttypes.h>
#include <stdio.h>
#include "list_int64_vector.h"
#include "macros.h"

void list_int64_vector_init(list_int64_vector_ptr list,
    unsigned int vector_dim)
{
  memset(list, 0, sizeof(*list));
  list->vector_dim = vector_dim;
}

/*
 * Increase the length of the list, and if needed, the number of possible
 * allocated cells of the list.
 *
 * list: the list.
 * index: index in which we hope to write an element.
 */
static void list_int64_vector_prepare_write(list_int64_vector_ptr list,
    unsigned int index)
{
  if (index >= list->alloc) {
    list->alloc = index + 4 + list->alloc / 4;
    list->v = (int64_vector_t *) realloc(list->v, list->alloc *
        sizeof(int64_vector_t));
    for (unsigned int i = list->length; i < list->alloc; i++) {
      int64_vector_init(list->v[i], list->vector_dim);
    }
  }
  if (list->length <= index) {
    list->length = index + 1;
  }
}

unsigned int list_int64_vector_add_int64_vector(list_int64_vector_ptr list,
    int64_vector_srcptr v)
{
  ASSERT(v->dim == list->vector_dim);

  list_int64_vector_prepare_write(list, list->length);
  int64_vector_set(list->v[list->length - 1], v);
  return list->length - 1;
}

void list_int64_vector_clear(list_int64_vector_ptr list)
{
  for (unsigned int i = 0; i < list->alloc; i++) {
    ASSERT(list->v[i]->dim != 0);
    int64_vector_clear(list->v[i]);
  }
  free(list->v);
  memset(list, 0, sizeof(*list));
}

void list_int64_vector_delete_elements(list_int64_vector_ptr list)
{
  list->length = 0;
}

void list_int64_vector_fprintf(FILE * file, list_int64_vector_srcptr list)
{
  fprintf(file, "[");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      int64_vector_fprintf(file, list->v[i]);
    }
  }
  fprintf(file, "]\n");
}

void list_int64_vector_fprintf_comment(FILE * file,
    list_int64_vector_srcptr list)
{
  fprintf(file, "# [");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "# ");
      int64_vector_fprintf(file, list->v[i]);
    }
    fprintf(file, "# ");
  }
  fprintf(file, "]\n");
}

void list_int64_vector_extract_mat_int64(list_int64_vector_ptr list,
    mat_int64_srcptr matrix)
{
  ASSERT(list->length == 0);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, matrix->NumRows);

  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    mat_int64_extract_vector(v_tmp, matrix, i);
    list_int64_vector_add_int64_vector(list, v_tmp);
  }

  int64_vector_clear(v_tmp);
}

int int64_vector_in_polytop_list_int64_vector(int64_vector_srcptr vec,
    list_int64_vector_srcptr list)
{
  ASSERT(list->length > 2);
  ASSERT(vec->dim > 1);
#ifndef NDEBUG
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(vec->dim == list->v[i]->dim);
    for (unsigned int j = 2; j < vec->dim; j++) {
      ASSERT(vec->c[j] == list->v[i]->c[j]);
    }
  }
#endif // NDEBUG

  for (unsigned int i = 0; i < list->length; i++) {
    if (!int64_vector_equal(vec, list->v[i])) {
      ASSERT(int64_vector_equal(vec, list->v[i]) == 0);

      return 1;
    }
  }

  int c = 0;
  unsigned i = 0;
  unsigned j = 0;
  for (i = 0, j = list->length - 1; i < list->length; j = i++) {
    if ( ((list->v[i]->c[1] > vec->c[1]) != (list->v[j]->c[1] > vec->c[1])) &&
        (vec->c[0] < (list->v[j]->c[0] - list->v[i]->c[0]) *
         (vec->c[1] - list->v[i]->c[1]) /
         (list->v[j]->c[1] - list->v[i]->c[1]) + list->v[i]->c[0]) ) {
      c = !c;
    }
  }

  return c;
}

static int compare_last(const void * p0, const void * p1)
{
  const int64_vector_t * v0 = (const int64_vector_t * ) p0;
  const int64_vector_t * v1 = (const int64_vector_t * ) p1;

  return (int)((*v0)->c[(*v0)->dim - 1] - (*v1)->c[(*v1)->dim - 1]);
}

void list_int64_vector_sort_last(list_int64_vector_ptr list)
{
  qsort(list->v, list->length, sizeof(list->v[0]), compare_last);
}

/*
 * Assume that list is sorted. Return the index of the first vector with a last
 *  coordinate less or equal to val.
 */
static unsigned int find_index(list_int64_vector_ptr list, int64_t val)
{
  unsigned int i = 0;

  for (i = list->length - 1; i != UINT_MAX; i--) {
    if (val >= list->v[i]->c[list->vector_dim - 1]) {
      return i;
    }
  }
  return 0;
}

void list_int64_vector_delete_last_coordinate(list_int64_vector_ptr list,
    int64_t val)
{
  if (list->length == 0) {
    return;
  }
  unsigned int i = find_index(list, val);

  for (unsigned int j = i + 1; j < list->length; j++) {
    int64_vector_set(list->v[j - (i + 1)], list->v[j]);
  }
  list->length = list->length - (i + 1);
}

//TODO: too naive.
void list_int64_vector_remove_duplicate(list_int64_vector_ptr list)
{
  if (list->length == 0) {
    return;
  }

  for (unsigned int i = 0; i < list->length - 1; i++) {
    for (unsigned int j = i + 1; j < list->length; j++) {
      if (int64_vector_equal(list->v[i], list->v[j])) {
        list_int64_vector_delete_int64_vector(list, j);
      }
    }
  }
}

void list_int64_vector_delete_int64_vector(list_int64_vector_ptr list,
    unsigned int pos)
{
  for (unsigned int i = pos + 1; i < list->length; i++) {
    int64_vector_set(list->v[i - 1], list->v[i]);
  }
  list->length--;
}
