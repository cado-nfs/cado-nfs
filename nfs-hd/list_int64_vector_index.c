#include "cado.h" // IWYU pragma: keep
#include "list_int64_vector_index.h"

void list_int64_vector_index_init(list_int64_vector_index_ptr list,
    unsigned int vector_dim)
{
  memset(list, 0, sizeof(*list));
  list->vector_dim = vector_dim;
}

static void list_int64_vector_index_prepare_write(list_int64_vector_index_ptr list,
    unsigned int index)
{
  if (index >= list->alloc) {
    list->alloc = index + 4 + list->alloc / 4;
    list->v = (int64_vector_index_t * ) realloc(list->v, list->alloc *
        sizeof(int64_vector_index_t));
    for (unsigned int i = list->length; i < list->alloc; i++) {
      int64_vector_init(list->v[i]->vec, list->vector_dim);
      list->v[i]->index = 0;
    }
  }
  if (list->length <= index) {
    list->length = index + 1;
  }
}

void list_int64_vector_index_add_int64_vector_index(
    list_int64_vector_index_ptr list, int64_vector_srcptr v, uint64_t index)
{
  ASSERT(v->dim == list->vector_dim);

  list_int64_vector_index_prepare_write(list, list->length);
  int64_vector_set(list->v[list->length - 1]->vec, v);
  list->v[list->length - 1]->index = index;
}

void list_int64_vector_index_clear(list_int64_vector_index_ptr list)
{
  for (unsigned int i = 0; i < list->alloc; i++) {
    ASSERT(list->v[i]->vec->dim != 0);
    int64_vector_clear(list->v[i]->vec);
  }
  free(list->v);
  memset(list, 0, sizeof(*list));
}

void list_int64_vector_index_fprintf(FILE * file,
    list_int64_vector_index_srcptr list)
{
  fprintf(file, "[");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "%" PRIu64 ": ", list->v[i]->index);
      int64_vector_fprintf(file, list->v[i]->vec);
    }
  }
  fprintf(file, "]\n");
}

void list_int64_vector_index_fprintf_comment(FILE * file,
    list_int64_vector_index_srcptr list)
{
  fprintf(file, "# [");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "# %" PRIu64 ": ", list->v[i]->index);
      int64_vector_fprintf(file, list->v[i]->vec);
    }
    fprintf(file, "# ");
  }
  fprintf(file, "]\n");
}

static int compare_vector_last(const void * p0, const void * p1)
{
  const int64_vector_index_t * v0 = (const int64_vector_index_t * ) p0;
  const int64_vector_index_t * v1 = (const int64_vector_index_t * ) p1;
  
  return (int)((*v0)->vec->c[(*v0)->vec->dim - 1] -
      (*v1)->vec->c[(*v1)->vec->dim - 1]);
}

void list_int64_vector_index_sort_last(list_int64_vector_index_ptr list)
{
  qsort(list->v, list->length, sizeof(list->v[0]), compare_vector_last);
}

void list_int64_vector_index_delete_int64_vector_index(
    list_int64_vector_index_ptr list, unsigned int pos)
{
  for (unsigned int i = pos + 1; i < list->length; i++) {
    int64_vector_set(list->v[i - 1]->vec, list->v[i]->vec);
    list->v[i - 1]->index = list->v[i]->index;
  }
  list->length--;
}

//TODO: too naive
void list_int64_vector_index_remove_duplicate(list_int64_vector_index_ptr list)
{
  for (unsigned int i = 0; i < list->length - 1; i++) {
    for (unsigned int j = i + 1; j < list->length; j++) {
      if (int64_vector_equal(list->v[i]->vec, list->v[j]->vec)) {
        list_int64_vector_index_delete_int64_vector_index(list, j);
      }
    }
  }
}

void list_int64_vector_index_remove_duplicate_sort(
    list_int64_vector_index_ptr list)
{
  list_int64_vector_index_sort_last(list);

  for (unsigned int i = 0; i < list->length - 1; i++) {
    unsigned int current = i + 1;
    while (current < list->length && list->v[i]->vec->c[list->vector_dim - 1]
        == list->v[current]->vec->c[list->vector_dim - 1]) {
      if (int64_vector_equal(list->v[i]->vec, list->v[current]->vec)) {
        list_int64_vector_index_delete_int64_vector_index(list, current);
      }
      current++;
    }
  }
}

void list_int64_vector_index_set(
    list_int64_vector_index_ptr list_new,
    list_int64_vector_index_srcptr list_old)
{
  ASSERT(list_new->vector_dim == list_old->vector_dim);

  list_new->length = list_old->length;
  list_new->alloc = list_old->alloc;
  list_new->v = (int64_vector_index_t * ) realloc(list_new->v,
      list_new->alloc * sizeof(int64_vector_index_t));
  for (unsigned int i = 0; i < list_new->length; i++) {
    int64_vector_init(list_new->v[i]->vec, list_new->vector_dim);
    int64_vector_set(list_new->v[i]->vec, list_old->v[i]->vec);
    list_new->v[i]->index = list_old->v[i]->index;
  }
  for (unsigned int i = list_new->length; i < list_new->alloc; i++) {
    int64_vector_init(list_new->v[i]->vec, list_new->vector_dim);
  }
}
