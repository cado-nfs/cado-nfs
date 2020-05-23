#include "cado.h" // IWYU pragma: keep
#include "int64_vector.h"
#include "mat_int64.h"
#include "list_double_vector.h"

void list_double_vector_init(list_double_vector_ptr list)
{
  list->length = 0;
  list->v = (double_vector_t * ) malloc(sizeof(double_vector_t) *
      DEFAULT_LENGTH_LIST_DOUBLE_VECTOR);
}

void list_double_vector_add_double_vector(list_double_vector_ptr list, double_vector_srcptr v)
{
  if ((list->length % DEFAULT_LENGTH_LIST_DOUBLE_VECTOR) == 0 && list->length != 0) {
    list->v = realloc(list->v, sizeof(double_vector_t) * (list->length +
          DEFAULT_LENGTH_LIST_DOUBLE_VECTOR));
  }
  double_vector_init(list->v[list->length], v->dim);
  double_vector_set(list->v[list->length], v);
  list->length++;
}

void list_double_vector_clear(list_double_vector_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->dim != 0);
    double_vector_clear(list->v[i]);
  }
  free(list->v);
  memset(list, 0, sizeof(*list));
}

void list_double_vector_fprintf(FILE * file, list_double_vector_srcptr list)
{
  fprintf(file, "[");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      double_vector_fprintf(file, list->v[i]);
    }
  }
  fprintf(file, "]\n");
}

void list_double_vector_extract_mat_int64(list_double_vector_ptr list,
    mat_int64_srcptr matrix)
{
  ASSERT(list->length == 0);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, matrix->NumRows);
  double_vector_t vd_tmp;
  double_vector_init(vd_tmp, v_tmp->dim);

  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    mat_int64_extract_vector(v_tmp, matrix, i);
    int64_vector_to_double_vector(vd_tmp, v_tmp);
    list_double_vector_add_double_vector(list, vd_tmp);
  }

  int64_vector_clear(v_tmp);
  double_vector_clear(vd_tmp);
}

