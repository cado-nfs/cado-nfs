#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#include <stdlib.h>
#include "factor_base.h"
#include "macros.h"

void factor_base_init(factor_base_ptr factor_base, uint64_t number_element_1,
    uint64_t number_element_u, uint64_t number_element_pr)
{
  ASSERT(number_element_1 > 1);
  ASSERT(number_element_u > 1);

  factor_base->number_element_1 = number_element_1;
  factor_base->factor_base_1 = (ideal_1_t *)
    malloc ((number_element_1) * sizeof(ideal_1_t));

  factor_base->number_element_u = number_element_u;
  factor_base->factor_base_u = (ideal_u_t *)
    malloc ((number_element_u) * sizeof(ideal_u_t));

  factor_base->number_element_pr = number_element_pr;
  factor_base->factor_base_pr = (ideal_pr_t *)
    malloc ((number_element_pr) * sizeof(ideal_pr_t));
}

void factor_base_realloc(factor_base_ptr factor_base,
    uint64_t number_element_1, uint64_t number_element_u,
    uint64_t number_element_pr)
{
  ASSERT(number_element_1 <= factor_base->number_element_1);
  ASSERT(number_element_u <= factor_base->number_element_u);
  ASSERT(number_element_pr <= factor_base->number_element_pr);

  factor_base->factor_base_1 =
    (ideal_1_t * ) realloc(factor_base->factor_base_1,
        number_element_1 * sizeof(ideal_1_t));
  factor_base->number_element_1 = number_element_1;

  factor_base->factor_base_u =
    (ideal_u_t * ) realloc(factor_base->factor_base_u,
        number_element_u * sizeof (ideal_u_t));
  factor_base->number_element_u = number_element_u;

  factor_base->factor_base_pr =
    (ideal_pr_t * ) realloc(factor_base->factor_base_pr,
        number_element_pr * sizeof (ideal_pr_t));
  factor_base->number_element_pr = number_element_pr;
}

void factor_base_clear(factor_base_ptr factor_base, unsigned int t)
{
  for (uint64_t i = 0; i < factor_base->number_element_1; i++) {
    ideal_1_clear(factor_base->factor_base_1[i], t);
  }
  free(factor_base->factor_base_1);

  for (uint64_t i = 0; i < factor_base->number_element_u; i++) {
    ideal_u_clear(factor_base->factor_base_u[i], t);
  }
  free(factor_base->factor_base_u);

  for (uint64_t i = 0; i < factor_base->number_element_pr; i++) {
    ideal_pr_clear(factor_base->factor_base_pr[i], t);
  }
  free(factor_base->factor_base_pr);
}

void factor_base_fprintf(FILE * file, factor_base_srcptr factor_base,
    unsigned int t)
{
  for (uint64_t i = 0; i < factor_base->number_element_1; i++) {
    ideal_1_fprintf(file, factor_base->factor_base_1[i], t);
  }

  for (uint64_t i = 0; i < factor_base->number_element_u; i++) {
    ideal_u_fprintf(file, factor_base->factor_base_u[i], t);
  }

  for (uint64_t i = 0; i < factor_base->number_element_pr; i++) {
    ideal_pr_fprintf(file, factor_base->factor_base_pr[i], t);
  }
}
