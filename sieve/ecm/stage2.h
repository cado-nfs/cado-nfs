#ifndef STAGE2_H
#define STAGE2_H

#include <limits.h>

#define PAIR_INCR_V (UINT_MAX-1)
#define PAIR_END UINT_MAX

#ifdef __cplusplus
extern "C" {
#endif

/* See the comments above the stage2_one_w function in stage2.c for details */
typedef struct {
  unsigned int B2;
  unsigned int w;         /* Parameter for the Baby-Step-Giant-Step algorithm */
  unsigned int U_len;     /* len of U array */
  unsigned int *U;        /* list of integer coprime to w below w/2 */
  unsigned int vmin, vmax; /* Upper and lower bound on V */
  unsigned int *pairs;    /* This array contains the indexes of the U values
                           * that should be use in the product. We start with
                           * v=vmin; each time PAIR_INCR_V is seen, we increase
                           * v by 1. The special value PAIR_END is used as a end
                           * marker
                           */
} stage2_plan_t;


/* To compute the cost of stage2, a modular multiplication will be counted as
 * having cost 1, so the cost of the doubling and the differential addition must
 * be given with respect to the cost of one modular multiplication.
 */
struct stage2_cost_s
{
  int is_ecm; /* 1 if we are doing stage2 for ECM, 0 for P-/+1.
               * If 1 it means we need to homogenize the coordinates of the
               * points (all projective points must have the same Z-coordinate)
               * which induces an additional cost.
               */
  double dbl; /* cost of a doubling (for an additive group) */
  double dadd; /* cost of a differential addition (for an additive group) */
};

typedef struct stage2_cost_s stage2_cost_t;


void stage2_make_plan (stage2_plan_t *, unsigned int, unsigned int,
                       const stage2_cost_t *, int);
void stage2_clear_plan (stage2_plan_t *);


#ifdef __cplusplus
}
#endif

#endif  /* STAGE2_H */
