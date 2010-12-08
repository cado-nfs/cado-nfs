#ifndef STAGE2_H
#define STAGE2_H

#define NEXT_D 254
#define NEXT_PASS 255

typedef struct {
  unsigned int B2;
  unsigned int d;          /* The d parameter for stage 2 */
  unsigned int s1;
  unsigned int i0, i1;     /* Stage 2 needs V_{id}(x+1/x) for i0 <= i < i1 */
  int *S1;                 /* Contains eulerphi(d)/2 elements, for each 
                              k1 in S1 we precompute V_k(x+1/x) */
  unsigned char *pairs;    /* For each element in S2, a list of i,j values 
                              for which to multiply V_{id}(x+1/x) - V_{j}(x+1/x)
                              to accumulator */
} stage2_plan_t;

void stage2_make_plan (stage2_plan_t *, unsigned int, unsigned int, int);
void stage2_clear_plan (stage2_plan_t *);

#endif  /* STAGE2_H */
