#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "modredc_3ul.h"
#include "stage2.h"

typedef struct {
  unsigned long *E;        /* The exponent for stage 1 */
  unsigned long E_mask;    /* Mask (one bit set) of MSB in E[E_nrwords-1] */
  unsigned int E_nrwords;  /* Number of words in exponent */
  unsigned int exp2;       /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  stage2_plan_t stage2;
} pm1_plan_t;

/* Notes:
  If s_1 * s_2 = eulerphi(d), we need s_1 precomputed values, s_2 passes and
  ~(B2-B1)/d steps in each pass. Assuming we can compute V_{k_1}(x+1/x)
  and V_{k_2}(x+1/x) with about 1 multiply per value by using a common
  addition chain for S_1 \cup S_2, we can reduce the precomputation cost to
  only s_1 + s_2 instead of eulerphi(P) if we do only 1 pass. This allows for
  a larger d: with only 1 pass, sqrt(B2-B1) is optimal, with a flexible number
  of passes, (B2-B1)^{2/3} is optimal.
  
  In one pass we process one k_2 \in S_2 and thus primes p = i*d +- j with
  j = k_1 + k_2, k_1 \in S_1. We'd like to be able to bound these j by d/2 
  so that we get compact "block" with no overlap between consecutive blocks.
  This prevents having to start at lower i/continue to larger i than 
  necessary, to be able to write all p, B_1 < p <= B_2, as 
  p = i*d +- (k_1 + k_2). However, there seems to be no way to write the 
  smallest positive representatives <d/2 of units mod d as a sum of two sets:
  e.g. for d=30, we'd like {1,7,11,13}.

  
*/


void pm1_stage1_ul (residueredcul_t, const unsigned long *, const int, 
                 const modulusredcul_t);
int pm1_ul (modintredcul_t, const modulusredcul_t, const pm1_plan_t *);
void pm1_stage1_15ul (residueredc15ul_t, const unsigned long *, const int, 
                    const modulusredc15ul_t);
int pm1_15ul (modintredc15ul_t, const modulusredc15ul_t, const pm1_plan_t *);

void pm1_stage1_2ul2 (residueredc2ul2_t, const unsigned long *, const int, 
                    const modulusredc2ul2_t);
int pm1_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const pm1_plan_t *);
void pm1_stage1_3ul (residueredc3ul_t, const unsigned long *, const int, 
                    const modulusredc3ul_t);
int pm1_3ul (modintredc3ul_t, const modulusredc3ul_t, const pm1_plan_t *);
void pm1_make_plan (pm1_plan_t *, const unsigned int, const unsigned int, int);
void pm1_clear_plan (pm1_plan_t *);
