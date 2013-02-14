#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ASSERT(x)       assert(x)
#define croak__(x,y) do {                                               \
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",                     \
                (x),__func__,__FILE__,__LINE__,(y));                    \
    } while (0)
#define ASSERT_ALWAYS(x)                                                \
    do {                                                                \
        if (!(x)) {                                                     \
            croak__("code BUG() : condition " #x " failed",             \
                    "Abort");                                           \
            abort();                                                    \
        }                                                               \
    } while (0)

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#include "methods.h"


prior_ptr get_prior(int lo, int hi, int fbb);

static const float zero_hist[60] = {0,}; // fill in with zeros.

int main(int argc, char **argv)
{
    int cofac_range[2] = {100, 110};
    int fbb = 22;
    prior_srcptr prior = get_prior(cofac_range[0], cofac_range[1], fbb);

    float acc_failure1[60], acc_failure5[60];
    float acc_failure7[60], acc_failure11[60] ;
    for (int i = 0; i < 60; ++i) {
        acc_failure1[i] = 1;
        acc_failure5[i] = 1;
        acc_failure7[i] = 1;
        acc_failure11[i] = 1;
    }

    ppm1_history_t ppm1_history, old_ppm1_history;
    ppm1_history->pm1_success1 = zero_hist;
    ppm1_history->pm1_success5 = zero_hist;
    ppm1_history->pm1_success7 = zero_hist;
    ppm1_history->pm1_success11 = zero_hist;
    ppm1_history->pp1_success1 = zero_hist;
    ppm1_history->pp1_success5 = zero_hist;
    ppm1_history->pp1_success7 = zero_hist;
    ppm1_history->pp1_success11 = zero_hist;

    for (int i = 0; i < 200; ++i) {
        cofac_method_srcptr method;
        old_ppm1_history->pm1_success1 = ppm1_history->pm1_success1;
        old_ppm1_history->pm1_success5 = ppm1_history->pm1_success5;
        old_ppm1_history->pm1_success7 = ppm1_history->pm1_success7;
        old_ppm1_history->pm1_success11 = ppm1_history->pm1_success11;
        old_ppm1_history->pp1_success1 = ppm1_history->pp1_success1;
        old_ppm1_history->pp1_success5 = ppm1_history->pp1_success5;
        old_ppm1_history->pp1_success7 = ppm1_history->pp1_success7;
        old_ppm1_history->pp1_success11 = ppm1_history->pp1_success11;
        method = get_method_naive(cofac_range, prior,
                acc_failure1, acc_failure5, acc_failure7,
                acc_failure11, ppm1_history);
        int purged_limit = 0;
        for (int j = 0; j < 60; ++j) {
            if (method->type == ECM) {
                acc_failure1[j] *= (1-method->success1[j]);
                acc_failure5[j] *= (1-method->success5[j]);
                acc_failure7[j] *= (1-method->success7[j]);
                acc_failure11[j] *= (1-method->success11[j]);
            } else {
#define DO(a, b, c) do { if (c != 1) { a *= (1-b)/(1-c); } } while (0);
              DO(acc_failure1[j] , ppm1_history->pm1_success1[j], old_ppm1_history->pm1_success1[j]);
              DO(acc_failure5[j] , ppm1_history->pm1_success5[j], old_ppm1_history->pm1_success5[j]);
              DO(acc_failure7[j] , ppm1_history->pm1_success7[j], old_ppm1_history->pm1_success7[j]);
              DO(acc_failure11[j],ppm1_history->pm1_success11[j], old_ppm1_history->pm1_success11[j]);
              DO(acc_failure1[j] , ppm1_history->pp1_success1[j], old_ppm1_history->pp1_success1[j]);
              DO(acc_failure5[j] , ppm1_history->pp1_success5[j], old_ppm1_history->pp1_success5[j]);
              DO(acc_failure7[j] , ppm1_history->pp1_success7[j], old_ppm1_history->pp1_success7[j]);
              DO(acc_failure11[j],ppm1_history->pp1_success11[j], old_ppm1_history->pp1_success11[j]);
#undef DO
            }
            float fail = (acc_failure1[j] + acc_failure5[j] +
                    acc_failure7[j] + acc_failure11[j]) / 4;
            // take the largest bit size for which primes of that size
            // have been found with probability > 90%.
            if (fail < 0.1) 
                purged_limit = j;
        }
        printf("{ ");
        switch (method->type) {
            case ECM:
                printf("EC_METHOD");
                break;
            case PM1:
                printf("PM1_METHOD");
                break;
            case PP1_27:
                printf("PP1_27_METHOD");
                break;
            case PP1_65:
                printf("PP1_65_METHOD");
                break;
        }
        printf(",\t%d,%d, %d },\n", method->B1, method->B2, purged_limit);
    }
    // final Baysian thm
    printf("bits: prob_to_find_p(i) prob_still_exist_p(i)\n");
//    for (int i = 0; i < 60; ++i)
//        printf("%d: %f, %f\n", i, 1-acc_failure[i], prior->prob[i]*acc_failure[i]/(1-prior->prob[i]*(1-acc_failure[i])));
}

