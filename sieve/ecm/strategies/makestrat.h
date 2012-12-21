#ifndef __MAKESTRAT_H__
#define __MAKESTRAT_H__

// name of method types
#define PM1 0
#define PP1_27 1
#define PP1_65 2
#define ECM 3

typedef struct {
    int type;           // ecm, p-1, p+1
    int B1;
    int B2;
    int param;          // sigma, x0
    float ms[3];        // millisecs, for i words input
    float success1[60];  // proba to find an i bits prime
    float success5[60];  // proba to find an i bits prime
    float success7[60];  // proba to find an i bits prime
    float success11[60];  // proba to find an i bits prime
} cofac_method_struct;

typedef cofac_method_struct  cofac_method_t[1];
typedef cofac_method_struct *cofac_method_ptr;
typedef const cofac_method_struct *cofac_method_srcptr;


// A "prior" is the set of probabilities for the existence of a prime of
// a given length dividing the cofactor. 
// It depends on the size of the cofactor, and of the factor base bound
// (no prime below a certain size). 
// If we have run several unsuccessful methods, we can deduce the
// "knowledge", i.e. the probablity of existence of a prime factor of the
// given length by applying the Bayesian theorem to the prior and the
// accumulated failure probability of the previous methods.
// The prior changes when we find a factor.

typedef struct {
    int cofac_range[2]; // it is for cofac with sizes in this [.,.[ range;
    int fbb;            // and with factor base bound fbb.
    float prob1[60];     // proba that there is a prime factor of i bits.
    float prob5[60];     // proba that there is a prime factor of i bits.
    float prob7[60];     // proba that there is a prime factor of i bits.
    float prob11[60];     // proba that there is a prime factor of i bits.
} prior_struct;

typedef prior_struct  prior_t[1];
typedef prior_struct *prior_ptr;
typedef const prior_struct *prior_srcptr;

typedef struct {
    const float *pm1_success1;
    const float *pm1_success5;
    const float *pm1_success7;
    const float *pm1_success11;
    const float *pp1_success1;
    const float *pp1_success5;
    const float *pp1_success7;
    const float *pp1_success11;
} ppm1_history_struct;

typedef       ppm1_history_struct  ppm1_history_t[1];
typedef       ppm1_history_struct *ppm1_history_ptr;
typedef const ppm1_history_struct *ppm1_history_srcptr;


#endif   /* __MAKESTRAT_H__ */
