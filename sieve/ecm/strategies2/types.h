#ifndef __TYPES_H__
#define __TYPES_H__

//////////////////////////////////////////////////////////////////////////////
// The type for histograms of probabilities.
// Typically, the i-th coefficient is the probability related to primes
// with i-bits.

typedef float histogram_t [60];


//////////////////////////////////////////////////////////////////////////////
// Cofactorization methods
// ECM, p-1, p+1
// We distinguish 2 variants of p+1 depending on their starting point,
// and having different properties related to congruences of p mod 12.
// See Kruppa's PhD thesis.

#define PM1 0
#define PP1_27 1
#define PP1_65 2
#define ECM 3

typedef struct {
    int type;           // ecm, p-1, p+1
    int B1;
    int B2;
    float ms[3];            // millisecs, for i words input
    histogram_t success1;   // proba to find a prime == 1 mod 12
    histogram_t success5;   // proba to find a prime    5
    histogram_t success7;   // proba to find a prime    7
    histogram_t success11;  // proba to find a prime    11
} cofac_method_struct;

typedef cofac_method_struct  cofac_method_t[1];
typedef cofac_method_struct *cofac_method_ptr;
typedef const cofac_method_struct *cofac_method_srcptr;


//////////////////////////////////////////////////////////////////////////////
// P-1 and P+1 are more or less "one-shot". Redoing it with the same
// parameters is useless. However, it can make sense to redo them with
// a higher B1. Furthermore, they have a complicated interactions with
// each others, due to congruence properties mod 12.
// The following type allows to keep track of what happened before.

typedef struct {
    // Success probability of previously run p-1, for each congruence
    // class of prime. These are just pointers, since they will be copies
    // of data that come from a method.
    const float *pm1_success1;
    const float *pm1_success5;
    const float *pm1_success7;
    const float *pm1_success11;
    // The same for p+1. Here we talk about "genuine" p+1, i.e. p+1 for
    // which we know that the starting point is such that the
    // discriminant of the quadratic extension is not a square, so that
    // we are not redoing p-1. Hence, running a p+1 instance will
    // contribute to some of the data below, but also to some of the p-1
    // data.
    const float *pp1_success1;
    const float *pp1_success5;
    const float *pp1_success7;
    const float *pp1_success11;
} ppm1_history_struct;

typedef       ppm1_history_struct  ppm1_history_t[1];
typedef       ppm1_history_struct *ppm1_history_ptr;
typedef const ppm1_history_struct *ppm1_history_srcptr;


//////////////////////////////////////////////////////////////////////////////
// Current state in the cofactorization chain.
// In principle, the method to be chosen to continue depends only on the
// data within this structure (and some global parameters).

typedef struct {
    int cofac_size[2];            // A range of bit-sizes for the cofactor.
    const float *prior;           // Probabilities of presence of a prime
                                  // at the beginning.
    ppm1_history_t ppm1_history;  // What happened before with p+1 and p-1.
    histogram_t acc_failure1;     // Accumulated probabilities of failure
    histogram_t acc_failure5;     // of previous methods for each congruence
    histogram_t acc_failure7;     // class mod 12.
    histogram_t acc_failure11;
} cofac_state_struct;

typedef       cofac_state_struct  cofac_state_t[1];
typedef       cofac_state_struct *cofac_state_ptr;
typedef const cofac_state_struct *cofac_state_srcptr;


//////////////////////////////////////////////////////////////////////////////
// Node in the strategy tree.
// A strategy is a tree of nodes. Each node corresponds to a state (in
// the above sense), and a method to be run. There is also pointers to
// child nodes that correspond to how to continue the cofactorization
// depending on the outcome of running the method.
//
// The number of childs is bounded by a constant, to simplify our life
#define MAX_CHILD 10

// Child types:
#define CHILD_ABORT      0       // Case where proba of being smooth is
                                 // considered to small to continue.
#define CHILD_NOT_SMOOTH 1       // Case where we found a prime cofac that is
                                 // too large.
#define CHILD_SMOOTH     2       // Case where we have finished.
#define CHILD_CONTINUE   3       // Other cases, we continue the cofac strategy.

// Forward declaration of a pointer type:
struct __node_struct;
typedef struct __node_struct node_struct;

// Helper structure for child pointers.
typedef struct {
    int type;
    float proba;           // Probability that we go to that child.
    node_struct * node;    // NULL pointer unless type == CHILD_CONTINUE
} child_pointer_struct;

// The main node structure.
struct __node_struct {
    cofac_state_t state;         // State when entering the node
    cofac_method_srcptr method;  // Method to run within this node
    
    // Where to go, depending on the outcome of running the method:
    int n_childs;
    child_pointer_struct childs[MAX_CHILD];

    // Performance data for the subtree rooted at this node.
    float ms;                    // Expected runtime of the subtree.
    float success_proba;         // Success probability of the subtree.
};

typedef       node_struct  node_t[1];
typedef       node_struct *node_ptr;
typedef const node_struct *node_srcptr;





#endif   /* __TYPES_H__ */
