#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "priors.h"

// Warning: priors are supposed to be precomputed and statically stored,
// so that they don't have to be free.
// This is currently not the case -> memory leaks.
prior_ptr get_prior(int lo, int hi, int fbb)
{
    prior_ptr p;
    p = (prior_ptr) malloc (sizeof(prior_struct));
    assert (p != NULL);

    p->cofac_range[0] = lo;
    p->cofac_range[1] = hi;
    p->fbb = fbb;
    for (int i = 0; i <= fbb; ++i)
        p->prob1[i] = 0.0;
    for (int i = fbb+1; i < 60; ++i) 
        // Heurisitic formula for the moment...
        p->prob1[i] = 0.5*powf((float)i, -0.8);
    for (int i = 0; i < 60; ++i) {
        p->prob5[i] = p->prob1[i];
        p->prob7[i] = p->prob1[i];
        p->prob11[i] = p->prob1[i];
    }
    return p;
}



