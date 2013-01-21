#include "types.h"

// The list of available methods is stored in an external file, that is
// automatically generated with a script that run benchs.
#define NB_METHOD 796
static const cofac_method_t method_list[NB_METHOD] = {
#include "list_methods.i"
};


// Select the best method for the given state, according to some naive
// notion. 
cofac_method_srcptr get_method_naive(cofac_state_srcptr state);

// Compute the new state after the given method is run. This assumes the
// method failed. Otherwise, the prior and the cofac_size fields must be
// changed.
void update_state(cofac_state_ptr newstate, cofac_state_srcptr state,
        cofac_method_srcptr method);


// Let x be the proba that there exist a prime of a given size.
// Let y be the proba that all previously run methods failed to find such
// a prime.
// Then the proba that there still exist a prime of this size is given by
// the Bayesian Theorem, hence the following formula.
static inline bayes(float x, float y)
{
    return x*y/(1-x*(1-y));
}


static const histogram_t zero_sub = {0,}; // fill in with zeros.


cofac_method_srcptr get_method_naive(cofac_state_srcptr state);
{
    int best_meth = -1;
    float best_res = 0;
    int nb_bits = state->cofac_size[1]-1;
    int nb_words = 1 + ((nb_bits-1) / 64);

    histogram_t prob1, prob5, prob7, prob11;

    // Get the probability of presence of a factor of a given size:
    for (int i = 0; i < 60; ++i) {
        prob1[i] = bayes(state->prior[i], state->acc_failure1[i]);
        prob5[i] = bayes(state->prior[i], state->acc_failure5[i]);
        prob7[i] = bayes(state->prior[i], state->acc_failure7[i]);
        prob11[i] = bayes(state->prior[i], state->acc_failure11[i]);
    }

    // Evaluate all available method: for each one, compute the
    // probability of finding a factor.
    for (int i = 0; i < NB_METHOD; ++i) {
        // For p-1 and p+1, one has to subtract the success probabilities
        // of previously run p-1 and p+1.
        const float *ppm1_sub1 = zero_sub;
        const float *ppm1_sub5 = zero_sub;
        const float *ppm1_sub7 = zero_sub;
        const float *ppm1_sub11= zero_sub;
        if (method_list[i]->type == PM1) {
            ppm1_sub1  = state->ppm1_history->pm1_success1;
            ppm1_sub5  = state->ppm1_history->pm1_success5;
            ppm1_sub7  = state->ppm1_history->pm1_success7;
            ppm1_sub11 = state->ppm1_history->pm1_success11;
        }
        if (method_list[i]->type == PP1_27) {
            ppm1_sub1  = state->ppm1_history->pm1_success1;
            ppm1_sub5  = state->ppm1_history->pp1_success5;
            ppm1_sub7  = state->ppm1_history->pm1_success7;
            ppm1_sub11 = state->ppm1_history->pp1_success11;
        }
        if (method_list[i]->type == PP1_65) {
            ppm1_sub1  = state->ppm1_history->pm1_success1;
            ppm1_sub5  = state->ppm1_history->pm1_success5;
            ppm1_sub7  = state->ppm1_history->pp1_success7;
            ppm1_sub11 = state->ppm1_history->pp1_success11;
        }

        // We can now get the success probability of the method:
        float res = 0.0;
        for (int j = 0; j < 60; ++j) {
            float r;
            r  = prob1[j]  * (method_list[i]->success1[j] - ppm1_sub1[j]);
            r += prob5[j]  * (method_list[i]->success5[j] - ppm1_sub5[j]);
            r += prob7[j]  * (method_list[i]->success7[j] - ppm1_sub7[j]);
            r += prob11[j] * (method_list[i]->success11[j]- ppm1_sub11[j]);
            r /= 4;
            res += r;
        }
        // We divide this probability by the cost of running the method,
        // and keep the best according to this naive measure.
        res /= method_list[i]->ms[nb_words];
        if ((best_meth == -1) || (res > best_res)) {
            best_meth = i;
            best_res = res;
        }
    }

    return &method_list[best_meth][0];
}



const float * max_prob(const float *p1, const float *p2) {
    float s1, s2;
    for(int i = 0; i < 60; ++i) {
        s1 += p1[i];
        s2 += p2[i];
    }
    return (s1 < s2) ? p2 : p1;
}

void update_state(cofac_state_ptr newstate, cofac_state_srcptr state,
        cofac_method_srcptr method)
{
  // We assume failure, so that cofac_size and prior are unchanged.
  // In case of success of the methods, these are the two fields to be
  // changed:
  newstate->cofac_size[0] = state->cofac_size[0];
  newstate->cofac_size[1] = state->cofac_size[1];
  newstate->prior = state->prior;

  // Accumulated failure is just the product of the failure
  // probabilities.
  for (int i = 0; i < 60; ++i) {
    newstate->acc_failure1[i] =(1-method->success1[i])*state->acc_failure1[i];
    newstate->acc_failure5[i] =(1-method->success5[i])*state->acc_failure5[i];
    newstate->acc_failure7[i] =(1-method->success7[i])*state->acc_failure7[i];
    newstate->acc_failure11[i]=(1-method->success11[i])*state->acc_failure11[i];
  }

  // Update the history for p-1 and p+1.
  // First copy the previous one:
  newstate->ppm1_history->pm1_success1 = state->ppm1_history->pm1_success1;
  newstate->ppm1_history->pm1_success5 = state->ppm1_history->pm1_success5;
  newstate->ppm1_history->pm1_success7 = state->ppm1_history->pm1_success7;
  newstate->ppm1_history->pm1_success11= state->ppm1_history->pm1_success11;
  newstate->ppm1_history->pp1_success1 = state->ppm1_history->pp1_success1;
  newstate->ppm1_history->pp1_success5 = state->ppm1_history->pp1_success5;
  newstate->ppm1_history->pp1_success7 = state->ppm1_history->pp1_success7;
  newstate->ppm1_history->pp1_success11= state->ppm1_history->pp1_success11;
  // Then update:
  if (method->type == PM1) {
    newstate->ppm1_history->pm1_success1 = method->success1;
    newstate->ppm1_history->pm1_success5 = method->success5;
    newstate->ppm1_history->pm1_success7 = method->success7;
    newstate->ppm1_history->pm1_success11= method->success11;
  } else if (method->type == PP1_27) {
    newstate->ppm1_history->pm1_success1 = max_prob(method->success1, state->ppm1_history->pm1_success1);
    newstate->ppm1_history->pp1_success5 = max_prob(method->success5, state->ppm1_history->pp1_success5);
    newstate->ppm1_history->pm1_success7 = max_prob(method->success7, state->ppm1_history->pm1_success7);
    newstate->ppm1_history->pp1_success11= max_prob(method->success11,state->ppm1_history->pp1_success11);
  } else if (method->type == PP1_65) {
    newstate->ppm1_history->pm1_success1 = max_prob(method->success1, state->ppm1_history->pm1_success1);
    newstate->ppm1_history->pm1_success5 = max_prob(method->success5, state->ppm1_history->pm1_success5);
    newstate->ppm1_history->pp1_success7 = max_prob(method->success7, state->ppm1_history->pp1_success7);
    newstate->ppm1_history->pp1_success11= max_prob(method->success11,state->ppm1_history->pp1_success11);
  } 
}


