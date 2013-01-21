#ifndef __METHODS_H__
#define __METHODS_H__

// Select the best method for the given state, according to some naive
// notion. 
cofac_method_srcptr get_method_naive(cofac_state_srcptr state);

// Compute the new state after the given method is run. This assumes the
// method failed. Otherwise, the prior and the cofac_size fields must be
// changed.
void update_state(cofac_state_ptr newstate, cofac_state_srcptr state,
        cofac_method_srcptr method);

#endif   /* __METHODS_H__ */
