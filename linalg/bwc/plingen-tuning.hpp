#ifndef PLINGEN_TUNING_HPP_
#define PLINGEN_TUNING_HPP_

#include "params.h"
#include "plingen.hpp"
#include "select_mpi.h"
#include <map>								\

/* This object is passed as a companion info to a call
 * of bw_biglingen_recursive ; it is computed by the code in
 * plingen-tuning.cpp
 * but once tuning is over, it is essentially fixed.
 */

struct lingen_substep_scheduling_characteristics {/*{{{*/
    /* output characteristics -- the ones we have to choose */

    /* batch, shrink0, shrink2. Used to be static parameters, now they're
     * dynamic.
     *
     * batch: with batch=b, schedule b times more simultaneous
     *  transforms, so that we can have more parallelism. Costs b times
     *  more RAM
     * shrink0: divide the number of transforms on dimension 0 (m/n for
     *  MP, (m+n)/r for MUL) so as to save memory
     * shrink2: divide the number of transforms on dimension 2 (m+n)/r
     *  for both MP and MUL) so as to save memory
     */

    /* shrink0 and shrink2 are important. Here, we restrict our operation
     * to shrink0*shrink2 multiplications of matrices of size
     * (n0/shrink0) times (n2/shrink2).
     */
    unsigned int shrink0 = 1, shrink2 = 1;

    /* batch is not used much, because it goes in the direction of
     * spending _more_ memory, which is rarely what we're interested in.
     * What this does is that several steps of the outer loop (off size
     * n1/r) are done simultaneously: more transforms are kept live.
     */
    unsigned int batch = 1;

};/*}}}*/

struct lingen_call_companion {
    bool recurse;
    unsigned int weight;
    double ttb;
    struct {
        lingen_substep_scheduling_characteristics S;
        double tt;
        size_t ram;
    } mp, mul;
    struct key {
        int depth;
        size_t L;
        bool operator<(key const& a) const {
            if (depth < a.depth) return true;
            if (depth > a.depth) return false;
            return L < a.L;
        }
    };
};

typedef std::map<lingen_call_companion::key, lingen_call_companion> lingen_hints_t;

void plingen_tuning_decl_usage(cxx_param_list & pl);
void plingen_tuning_lookup_parameters(cxx_param_list & pl);
lingen_hints_t plingen_tuning(dims * d, MPI_Comm comm, cxx_param_list & pl);

#endif	/* PLINGEN_TUNING_HPP_ */
