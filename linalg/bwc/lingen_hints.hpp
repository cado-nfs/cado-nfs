#ifndef LINGEN_HINTS_HPP_
#define LINGEN_HINTS_HPP_

#include <map>
#include "select_mpi.h"
#include "lingen_call_companion.hpp"

struct lingen_hints : public std::map<lingen_call_companion::key, lingen_call_companion> {
    typedef std::map<lingen_call_companion::key, lingen_call_companion> super;
    double tt_scatter_per_unit;
    double tt_gather_per_unit;
    int ipeak;
    size_t peak;
    void share(int root, MPI_Comm comm);
};


#endif	/* LINGEN_HINTS_HPP_ */
