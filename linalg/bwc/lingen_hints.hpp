#ifndef CADO_LINGEN_HINTS_HPP
#define CADO_LINGEN_HINTS_HPP

#include <cstddef>
#include <map>
#include <iosfwd>

#include "select_mpi.h"
#include "lingen_call_companion.hpp"

struct lingen_hints : public std::map<lingen_call_companion::key, lingen_call_companion> {
    typedef std::map<lingen_call_companion::key, lingen_call_companion> super;
    double tt_scatter_per_unit;
    double tt_gather_per_unit;
    int ipeak;
    size_t peak;
    void share(int root, MPI_Comm comm);
    std::istream& unserialize(std::istream& is);
    std::ostream& serialize(std::ostream& os) const;
};

inline std::ostream& operator<<(std::ostream& os, lingen_hints const & c) {
    return c.serialize(os);
}

inline std::istream& operator>>(std::istream& is, lingen_hints & c) {
    return c.unserialize(is);
}



#endif	/* LINGEN_HINTS_HPP_ */
