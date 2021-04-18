#include "cado.h" // IWYU pragma: keep
#include <vector>
#include <ostream>       // for istream, operator<<, ostream, basic_ostream
#include <type_traits>   // for enable_if, is_trivially_copyable
#include <utility>       // for pair, move
#include "select_mpi.h"
#include "lingen_hints.hpp"

template<typename T>
    typename std::enable_if<std::is_trivially_copyable<typename T::mapped_type>::value, void>::type
share(T & m, int root, MPI_Comm comm)
{
    typedef typename T::key_type K;
    typedef typename T::mapped_type M;
    typedef std::pair<K, M> V;

    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<V> serial;
    serial.insert(serial.end(), m.begin(), m.end());
    unsigned long ns = serial.size();
    MPI_Bcast(&ns, 1, MPI_UNSIGNED_LONG, root, comm);
    if (rank) serial.assign(ns, V());
    MPI_Bcast(serial.data(), ns * sizeof(V), MPI_BYTE, 0, comm);
    if (rank) m.insert(serial.begin(), serial.end());
}

void lingen_hints::share(int root, MPI_Comm comm)
{
    ::share((super&) *this, root, comm);
    MPI_Bcast(&tt_gather_per_unit, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&tt_scatter_per_unit, 1, MPI_DOUBLE, root, comm);
}

std::istream& lingen_hints::unserialize(std::istream& is) {
    size_t n;
    for(is >> n ; is && n-- ; ) {
        key_type K;
        is >> K;
        auto it = find(K);
        if (!is.good()) return is;
        if (it != end()) {
            is >> it->second;
        } else {
            /* This is very dangerous, as M is not complete at this point
             * */
            mapped_type M;
            is >> M;
            M.complete = false;
            (*this)[K] = std::move(M);
        }
    }
    return is;
}

std::ostream& lingen_hints::serialize(std::ostream& os) const {
    os << size() << "\n";
    for(auto const & x : *this) {
        os << x.first << " " << x.second << "\n";
    }
    return os;
}
