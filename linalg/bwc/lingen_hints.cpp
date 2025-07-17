#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <ostream>
#include <istream>
#include <type_traits>
#include <utility>
#include <vector>

#include "lingen_hints.hpp"
#include "select_mpi.h"

template <typename T>
static void
share(T & m, int root, MPI_Comm comm)
    requires std::is_trivially_copyable_v<typename T::mapped_type>
{
    using K = typename T::key_type;
    using M = typename T::mapped_type;
    using V = std::pair<K, M>;

    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<V> serial;
    serial.insert(serial.end(), m.begin(), m.end());
    unsigned long ns = serial.size();
    MPI_Bcast(&ns, 1, MPI_UNSIGNED_LONG, root, comm);
    if (rank)
        serial.assign(ns, V());
    MPI_Bcast(serial.data(), ns * sizeof(V), MPI_BYTE, 0, comm);
    if (rank)
        m.insert(serial.begin(), serial.end());
}

void lingen_hints::share(int root, MPI_Comm comm)
{
    ::share((super &)*this, root, comm);
    MPI_Bcast(&tt_gather_per_unit, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&tt_scatter_per_unit, 1, MPI_DOUBLE, root, comm);
}

std::istream & lingen_hints::unserialize(std::istream & is)
{
    size_t n;
    for (is >> n; is && n--;) {
        key_type K;
        is >> K;
        auto it = find(K);
        if (!is.good())
            return is;
        if (it != end()) {
            is >> it->second;
        } else {
            /* This is very dangerous, as M is not complete at this point
             * */
            mapped_type M;
            is >> M;
            M.complete = false;
            (*this)[K] = M;
        }
    }
    return is;
}

std::ostream & lingen_hints::serialize(std::ostream & os) const
{
    os << size() << "\n";
    for (auto const & x: *this) {
        os << x.first << " " << x.second << "\n";
    }
    return os;
}
