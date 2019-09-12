#include "cado.h"
#include "select_mpi.h"
#include "lingen_hints.hpp"
#include <vector>

#if defined(__GNUC__) && ! GNUC_VERSION_ATLEAST(5, 0, 0)
/* ugly workaround. g++ < 5 has no standards-conforming is_trivially_copyable
 * https://stackoverflow.com/questions/25123458/is-trivially-copyable-is-not-a-member-of-std
 *
 * (I think that it's the only use in the cado-nfs source tree thus far,
 * so it's fine to have it here. Otherwise it would escalate to cado.h,
 * probably)
 */
namespace std {
    template<typename T>
    struct is_trivially_copyable {
        static constexpr const bool value = __has_trivial_copy(T);
    };
}
#endif

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
    MPI_Bcast(&serial.front(), ns * sizeof(V), MPI_BYTE, 0, comm);
    if (rank) m.insert(serial.begin(), serial.end());
}

void lingen_hints::share(int root, MPI_Comm comm)
{
    ::share((super&) *this, root, comm);
    MPI_Bcast(&tt_gather_per_unit, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&tt_scatter_per_unit, 1, MPI_DOUBLE, root, comm);
}

