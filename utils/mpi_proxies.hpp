#ifndef CADO_MPI_PROXIES_HPP
#define CADO_MPI_PROXIES_HPP

#include "cado_config.h"

#include <cstddef>

#include <vector>
#include <array>
#include <type_traits>

#include "select_mpi.h"

/* here we have stuff that should ease mpi-sending and receiving of some
 * standard containers. We only do fairly easy stuff at the moment.
 *
 * Note that of course there's boost::mpi, so the point of this file
 * becomes completely moot if it ever grows beyond trivial size.
 */

namespace cado_mpi {

template<typename T> struct type_tag {};
template<> struct type_tag<unsigned int> { static constexpr int value = MPI_UNSIGNED; };
#ifndef    UNSIGNED_LONG_IS_EXACTLY_UNSIGNED
template<> struct type_tag<unsigned long> { static constexpr int value = MPI_UNSIGNED_LONG; };
#endif
#ifndef    UNSIGNED_LONG_LONG_IS_EXACTLY_UNSIGNED_LONG
template<> struct type_tag<unsigned long long> { static constexpr int value = MPI_UNSIGNED_LONG_LONG; };
#endif
#if !defined(UINT32_T_IS_EXACTLY_UNSIGNED) && !defined(UINT32_T_IS_EXACTLY_UNSIGNED_LONG)
template<> struct type_tag<uint32_t> { static constexpr int value = UINT32_T; };
#endif
#if !defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG) && !defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG_LONG)
template<> struct type_tag<uint64_t> { static constexpr int value = UINT64_T; };
#endif
template<> struct type_tag<int> { static constexpr int value = MPI_INT; };
#ifndef    LONG_IS_EXACTLY_INT
template<> struct type_tag<long> { static constexpr int value = MPI_LONG; };
#endif
#ifndef    LONG_LONG_IS_EXACTLY_LONG
template<> struct type_tag<long long> { static constexpr int value = MPI_LONG_LONG; };
#endif
#if !defined(INT32_T_IS_EXACTLY_INT) && !defined(INT32_T_IS_EXACTLY_LONG)
template<> struct type_tag<int32_t> { static constexpr int value = INT32_T; };
#endif
#if !defined(INT64_T_IS_EXACTLY_LONG) && !defined(INT64_T_IS_EXACTLY_LONG_LONG)
template<> struct type_tag<int64_t> { static constexpr int value = INT64_T; };
#endif
/* we might want to add more aliases, but pay attention to the fact that
 * we must have unambiguous resolution of the template structs. See also
 * config/check_types.cmake and cado_config.h.in
 */

/*
template<typename T>
void recv(T & ps, int mpi_peer, int tag, MPI_Comm comm);
template<typename T>
void send(T const & ps, int mpi_root, int tag, MPI_Comm comm);
*/

template<typename T>
void recv(std::vector<T> & ps, int mpi_peer, int tag, MPI_Comm comm)
requires std::is_scalar_v<T>
{
    constexpr int mpi_type_tag = type_tag<T>::value;
    MPI_Status status;
    MPI_Probe(mpi_peer, tag, comm, &status);
    int count;
    MPI_Get_count(&status, mpi_type_tag, &count);
    size_t oldsize = ps.size();
    ps.insert(ps.end(), count, T());
    auto tail = ps.data() + oldsize;
    MPI_Recv(tail, count, mpi_type_tag,
            status.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE);
}

template<typename T, std::size_t N>
void
recv(std::vector<std::array<T, N>> & ps, int mpi_peer, int tag, MPI_Comm comm)
requires std::is_scalar_v<T>
{
    constexpr int mpi_type_tag = type_tag<T>::value;
    MPI_Status status;
    MPI_Probe(mpi_peer, tag, comm, &status);
    int count;
    MPI_Get_count(&status, mpi_type_tag, &count);
    size_t oldsize = ps.size();
    ps.insert(ps.end(), count/N, std::array<T, N>());
    auto tail = ps.data() + oldsize;
    MPI_Recv(tail, count, mpi_type_tag,
            status.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE);
}

template<typename T>
void
send(std::vector<T> const & ps, int mpi_root, int tag, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Send(ps.data(),
            ps.size(),
            type_tag<T>::value,
            mpi_root, tag, comm);
}

template<typename T, std::size_t N>
void
send(std::vector<std::array<T, N>> const & ps, int mpi_root, int tag, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Send(ps.data(),
            N * ps.size(),
            type_tag<T>::value,
            mpi_root, tag, comm);
}

template<typename T>
void
allgather(std::vector<T> const & in, std::vector<T> & out, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    int size;
    MPI_Comm_size(comm, &size);
    out.assign(in.size() * size, T());
    MPI_Allgather(in.data(), in.size(), type_tag<T>::value,
                  out.data(), in.size(), type_tag<T>::value,
                  comm);
}

} // end of namespace cado_mpi


#endif	/* CADO_MPI_PROXIES_HPP */
