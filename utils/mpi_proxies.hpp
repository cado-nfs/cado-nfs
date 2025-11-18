#ifndef CADO_MPI_PROXIES_HPP
#define CADO_MPI_PROXIES_HPP

#include "cado_config.h"

#include <cstddef>
#include <cstdlib>

#include <vector>
#include <array>
#include <type_traits>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "select_mpi.h"
#include "runtime_numeric_cast.hpp"

/* here we have stuff that should ease mpi-sending and receiving of some
 * standard containers. We only do fairly easy stuff at the moment.
 *
 * Note that of course there's boost::mpi, so the point of this file
 * becomes completely moot if it ever grows beyond trivial size.
 */

namespace cado_mpi {

    /* We'd like the following to be constexpr, but they can't because
     * openmpi defines MPI_UNSIGNED and friends as results of dirty casts
     * from void*, which can't be constexpr.
     */
template<typename T> struct type_tag {};
template<> struct type_tag<unsigned int> { static MPI_Datatype value() { return MPI_UNSIGNED; } };
#ifndef    UNSIGNED_LONG_IS_EXACTLY_UNSIGNED
template<> struct type_tag<unsigned long> { static MPI_Datatype value() { return MPI_UNSIGNED_LONG; } };
#endif
#ifndef    UNSIGNED_LONG_LONG_IS_EXACTLY_UNSIGNED_LONG
template<> struct type_tag<unsigned long long> { static MPI_Datatype value() { return MPI_UNSIGNED_LONG_LONG; } };
#endif
#if !defined(UINT32_T_IS_EXACTLY_UNSIGNED) && !defined(UINT32_T_IS_EXACTLY_UNSIGNED_LONG)
template<> struct type_tag<uint32_t> { static MPI_Datatype value() { return UINT32_T; } };
#endif
#if !defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG) && !defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG_LONG)
template<> struct type_tag<uint64_t> { static MPI_Datatype value() { return UINT64_T; } };
#endif
template<> struct type_tag<int> { static MPI_Datatype value() { return MPI_INT; } };
#ifndef    LONG_IS_EXACTLY_INT
template<> struct type_tag<long> { static MPI_Datatype value() { return MPI_LONG; } };
#endif
#ifndef    LONG_LONG_IS_EXACTLY_LONG
template<> struct type_tag<long long> { static MPI_Datatype value() { return MPI_LONG_LONG; } };
#endif
#if !defined(INT32_T_IS_EXACTLY_INT) && !defined(INT32_T_IS_EXACTLY_LONG)
template<> struct type_tag<int32_t> { static MPI_Datatype value() { return INT32_T; } };
#endif
#if !defined(INT64_T_IS_EXACTLY_LONG) && !defined(INT64_T_IS_EXACTLY_LONG_LONG)
template<> struct type_tag<int64_t> { static MPI_Datatype value() { return INT64_T; } };
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
    constexpr MPI_Datatype mpi_type_tag = type_tag<T>::value();
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
    constexpr MPI_Datatype mpi_type_tag = type_tag<T>::value();
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
            type_tag<T>::value(),
            mpi_root, tag, comm);
}

template<typename T, std::size_t N>
void
send(std::vector<std::array<T, N>> const & ps, int mpi_root, int tag, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Send(ps.data(),
            N * ps.size(),
            type_tag<T>::value(),
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
    MPI_Allgather(in.data(), in.size(), type_tag<T>::value(),
                  out.data(), in.size(), type_tag<T>::value(),
                  comm);
}

template<typename T>
void
allreduce(T & in, MPI_Op op, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Allreduce(MPI_IN_PLACE, &in, 1, type_tag<T>::value(), op, comm);
}

template<typename T>
void
broadcast(T & in, int root, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Bcast(&in, 1, type_tag<T>::value(), root, comm);
}

template<typename T>
void
broadcast(T * in, size_t n, int root, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    MPI_Bcast(in, runtime_numeric_cast<int>(n), type_tag<T>::value(), root, comm);
}

static inline void
broadcast(mpz_ptr z, int root, MPI_Comm comm)
{
    int k;
    MPI_Comm_rank(comm, &k);
    static_assert(std::is_signed_v<mp_size_t>);
    /* it is very important that nlimbs is actuall z->_mp_size (with the
     * sign bit), otherwise broadcasting would erase the sign!
     */
    mp_size_t nlimbs = mpz_size(z) * mpz_sgn(z);
    broadcast(nlimbs, root, comm);
    mp_limb_t * ptr = mpz_limbs_write(z, std::abs(nlimbs));
    broadcast(ptr, std::abs(nlimbs), root, comm);
    mpz_limbs_finish(z, nlimbs);
    // broadcast(z->_mp_size, root, comm);
}

static inline void
broadcast(cxx_mpz & z, int root, MPI_Comm comm)
{
    broadcast((mpz_ptr) z, root, comm);
}



template<typename T>
bool
mpi_data_agrees(T const & in, MPI_Comm comm)
    requires std::is_scalar_v<T>
{
    /* in principle it should also be possible to do min and max */
    T b_and, b_or;
    // NOLINTBEGIN(google-readability-casting)
    MPI_Allreduce((void *) (&in), (void *) (&b_and), sizeof(T), MPI_BYTE, MPI_BAND, comm);
    MPI_Allreduce((void *) (&in), (void *) (&b_or),  sizeof(T), MPI_BYTE, MPI_BOR, comm);
    // NOLINTEND(google-readability-casting)
    return b_and == b_or;
}/**/

} // end of namespace cado_mpi


#endif	/* CADO_MPI_PROXIES_HPP */
