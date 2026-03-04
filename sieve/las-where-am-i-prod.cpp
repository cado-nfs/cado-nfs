#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#ifdef TRACE_K
#error "This file *must not* be compiled with TRACE_K defined"
#endif

#include <cstdint>

#include "las-where-am-i-proxy.hpp"   // for where_am_I, where_am_I::pimpl_t

/* las-where-am-i-prod.hpp defines the (empty) struct where_am_I::impl */
#include "las-where-am-i-prod.hpp" // IWYU pragma: keep


struct cxx_param_list; // IWYU pragma: keep
class nfs_work; // IWYU pragma: keep


int extern_trace_on_spot_ab(cxx_mpz const &, cxx_mpz const &) {
    return 0;
}

int extern_trace_on_spot_ab(int64_t, uint64_t) {
    return 0;
}

/* The ctors / dtors in production code are trivial. We provide external
 * definitions because of !178, though.
 *
 * We do _not_ want the class to be trivially destructible/copyable
 * externally, because there's a minimal runtime and compilation time
 * savings by doing it via external symbol resolutions.
 */

where_am_I::where_am_I() = default;
where_am_I::~where_am_I() = default;    // NOLINT
where_am_I::where_am_I(where_am_I const & x) = default;
where_am_I & where_am_I::operator=(where_am_I const & x) = default;

void where_am_I::decl_usage(cxx_param_list &)
{
}

void where_am_I::interpret_parameters(cxx_param_list &)
{
}

/* This fills all the trace_* structures from the main one. The main
 * structure is the one for which a non-NULL pointer is passed.
 */
void where_am_I::begin_special_q(
        nfs_work const &,
        special_q_data_class auto const &)
{
}

template void where_am_I::begin_special_q(nfs_work const &, qlattice_basis const &);
template void where_am_I::begin_special_q(nfs_work const &, siqs_special_q_data const &);

/* }}} */
