#include "cado.h" // IWYU pragma: keep

#ifdef TRACE_K
#error "This file *must not* be compiled with TRACE_K defined"
#undef TRACE_K
#endif

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

// IWYU pragma: no_include "las-where-am-i.hpp"
#include <cstdint>                    // for int64_t, uint64_t
#include <memory>                     // for unique_ptr
#include "las-where-am-i-proxy.hpp"   // for where_am_I, where_am_I::pimpl_t
#include "las-where-am-i-prod.hpp"    // for where_am_I::impl
struct cxx_param_list; // IWYU pragma: keep
class nfs_work; // IWYU pragma: keep


int extern_trace_on_spot_ab(int64_t, uint64_t) {
    return 0;
}

/* The ctors / dtors in production code are trivial

where_am_I::where_am_I() : pimpl{ new impl{} } { }
where_am_I::~where_am_I() = default;
where_am_I::where_am_I(where_am_I const & x) : pimpl(new impl(*x.pimpl)) {
}
where_am_I & where_am_I::operator=(where_am_I const & x) {
    *pimpl = *x.pimpl;
    return *this;
}

 */

void where_am_I::decl_usage(cxx_param_list &)
{
}

void where_am_I::interpret_parameters(cxx_param_list &)
{
}

/* This fills all the trace_* structures from the main one. The main
 * structure is the one for which a non-NULL pointer is passed.
 */
void where_am_I::begin_special_q(nfs_work const &)
{
    return;
}

/* }}} */
