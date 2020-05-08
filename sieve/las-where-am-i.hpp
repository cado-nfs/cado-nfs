
/* This header file is a proxy to the ehader that defines the real
 * implementation of the where_am_I structure. It varies depending on the
 * compilation flag TRACE_K (as well as TRACK_CODE_PATH, but IMO we
 * should only have only one debug flavor).
 */

#ifdef TRACE_K
#include "las-where-am-i-debug.hpp" // IWYU pragma: keep
#else
#include "las-where-am-i-prod.hpp"  // IWYU pragma: keep
#endif
