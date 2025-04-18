#ifndef CADO_LAS_MULTIOBJ_GLOBALS_HPP
#define CADO_LAS_MULTIOBJ_GLOBALS_HPP

/* This interface provides link time resolution of what would otherwise
 * be compile-time constants. The goal is to avoid recompilation in cases
 * where it can be avoided at zero cost.
 */

extern const int support_large_q;
extern const int dlp_descent;

#endif	/* CADO_LAS_MULTIOBJ_GLOBALS_HPP */
