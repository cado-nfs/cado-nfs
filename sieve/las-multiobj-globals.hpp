#ifndef LAS_MULTIOBJ_GLOBALS_HPP_
#define LAS_MULTIOBJ_GLOBALS_HPP_

/* This interface provides link time resolution of what would otherwise
 * be compile-time constants. The goal is to avoid recompilation in cases
 * where it can be avoided at zero cost.
 */

extern int support_large_q;
extern int dlp_descent;

#endif	/* LAS_MULTIOBJ_GLOBALS_HPP_ */
