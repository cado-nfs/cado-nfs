#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include "las-multiobj-globals.hpp"     // IWYU pragma: keep
#include "las-config.hpp" // IWYU pragma: keep

#ifndef SUPPORT_LARGE_Q
const int support_large_q = 0;
#else
const int support_large_q = 1;
#endif
#ifndef DLP_DESCENT
const int dlp_descent = 0;
#else
const int dlp_descent = 1;
#endif
