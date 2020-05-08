#include "cado.h" // IWYU pragma: keep
#include "las-multiobj-globals.hpp"     // IWYU pragma: keep
#include "las-config.h" // IWYU pragma: keep

#ifndef SUPPORT_LARGE_Q
int support_large_q = 0;
#else
int support_large_q = 1;
#endif
#ifndef DLP_DESCENT
int dlp_descent = 0;
#else
int dlp_descent = 1;
#endif
