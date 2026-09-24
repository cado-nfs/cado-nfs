#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "mpqs.h"       // for mpqs_15ul
#include "arith/modredc_15ul_default.h" // IWYU pragma: keep
#define mpqs mpqs_15ul
// scan-headers: skip
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"       // NOLINT(bugprone-suspicious-include)
