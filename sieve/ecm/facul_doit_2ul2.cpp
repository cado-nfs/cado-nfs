#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "facul_ecm.h"  // for ecm_2ul2
// IWYU pragma: no_include "mpqs.h"       // for mpqs_2ul2
#include "arith/modredc_2ul2_default.h" // IWYU pragma: keep
#define ecm ecm_2ul2
#define mpqs mpqs_2ul2
// scan-headers: skip
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"       // NOLINT(bugprone-suspicious-include)
