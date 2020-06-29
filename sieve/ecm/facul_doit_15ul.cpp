#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "facul_ecm.h"  // for ecm_15ul
// IWYU pragma: no_include "mpqs.h"       // for mpqs_15ul
// IWYU pragma: no_include "pm1.h"        // for pm1_15ul
// IWYU pragma: no_include "pp1.h"        // for pp1_27_15ul, pp1_65_15ul
#include "modredc_15ul_default.h" // IWYU pragma: keep
#define pm1 pm1_15ul
#define pp1_27 pp1_27_15ul
#define pp1_65 pp1_65_15ul
#define ecm ecm_15ul
#define mpqs mpqs_15ul
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"
