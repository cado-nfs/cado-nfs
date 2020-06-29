#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "facul_ecm.h"  // for ecm_2ul2
// IWYU pragma: no_include "mpqs.h"       // for mpqs_2ul2
// IWYU pragma: no_include "pm1.h"        // for pm1_2ul2
// IWYU pragma: no_include "pp1.h"        // for pp1_27_2ul2, pp1_65_2ul2
#include "modredc_2ul2_default.h" // IWYU pragma: keep
#define pm1 pm1_2ul2
#define pp1_27 pp1_27_2ul2
#define pp1_65 pp1_65_2ul2
#define ecm ecm_2ul2
#define mpqs mpqs_2ul2
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"
