#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "facul_ecm.h"  // for ecm_mpz
// IWYU pragma: no_include "mpqs.h"       // for mpqs_mpz
// IWYU pragma: no_include "pm1.h"        // for pm1_mpz
// IWYU pragma: no_include "pp1.h"        // for pp1_27_mpz, pp1_65_mpz
#include "modredc_ul_default.h" // IWYU pragma: keep
#define pm1 pm1_ul
#define pp1_27 pp1_27_ul
#define pp1_65 pp1_65_ul
#define ecm ecm_ul
#define mpqs mpqs_ul
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"
