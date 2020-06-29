#include "cado.h" // IWYU pragma: keep
// facul_doit.cpp will include these files anyway, we don't need to
// include them here.
// IWYU pragma: no_include "facul_ecm.h"  // for ecm_mpz
// IWYU pragma: no_include "mpqs.h"       // for mpqs_mpz
// IWYU pragma: no_include "pm1.h"        // for pm1_mpz
// IWYU pragma: no_include "pp1.h"        // for pp1_27_mpz, pp1_65_mpz
#include "mod_mpz_default.h" // IWYU pragma: keep
#define pm1 pm1_mpz
#define pp1_27 pp1_27_mpz
#define pp1_65 pp1_65_mpz
#define ecm ecm_mpz
#define mpqs mpqs_mpz
#define FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#include "facul_doit.cpp"
