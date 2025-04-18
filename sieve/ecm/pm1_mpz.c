#include "cado.h" // IWYU pragma: keep
#include "arith/mod_mpz_default.h" // IWYU pragma: keep
#define pm1 pm1_mpz
#define pm1_stage1 pm1_stage1_mpz
#define pp1_stage2 pp1_stage2_mpz
// IWYU pragma: no_include "pm1.h"
// IWYU pragma: no_include "pp1.h"
// scan-headers: skip
#include "pm1.c"        // NOLINT(bugprone-suspicious-include)
