#include "cado.h" // IWYU pragma: keep
#include "arith/modredc_15ul_default.h" // IWYU pragma: keep
#define pm1 pm1_15ul
#define pm1_stage1 pm1_stage1_15ul
#define pp1_stage2 pp1_stage2_15ul
// IWYU pragma: no_include "pm1.h"
// IWYU pragma: no_include "pp1.h"
// scan-headers: skip
#include "pm1.c"        // NOLINT(bugprone-suspicious-include)
