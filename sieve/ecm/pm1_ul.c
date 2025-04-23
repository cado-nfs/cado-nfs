#include "cado.h" // IWYU pragma: keep
#include "arith/modredc_ul_default.h" // IWYU pragma: keep
#define pm1 pm1_ul
#define pm1_stage1 pm1_stage1_ul
#define pp1_stage2 pp1_stage2_ul
// IWYU pragma: no_include "pm1.h"
// IWYU pragma: no_include "pp1.h"
// scan-headers: skip
#include "pm1.c"        // NOLINT(bugprone-suspicious-include)
