#include "cado.h" // IWYU pragma: keep
#include "tests_common.h"

int main(int argc, char const * argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_clear();
}
