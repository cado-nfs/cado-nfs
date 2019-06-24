# posix_memalign is buggy on openbsd-59-amd64
# both with gcc 4.2.1 and egcc (based on gcc 4.9.3)
# (see https://marc.info/?l=openbsd-bugs&m=155377511718586)
  
include(CheckCSourceRuns)

set(check_posix_memalign_code "
#include <stdlib.h>
int main()
{
    void *res = NULL;
    for (size_t i = 1; i <= 2340; i++)
        {
            int rc = posix_memalign (&res, 64, i);
            if (rc != 0 || (((unsigned long) res) % 64) != 0)
                __builtin_abort ();
        }
    return 0;
}
")

CHECK_C_SOURCE_RUNS("${check_posix_memalign_code}" HAVE_POSIX_MEMALIGN)

if(HAVE_POSIX_MEMALIGN)
    message(STATUS "Testing whether posix_memalign exists and works -- yes")
else()
    message(STATUS "Testing whether posix_memalign exists and works -- no")
endif()
