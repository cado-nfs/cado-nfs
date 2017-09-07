# test if rand() gives true alea
message(STATUS "Testing if rand() gives true alea")
try_run(test-rand_runs test-rand_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/test-rand.c)
if(test-rand_runs EQUAL 0)
   set(HAVE_RAND_BUG 0)
   message(STATUS "Testing if rand() gives true alea -- Probably")
else()
   set(HAVE_RAND_BUG 1)
   message(STATUS "Testing if rand() gives true alea -- Buggy")
endif()

message(STATUS "Testing if srand() yields a deterministic sequence")
include(CheckCSourceRuns)
CHECK_C_SOURCE_RUNS("
    #include <stdlib.h>
    int main()
    {
        srand(1);
        unsigned long a = rand();
        srand(1);
        unsigned long b = rand();
        return a == b ? EXIT_SUCCESS : EXIT_FAILURE;
    }
" SRAND_CHECK)

if (SRAND_CHECK_EXITCODE EQUAL 0)
    message(STATUS "Testing if srand() yields a deterministic sequence -- Yes")
    set(HAVE_USUAL_SRAND_DETERMINISTIC_BEHAVIOR 1)
else()
    # From "man srand" in openbsd:
    #
    #   Standards insist that this interface return deterministic results.
    #   Unsafe usage is very common, so OpenBSD changed the subsystem to return
    #   non-deterministic results by default.
    #
    message(STATUS "Testing if srand() yields a deterministic sequence -- No")
    set(HAVE_USUAL_SRAND_DETERMINISTIC_BEHAVIOR 0)
    CHECK_FUNCTION_EXISTS(srand_deterministic	HAVE_SRAND_DETERMINISTIC)
    if(NOT HAVE_SRAND_DETERMINISTIC)
        message(WARNING "tests in cado-nfs will not be reproducible because you do not have deterministic rand()")
    endif()
endif()
