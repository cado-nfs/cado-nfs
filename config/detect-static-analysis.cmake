#
# There are several ways to run static analysis on the cado-nfs code.
#
# The first one is done automatically in the git pipelines with the
# coverity branch, which run cado-nfs through the synopsis coverity tool.
#
# The second one is to work locally with llvm scan-build, e.g. as:
#
# DEBUG=1 CLANG=1 scan-build make cmake
# DEBUG=1 CLANG=1 scan-build make -j8 utils
#
# We also know of the CodeChecker tool: https://codechecker.readthedocs.io/en/latest/
#
# In both cases, the STATIC_ANALYSIS config flag is set in cado_config.h,
# so that asserts such that ASSERT_FOR_STATIC_ANALYZER are expanded.

message(STATUS "Are we running some kind of static analysis tool?")
string(REGEX MATCH "scan-build" STATIC_ANALYSIS ${CMAKE_C_COMPILER})
if(NOT STATIC_ANALYSIS)
    include(CheckCSourceCompiles)
    set(OLD_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
    set(CMAKE_REQUIRED_QUIET 1)
    CHECK_C_SOURCE_COMPILES("
#if !defined(__COVERITY__)
    choke me
#endif
    int main() { return 0; }
    " RUNNING_COVERITY)
    if (RUNNING_COVERITY)
        set(STATIC_ANALYSIS "coverity")
    endif()
    set(CMAKE_REQUIRED_QUIET "${OLD_CMAKE_REQUIRED_QUIET}")
endif()
if(NOT STATIC_ANALYSIS)
    string(COMPARE NOTEQUAL "$ENV{CC_DATA_FILES_DIR}" "" RUNNING_CODECHECKER)
    if (RUNNING_CODECHECKER)
        set(STATIC_ANALYSIS "CodeChecker")
    endif()
endif()

if(STATIC_ANALYSIS)
    message(STATUS "Are we running some kind of static analysis tool? -- Yes, ${STATIC_ANALYSIS}")
else()
    message(STATUS "Are we running some kind of static analysis tool? -- No")
endif()
