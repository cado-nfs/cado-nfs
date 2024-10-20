#
# Try to detect the fact that we're running coverage tests.
#
# It seems that some tests, esp. in debug mode, are overly sensitive to
# coverage checking, and tend to take a very long time. By detecting
# coverage at the cmake level, we have the option to special-case these
# tests in coverage mode (e.g., discard them, or make them less
# concurrency-tolerant)
#

message(STATUS "Are we running some kind of coverage testing tool?")
string(REGEX MATCH "--coverage"
    COVERAGE_TEST
    ${CMAKE_C_FLAGS} ${CMAKE_CXX_FLAGS}
    ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}
)

if (COVERAGE_TEST)
    message(STATUS "Are we running some kind of coverage testing tool? -- Yes")
else()
    message(STATUS "Are we running some kind of coverage testing tool? -- No")
endif()
