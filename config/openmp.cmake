
if(MINGW)
    message(STATUS "OpenMP disabled for MINGW")
else()
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (HAVE_OPENMP 1)
        # Don't set unconditionally !
        # set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        # OpenMP.cmake leaves crap around with clang.
        set(OpenMP_C_FLAGS CACHE STRING "C compiler flags for OpenMP parallelization" FORCE)
        set(OpenMP_CXX_FLAGS CACHE STRING "CXX compiler flags for OpenMP parallelization" FORCE)
    endif()
endif()
