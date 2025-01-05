
include(FindOpenMP)

if (OPENMP_FOUND)
    set (HAVE_OPENMP 1)
    message(STATUS "OpenMP is enabled")
    message(STATUS "C flags for OpenMP: ${OpenMP_C_FLAGS}")
    message(STATUS "C libraries for OpenMP: ${OpenMP_C_LIBRARIES}")
    message(STATUS "C++ flags for OpenMP: ${OpenMP_CXX_FLAGS}")
    message(STATUS "C++ libraries for OpenMP: ${OpenMP_CXX_LIBRARIES}")
    # FindOpenMP produces a space-separated list... If more than one
    # option is required, this will get in our way.
    separate_arguments(OpenMP_C_FLAGS)
    separate_arguments(OpenMP_CXX_FLAGS)
    separate_arguments(OpenMP_C_LIBRARIES)
    separate_arguments(OpenMP_CXX_LIBRARIES)
else()
    # OpenMP.cmake leaves crap around with clang.
    set(OpenMP_C_FLAGS CACHE STRING "C compiler flags for OpenMP parallelization" FORCE)
    set(OpenMP_CXX_FLAGS CACHE STRING "CXX compiler flags for OpenMP parallelization" FORCE)
    message(STATUS "OpenMP is disabled")
endif()

macro(mark_targets_as_openmp)
    if(OPENMP_FOUND)
        foreach(t ${ARGN})
            set_property(TARGET ${t} APPEND PROPERTY COMPILE_OPTIONS ${OpenMP_C_FLAGS})
            set_property(TARGET ${t} APPEND PROPERTY LINK_OPTIONS ${OpenMP_C_FLAGS})
            # get_target_property(tt ${t} TYPE)
            # if (target_type STREQUAL "LIBRARY")
            #     target_link_libraries(${t} ${OpenMP_C_LIBRARIES})
            # else()
            #     set_property(TARGET ${t} APPEND PROPERTY LIBRARIES ${OpenMP_C_LIBRARIES})
            # endif()
            target_link_libraries(${t} ${OpenMP_C_LIBRARIES})
        endforeach()
    endif()
endmacro()

