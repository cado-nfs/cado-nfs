
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/bblas)

if(HAVE_M4RI)
link_directories (${CADO_NFS_BINARY_DIR}/linalg/m4ri/.libs)
include_directories(${CADO_NFS_BINARY_DIR}/linalg/m4ri)
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/m4ri)
endif()
if(HAVE_M4RIE)
link_directories (${CADO_NFS_BINARY_DIR}/linalg/m4rie/.libs)
include_directories(${CADO_NFS_BINARY_DIR}/linalg/m4rie)
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/m4rie)
endif()

set(test_bblas_sources
        test_bblas.cpp
        test_bblas_base.cpp
        test_bblas_level2.cpp
        test_bblas_level3.cpp
        test_bpack.cpp
        test_bblas_level4_pluq.cpp
        test_bblas_level4_ple.cpp
        test_bblas_level4_gauss.cpp
        test_bblas_level5.cpp
    )
if(HAVE_M4RI)
    set(test_bblas_sources ${test_bblas_sources} test_bblas_m4ri.cpp)
endif()
add_executable(test_bblas EXCLUDE_FROM_ALL ${test_bblas_sources})
target_link_libraries(test_bblas bblas utils)
if(HAVE_M4RI)
add_dependencies(test_bblas m4ri-build)
target_link_libraries(test_bblas m4ri)
endif()
if(HAVE_M4RIE)
add_dependencies(test_bblas m4rie-build)
target_link_libraries(test_bblas m4rie)
endif()
set_target_properties(test_bblas PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries(test_bblas "${OpenMP_CXX_FLAGS}")

cado_define_test(test_bblas_level2 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l2
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_bblas_level3 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l3
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

# easily times out in coverage mode
cado_define_test(test_bblas_level4 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l4
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas
    TIMEOUT 120
)

cado_define_test(test_bblas_level5 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l5
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_gauss.c
    LIBRARIES bblas ${gmp_libname})
