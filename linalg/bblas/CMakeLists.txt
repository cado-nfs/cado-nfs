
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft)
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft/mpfq)

add_library(bblas
        bblas_debug.cpp
        bblas_level2a.cpp
        bblas_level2b.cpp
        bblas_level3a.cpp
        bblas_level3a1.cpp
        bblas_level3b.cpp
        bblas_level3c.cpp
        bblas_level3d.cpp
        bpack.cpp
        bblas_level4_echelon.cpp
        bblas_level4_gauss.cpp
        bblas_level4_lup.cpp
        bblas_level4_pluq.cpp
        bblas_level5.cpp
        bblas_perm_matrix.cpp
        bblas_gauss.c
        bblas_mat8.cpp
        bblas_mat64.cpp
)

set_target_properties(bblas PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries(bblas "${OpenMP_CXX_FLAGS}")

add_dependencies(bblas gf2x-build)
