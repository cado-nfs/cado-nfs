
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft)
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft/mpfq)

add_library(bblas
        bblas_m4ri.cpp
        bblas_debug.cpp
        bblas_level2a.cpp
        bblas_level2b.cpp
        bblas_level3a.cpp
        bblas_level3a1.cpp
        bblas_level3b.cpp
        bblas_level3c.cpp
        bblas_level3d.cpp
        bblas_level4_echelon.cpp
        bblas_level4_gauss.cpp
        bblas_level4_lup.cpp
        bblas_level4_ple.cpp
        bblas_level4_pluq.cpp
        bblas_level5.cpp
        bblas_perm_matrix.cpp
        bblas_gauss.c
        bblas_mat8.cpp
)

add_dependencies(bblas gf2x-build)
