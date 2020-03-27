
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft)
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft/mpfq)

add_library(bblas
        bblas_m4ri.cpp
        debug.cpp
        level2a.cpp
        level2b.cpp
        level3a.cpp
        level3a1.cpp
        level3b.cpp
        level3c.cpp
        level3d.cpp
        level4_echelon.cpp
        level4_gauss.cpp
        level4_lup.cpp
        level4_ple.cpp
        level4_pluq.cpp
        level5.cpp
        perm_matrix.cpp
        gauss.c
)

add_dependencies(bblas gf2x-build)
