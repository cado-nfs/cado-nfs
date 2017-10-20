set (SUF_SRC ${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st/)

include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st)
include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/ecm)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve)

cado_define_test(test_fm.c ${SUF_SRC}fm.c ${SUF_SRC}tab_fm.cpp) 

cado_define_test(test_strategy.c
  ${SUF_SRC}tab_strategy.c ${SUF_SRC}strategy.c
  ${SUF_SRC}fm.c ${SUF_SRC}tab_fm.cpp)

set(more)
if(MINGW)
    set(more LIBRARIES gnurx)
endif()
cado_define_test(test_decomp.c ${SUF_SRC}tab_decomp.c ${SUF_SRC}decomp.c ${more})

cado_define_test(test_convex_hull.c
  ${SUF_SRC}convex_hull.c
  ${SUF_SRC}point.c
  ${SUF_SRC}tab_point.c
  NO_DEFAULT_RUN
)

cado_divert_test(test_convex_hull 1)
cado_divert_test(test_convex_hull 2 1000 0.001)
cado_divert_test(test_convex_hull 3 100000 0.00001)
cado_divert_test(test_convex_hull 4 100000 -0.00001)
