set (SUF_SRC ${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st/)

include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st)
include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/ecm)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve)

cado_define_test(test_fm.cpp ${SUF_SRC}fm.cpp ${SUF_SRC}tab_fm.cpp) 

cado_define_test(test_strategy.cpp
  ${SUF_SRC}tab_strategy.cpp ${SUF_SRC}strategy.cpp
  ${SUF_SRC}fm.cpp ${SUF_SRC}tab_fm.cpp)

cado_define_test(test_decomp.cpp ${SUF_SRC}tab_decomp.cpp
    ${SUF_SRC}decomp.cpp LIBRARIES utils)

cado_define_test(test_convex_hull.cpp
  ${SUF_SRC}convex_hull.cpp
  NO_DEFAULT_RUN
  LIBRARIES utils
)

cado_divert_test(test_convex_hull 10 N=10 seed=1 x0=-20 x1=20)
cado_divert_test(test_convex_hull 256 N=256 seed=1 x0=-1000 x1=1000)
