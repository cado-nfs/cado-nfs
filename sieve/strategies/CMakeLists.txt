include_directories (utils_st)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

set (SUF_UTILS utils_st/)
set (SRC_UTILS ${SUF_UTILS}point.c ${SUF_UTILS}tab_point.c
    ${SUF_UTILS}convex_hull.c ${SUF_UTILS}fm.c ${SUF_UTILS}tab_fm.cpp
    ${SUF_UTILS}strategy.c ${SUF_UTILS}tab_strategy.c)

set (LIB_EXEC m facul utils)

add_executable(gfm gfm_exec.cpp generate_factoring_method.cpp ${SRC_UTILS})
target_link_libraries(gfm ${LIB_EXEC})

add_executable(benchfm benchfm_exec.cpp generate_factoring_method.cpp ${SRC_UTILS})
target_link_libraries(benchfm ${LIB_EXEC})

add_executable(gst 
  gst_exec.cpp
  generate_strategies.cpp
  gen_decomp.c
  ${SUF_UTILS}/decomp.c
  ${SUF_UTILS}/tab_decomp.c
  ${SRC_UTILS})

target_link_libraries(gst ${LIB_EXEC})

add_executable(finalst 
  final_exec.c
  finding_good_strategy.cpp
  ${SRC_UTILS})

target_link_libraries(finalst ${LIB_EXEC})

add_executable(gen_decomp
  test_gen_decomp.c
  gen_decomp.c
  ${SUF_UTILS}/decomp.c
  ${SUF_UTILS}/tab_decomp.c)

if(MINGW)
    target_link_libraries(gen_decomp gnurx m)
else()
    target_link_libraries(gen_decomp m)
endif()


include(CMakeLists-nodist.txt OPTIONAL)
