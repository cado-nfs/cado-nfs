include_directories (utils_st)
include_directories (../ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

set (SUF_UTILS utils_st/)
set (SRC_UTILS
    ${SUF_UTILS}convex_hull.cpp ${SUF_UTILS}fm.cpp ${SUF_UTILS}tab_fm.cpp
    ${SUF_UTILS}strategy.cpp ${SUF_UTILS}tab_strategy.cpp)

set (LIB_EXEC facul utils)

add_executable(gfm gfm_exec.cpp generate_factoring_method.cpp ${SRC_UTILS})
target_link_libraries(gfm ${LIB_EXEC})

add_executable(benchfm benchfm_exec.cpp generate_factoring_method.cpp ${SRC_UTILS})
target_link_libraries(benchfm ${LIB_EXEC})

add_executable(gst 
  gst_exec.cpp
  generate_strategies.cpp
  gen_decomp.cpp
  ${SUF_UTILS}/decomp.cpp
  ${SUF_UTILS}/tab_decomp.cpp
  ${SRC_UTILS})

target_link_libraries(gst ${LIB_EXEC})

add_executable(finalst 
  final_exec.cpp
  finding_good_strategy.cpp
  ${SRC_UTILS})

target_link_libraries(finalst ${LIB_EXEC})

add_executable(gen_decomp
  gen_decomp_main.cpp
  gen_decomp.cpp
  ${SUF_UTILS}/decomp.cpp
  ${SUF_UTILS}/tab_decomp.cpp)

target_link_libraries(gen_decomp utils)


include(CMakeLists-nodist.txt OPTIONAL)
