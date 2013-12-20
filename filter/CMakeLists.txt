set (bin_list
        dup1 dup2 purge merge merge-dl replay replay-dl sm reconstructlog-dl
)

# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 utils_with_io)

add_executable(dup2 dup2.c)
target_link_libraries(dup2 utils_with_io)

add_executable(purge filter_common.c purge.c)
target_link_libraries(purge utils_with_io)

add_executable(merge merge.c filter_common.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge utils_with_io)

add_executable(merge-dl merge.c filter_common.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge-dl utils_with_io)
set_target_properties(merge-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(replay filter_common.c replay.c sparse.c)
target_link_libraries(replay utils_with_io)

add_executable(replay-dl filter_common.c replay.c sparse.c)
target_link_libraries(replay-dl utils_with_io)
set_target_properties(replay-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(reconstructlog-dl reconstructlog.c)
target_link_libraries(reconstructlog-dl utils_with_io)

add_executable(sm sm.c)
target_link_libraries(sm utils_with_io)

add_executable(compute_cols EXCLUDE_FROM_ALL compute_cols.c)
target_link_libraries(compute_cols utils_with_io)

install(TARGETS ${bin_list} compute_cols RUNTIME DESTINATION bin/filter)

include(CMakeLists-nodist.txt OPTIONAL)
