set (bin_list
        dup1 dup2 purge merge merge-dl replay replay-dl sm sm_simple reconstructlog-dl
)

set (purge_sources purge_matrix.c singleton_removal.c
                   clique_removal.cpp purge.c)

# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(dup1)


add_executable(dup2 dup2.c)
target_link_libraries(dup2 utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(dup2)

add_executable(filter_galois filter_galois.cpp)
target_link_libraries(filter_galois utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(filter_galois)


add_executable(purge ${purge_sources})
target_link_libraries(purge utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(purge)



add_executable(merge merge.c merge_heap.c merge_bookkeeping.c
    merge_compute_weights.cpp
    read_purgedfile_in_parallel.c mst.c merge_replay_matrix.c sparse.c)
CADO_MARK_OPENMP_CAPABLE_TARGET(merge)
target_link_libraries(merge utils_with_io)

add_executable(merge-dl merge.c merge_heap.c merge_bookkeeping.c
    merge_compute_weights.cpp
    read_purgedfile_in_parallel.c mst.c merge_replay_matrix.c sparse.c)
CADO_MARK_OPENMP_CAPABLE_TARGET(merge-dl)
set_property(TARGET merge-dl APPEND PROPERTY COMPILE_DEFINITIONS FOR_DL)
target_link_libraries(merge-dl utils_with_io)

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(replay)


add_executable(replay-dl replay.c sparse.c)
target_link_libraries(replay-dl utils_with_io)
set_property(TARGET replay-dl APPEND PROPERTY COMPILE_DEFINITIONS FOR_DL)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(replay-dl)


add_executable(reconstructlog-dl reconstructlog.cpp)
target_link_libraries(reconstructlog-dl utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(reconstructlog-dl)

add_executable(sm_simple sm_simple.cpp)
target_link_libraries(sm_simple utils)

CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(sm_simple)


add_executable(singleton_simple singleton_simple.cpp)
target_link_libraries(singleton_simple utils_with_io)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(singleton_simple)


add_executable(shrink_rels shrink_rels.cpp)
target_link_libraries(shrink_rels utils)
CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(shrink_rels)

add_executable(sm sm.c)
target_link_libraries(sm utils_with_io)
CADO_MARK_OPENMP_CAPABLE_TARGET(sm)
CADO_MARK_MPI_CAPABLE_TARGET(sm)

add_executable(sm_append sm_append.cpp)
target_link_libraries(sm_append utils)
CADO_MARK_OPENMP_CAPABLE_TARGET(sm_append)
CADO_MARK_MPI_CAPABLE_TARGET(sm_append)

install(TARGETS ${bin_list} RUNTIME DESTINATION ${LIBSUFFIX}/filter)

include(CMakeLists-nodist.txt OPTIONAL)
