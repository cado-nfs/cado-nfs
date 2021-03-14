set (bin_list
        dup1 dup2 purge merge merge-dl replay replay-dl sm sm_simple reconstructlog-dl
)

set (purge_sources purge_matrix.c singleton_removal.c
                   clique_removal.cpp purge.c)

# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 utils_with_io)

add_executable(dup2 dup2.c)
target_link_libraries(dup2 utils_with_io)

add_executable(filter_galois filter_galois.cpp)
target_link_libraries(filter_galois utils_with_io)

add_executable(purge ${purge_sources})
target_link_libraries(purge utils_with_io)


add_executable(merge merge.c merge_heap.c merge_bookkeeping.c
    merge_compute_weights.c
    read_purgedfile_in_parallel.c mst.c merge_replay_matrix.c sparse.c)
target_link_libraries(merge utils_with_io)
set_target_properties(merge PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")

add_executable(merge-dl merge.c merge_heap.c merge_bookkeeping.c
    merge_compute_weights.c
    read_purgedfile_in_parallel.c mst.c merge_replay_matrix.c sparse.c)
target_link_libraries(merge-dl utils_with_io)
set_target_properties(merge-dl PROPERTIES
  COMPILE_FLAGS "-DFOR_DL ${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils_with_io)

add_executable(replay-dl replay.c sparse.c)
target_link_libraries(replay-dl utils_with_io)
set_target_properties(replay-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(reconstructlog-dl reconstructlog.cpp)
target_link_libraries(reconstructlog-dl utils_with_io)

add_executable(sm_simple sm_simple.cpp)
target_link_libraries(sm_simple utils)

add_executable(singleton_simple singleton_simple.cpp)
target_link_libraries(singleton_simple utils_with_io)

add_executable(shrink_rels shrink_rels.cpp)
target_link_libraries(shrink_rels utils)

if (${WITH_MPI})
    add_definitions(-DWITH_MPI)
endif()
add_executable(sm sm.c)
target_link_libraries(sm utils_with_io)
if (${WITH_MPI})
    set_target_properties(sm PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS}" LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")
else()
    set_target_properties(sm PROPERTIES
        COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")
endif()

add_executable(sm_append sm_append.cpp)
target_link_libraries(sm_append utils)
if (${WITH_MPI})
set_target_properties(sm_append PROPERTIES
  COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS}"
  LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}"
  )
else()
set_target_properties(sm_append PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}"
  LINK_FLAGS "${OpenMP_C_FLAGS}"
  )
endif()



install(TARGETS ${bin_list} RUNTIME DESTINATION ${LIBSUFFIX}/filter)

include(CMakeLists-nodist.txt OPTIONAL)
