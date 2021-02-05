add_executable(check_rels check_rels.c)
target_link_libraries(check_rels utils_with_io)
install(TARGETS check_rels RUNTIME DESTINATION ${LIBSUFFIX}/misc)

if (HAVE_GMPECM)
add_executable(descent_init_Fp descent_init_Fp.c smooth_detect.c)
target_link_libraries(descent_init_Fp utils ecm)
install(TARGETS descent_init_Fp RUNTIME DESTINATION ${LIBSUFFIX}/misc)
endif()

add_executable(debug_renumber EXCLUDE_FROM_ALL debug_renumber.cpp)
target_link_libraries(debug_renumber utils)

add_executable(extra_renumber EXCLUDE_FROM_ALL extra_renumber.cpp)
target_link_libraries(extra_renumber utils_with_io)

add_executable(renumber_convert EXCLUDE_FROM_ALL renumber_convert.cpp)
target_link_libraries(renumber_convert utils_with_io)

add_executable(indexed_relations_to_density EXCLUDE_FROM_ALL indexed_relations_to_density.cpp)
target_link_libraries(indexed_relations_to_density utils_with_io)
set_target_properties(indexed_relations_to_density PROPERTIES
  COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" LINK_FLAGS "${OpenMP_CXX_FLAGS}")

include(CMakeLists-nodist.txt OPTIONAL)
