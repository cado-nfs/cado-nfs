add_executable(check_rels check_rels.c)
target_link_libraries(check_rels utils_with_io)
install(TARGETS check_rels RUNTIME DESTINATION ${LIBSUFFIX}/misc)

if (HAVE_GMPECM)
add_executable(descent_init_Fp descent_init_Fp.cpp smooth_detect.cpp)
target_link_libraries(descent_init_Fp utils ecm)
install(TARGETS descent_init_Fp RUNTIME DESTINATION ${LIBSUFFIX}/misc)
endif()

add_executable(debug_renumber EXCLUDE_FROM_ALL debug_renumber.cpp)
target_link_libraries(debug_renumber utils)

add_executable(explain_indexed_relation EXCLUDE_FROM_ALL
    explain_indexed_relation.cpp)
target_link_libraries(explain_indexed_relation utils)

add_executable(indexed_relations_to_density EXCLUDE_FROM_ALL indexed_relations_to_density.cpp)
target_link_libraries(indexed_relations_to_density utils_with_io)

mark_targets_as_openmp(indexed_relations_to_density)

if(HAVE_GMPECM)
    add_executable(factor factor.cpp)
    target_link_libraries(factor utils ecm)
endif()

add_executable(cl_structure_from_group_order cl_structure_from_group_order.cpp)
target_link_libraries(cl_structure_from_group_order utils)

include(CMakeLists-nodist.txt OPTIONAL)
