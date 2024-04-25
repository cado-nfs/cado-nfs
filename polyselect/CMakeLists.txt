# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.

add_library(polyselect_common
    ${CADO_SHARED}
    auxiliary.c
    murphyE.c
    area.c
    size_optimization.c
    polyselect_alpha.c
    polyselect_norms.c
    )
target_link_libraries (polyselect_common utils)

# It _seems_ that several polyselect-like libraries share the same
# general working principle, and that these objects can be treated as a
# bundle. (we want them to be compiled only once, for efficiency)
add_library(polyselect_middle
    OBJECT
    polyselect_arith.c
    polyselect_collisions.c
    polyselect_data_series.c
    polyselect_special_q.c
    polyselect_thread.c
    polyselect_thread_team.c
    polyselect_thread_league.c
    polyselect_main_data.c
    polyselect_main_queue.c
    polyselect_match.c
    polyselect_poly_header.c
    polyselect_primes_table.c
    polyselect_priority_queue.c
    polyselect_proots.c
    polyselect_qroots.c
    polyselect_shash.c
    polyselect_stats.c
)
set_target_properties(polyselect_middle PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")

add_executable (polyselect
    polyselect.c
    $<TARGET_OBJECTS:polyselect_middle>
    )
target_link_libraries (polyselect polyselect_common ${pthread_libs})
set_target_properties(polyselect PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")

add_executable (polyselect_ropt
    ropt_main.c
    ropt.c
    ropt_arith.c
    ropt_io.c
    ropt_linear.c
    ropt_param.c
    ropt_quadratic.c
    ropt_stage1.c
    ropt_stage2.c
    ropt_str.c
    ropt_tree.c
)
target_link_libraries (polyselect_ropt polyselect_common ${pthread_libs})
set_target_properties(polyselect_ropt PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")


add_executable (polyselect_gfpn polyselect_gfpn.c gfpkdlpolyselect.c)
target_link_libraries (polyselect_gfpn  polyselect_common utils)

add_executable (polyselect3 polyselect3.c)
set_target_properties(polyselect3 PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")
target_link_libraries (polyselect3 polyselect_common utils)

if(ENABLE_SHARED)
    install(TARGETS polyselect_common DESTINATION ${LIBSUFFIX}/polyselect)
endif()
install(TARGETS polyselect RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
install(TARGETS polyselect_ropt RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
install(TARGETS polyselect3 RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)

# dlpolyselect
add_executable (dlpolyselect dlpolyselect.c)
set_target_properties(dlpolyselect PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")
target_link_libraries (dlpolyselect polyselect_common utils)
install(TARGETS dlpolyselect RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)

# skewness
add_executable (skewness skewness.c)
target_link_libraries (skewness polyselect_common utils)
install(TARGETS skewness RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
