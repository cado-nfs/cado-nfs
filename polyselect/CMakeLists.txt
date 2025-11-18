# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.

add_library(polyselect_common
    ${CADO_SHARED}
    auxiliary.c
    murphyE.cpp
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
mark_targets_as_openmp(polyselect_middle)

add_executable (polyselect
    polyselect.c
    $<TARGET_OBJECTS:polyselect_middle>
    )
target_link_libraries (polyselect polyselect_common ${pthread_libs})
mark_targets_as_openmp(polyselect)

add_executable (polyselect_ropt
    ropt_main.c
    best_polynomials_queue.cpp
    ropt_single_sublattice_priority_queue.cpp
    ropt_sublattice_priority_queue.cpp
    ropt_sublattice_crt.cpp
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
    ropt_sublattice_priority_queue.cpp
)
target_link_libraries (polyselect_ropt polyselect_common ${pthread_libs})
mark_targets_as_openmp(polyselect_ropt)


add_executable (polyselect_gfpn polyselect_gfpn.c gfpkdlpolyselect.c)
target_link_libraries (polyselect_gfpn  polyselect_common utils)

add_executable (polyselect3 polyselect3.c)
mark_targets_as_openmp(polyselect3)
target_link_libraries (polyselect3 polyselect_common utils)

if(ENABLE_SHARED)
    install(TARGETS polyselect_common DESTINATION ${LIBSUFFIX}/polyselect)
endif()
install(TARGETS polyselect RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
install(TARGETS polyselect_ropt RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
install(TARGETS polyselect3 RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)

# dlpolyselect
add_executable (dlpolyselect dlpolyselect.c)
mark_targets_as_openmp(dlpolyselect)
target_link_libraries (dlpolyselect polyselect_common utils)
install(TARGETS dlpolyselect RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)


include(CMakeLists-nodist.txt OPTIONAL)

# skewness
add_executable (skewness skewness.c)
target_link_libraries (skewness polyselect_common utils)
install(TARGETS skewness RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)

# score
add_executable (score score.c)
target_link_libraries (score polyselect_common utils)
install(TARGETS score RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
