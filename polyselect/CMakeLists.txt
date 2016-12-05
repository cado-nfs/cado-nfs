# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.

# Root OPTimization sources.
set (ropt_sources 
  ropt_linear.c ropt_quadratic.c ropt_str.c ropt_param.c
  ropt_tree.c ropt_arith.c ropt_io.c ropt_stage1.c ropt_stage2.c
  )

add_library(polyselect_common ${CADO_SHARED} auxiliary.c rho.c murphyE.c area.c size_optimization.c)
add_executable (polyselect polyselect.c polyselect_str.c polyselect_arith.c ropt.c ${ropt_sources})

add_executable (polyselect_ropt ropt_main.c ropt.c ${ropt_sources})
set_target_properties(polyselect_ropt PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS "${OpenMP_C_FLAGS}")
add_executable (polyselect_gfpn EXCLUDE_FROM_ALL polyselect_gfpn.c gfpkdlpolyselect.c ${ropt_sources})

# the utils library already includes -lgmp and -lm
target_link_libraries (polyselect_gfpn  polyselect_common utils)
target_link_libraries (polyselect_common utils)
target_link_libraries (polyselect polyselect_common utils ${pthread_libs})
target_link_libraries (polyselect_ropt polyselect_common utils ${pthread_libs})

if(ENABLE_SHARED)
    install(TARGETS polyselect_common DESTINATION ${LIBSUFFIX}/polyselect)
endif()
install(TARGETS polyselect RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
install(TARGETS polyselect_ropt RUNTIME DESTINATION ${LIBSUFFIX}/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
