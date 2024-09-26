
add_subdirectory(bblas)

include_directories(bwc) # for worker-threads.[ch]

include_directories(bblas)

# the bitlinalg library will almost surely go away at some point, or
# possibly be incorporated in bblas
add_library(bitlinalg blockmatrix.cpp)
target_link_libraries(bitlinalg bblas)

if(NOT HAVE_GF2X)
add_dependencies(bitlinalg gf2x-build)
endif()


add_executable(characters characters.cpp bwc/worker-threads.c)
target_link_libraries(characters bitlinalg utils_with_io pthread)
install(TARGETS characters RUNTIME DESTINATION ${LIBSUFFIX}/linalg)

if(HAVE_M4RI)
    set(M4RI_C_FLAGS ${CMAKE_C_FLAGS})
    set(M4RI_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    set(M4RI_LDFLAGS ${CMAKE_SHARED_LINKER_FLAGS})
    set(M4RI_CONFIGURE_FLAGS --disable-dependency-tracking CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${M4RI_C_FLAGS} CXXFLAGS=${M4RI_CXX_FLAGS} LDFLAGS=${M4RI_LDFLAGS} $ENV{M4RI_CONFIGURE_EXTRA_FLAGS})
    get_override(M4RI_CONFIGURE_FLAGS M4RI_CONFIGURE_FLAGS)
    add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/linalg/m4ri/Makefile
        COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/linalg/m4ri
        COMMAND echo Configuring linalg/m4ri with options ${M4RI_CONFIGURE_FLAGS}
        COMMAND cd m4ri && ${CADO_NFS_SOURCE_DIR}/linalg/m4ri/configure ${M4RI_CONFIGURE_FLAGS} < /dev/null
        COMMENT "Configuring linalg/m4ri"
        )
    add_custom_target(m4ri-build ALL
        DEPENDS ${CADO_NFS_BINARY_DIR}/linalg/m4ri/Makefile
        WORKING_DIRECTORY m4ri
        COMMAND "${CMAKE_MAKE_PROGRAM}"
        COMMENT "Building linalg/m4ri"
        )
else()
    message(STATUS "Found linalg/m4ri source tree, but no ./configure there. Run autoreconf ?")
endif()

if (HAVE_M4RIE)
    set(M4RIE_C_FLAGS ${CMAKE_C_FLAGS})
    set(M4RIE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    set(M4RIE_LDFLAGS ${CMAKE_SHARED_LINKER_FLAGS})
    set(M4RIE_CONFIGURE_FLAGS --disable-dependency-tracking CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${M4RIE_C_FLAGS} CXXFLAGS=${M4RIE_CXX_FLAGS} LDFLAGS=${M4RIE_LDFLAGS} $ENV{M4RIE_CONFIGURE_EXTRA_FLAGS})
    get_override(M4RIE_CONFIGURE_FLAGS M4RIE_CONFIGURE_FLAGS)
    add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/linalg/m4rie/Makefile
        COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/linalg/m4rie
        COMMAND echo Configuring linalg/m4rie with options ${M4RIE_CONFIGURE_FLAGS}
        COMMAND cd m4rie && ${CADO_NFS_SOURCE_DIR}/linalg/m4rie/configure ${M4RIE_CONFIGURE_FLAGS} < /dev/null
        COMMENT "Configuring linalg/m4rie"
        )
    add_custom_target(m4rie-build ALL
        DEPENDS ${CADO_NFS_BINARY_DIR}/linalg/m4rie/Makefile
        WORKING_DIRECTORY m4rie
        COMMAND "${CMAKE_MAKE_PROGRAM}"
        COMMENT "Building linalg/m4rie"
        )
endif()


include(CMakeLists-nodist.txt OPTIONAL)
