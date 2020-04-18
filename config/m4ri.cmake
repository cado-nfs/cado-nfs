if (IS_DIRECTORY ${CADO_NFS_SOURCE_DIR}/linalg/m4ri)
    if (EXISTS ${CADO_NFS_SOURCE_DIR}/linalg/m4ri/configure)
        message(STATUS "Found linalg/m4ri source tree, enabling external component")
        set(HAVE_M4RI 1)
    else()
        message(STATUS "Found linalg/m4ri source tree, but no ./configure there. Run autoreconf ?")
    endif()
endif()

if (HAVE_M4RI)
    if (IS_DIRECTORY ${CADO_NFS_SOURCE_DIR}/linalg/m4rie)
        if (EXISTS ${CADO_NFS_SOURCE_DIR}/linalg/m4rie/configure)
            message(STATUS "Found linalg/m4rie source tree, enabling external component")
            set(HAVE_M4RIE 1)
        else()
            message(STATUS "Found linalg/m4rie source tree, but no ./configure there. Run autoreconf ?")
        endif()
    endif()
endif()



