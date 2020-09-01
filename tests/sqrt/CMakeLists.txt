if (HAVE_CABSL AND HAVE_CLOG)
    # While crtalgsqrt works over mpi, actually testing the code base
    # over mpi is tricky, because the programs seems to auto-choose the
    # mpi setting based on the requested max ram.
    #
    #     set(mpiconfs mpi_square1 mpi_square2)
    #     if (NOT HAVE_MPI)
    #         set(mpiconfs)
    #     endif()
    set(mpiconfs)
    foreach(mpiconf nompi ${mpiconfs})
        set(suffix _${mpiconf})
        if(suffix STREQUAL _nompi)
            set(suffix)
        endif()
        cado_define_test(test_crtalgsqrt${suffix}
            SCRIPT
            ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/test_crtalgsqrt.sh
            ${CADO_NFS_BINARY_DIR}/sqrt/crtalgsqrt
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.dep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.ratdep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.poly
            TARGET_DEPENDENCIES crtalgsqrt)
    endforeach()
endif()
