    # While crtalgsqrt works over mpi, actually testing the code base
    # over mpi is tricky, because the programs seems to auto-choose the
    # mpi setting based on the requested max ram.
    #
    #     set(mpiconfs mpi_square1 mpi_square2)
    #     if (NOT HAVE_MPI)
    #         set(mpiconfs)
    #     endif()
    set(mpiconfs)
    if(HAVE_MPI)
        list(APPEND mpiconfs mpi_rect2)
    endif()
    set(crtalgsqrt_ram 3)
    set(crtalgsqrt_ram_mpi_rect2 0.000004)
    foreach(mpiconf nompi ${mpiconfs})
        set(suffix _${mpiconf})
        if(suffix STREQUAL _nompi)
            set(suffix)
        endif()
        cado_define_test(test_crtalgsqrt${suffix}
            SCRIPT
            ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/test_crtalgsqrt.sh
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.dep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.ratdep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.poly
            ram=${crtalgsqrt_ram${suffix}}
            TARGET_DEPENDENCIES crtalgsqrt)
    endforeach()
