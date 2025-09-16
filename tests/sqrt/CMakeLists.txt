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
            ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh
            --handle-mpirun-directly
            "${mpiconf}"
            ${CADO_NFS_BINARY_DIR}/sqrt/crtalgsqrt
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.dep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.ratdep.000
            ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.poly
            ram=${crtalgsqrt_ram${suffix}}
            TARGET_DEPENDENCIES crtalgsqrt
            # either 40462797324737803355716975649 or
            # 44371162641954939938390944368 or both.
            EXPECT_SHA1 e3115858599ff987d979a464e6339595a8fc8a6b,4ff2c48e27890c0977c0f1bd2f1482c7eaac23e8,80abc5bf07b5d24e7055879cd11090ad6442ff08,105dbae2b2b0ed5cb99b4f4ee8e115c223c1bec8
            SHA1_ON_REGEXP_LINES "^# [][0-9\\. ()+-]*$"
            SHA1_ON_SED_OUTPUT "s/^# \\[[0-9\\.]*\\]\\([0-9 ()+-]*\\)\\[\\([0-9]*\\)\\]$/\\2/"
        )
    endforeach()
