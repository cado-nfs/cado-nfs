
#############################################################
# mpi
# 20260816: use newer cmake's find_package.

# We want to retain the feature
# of reacting to the MPI env variable, though. As of now, the following
# should work:
#       MPI=1 , MPI=on, MPI=yes, etc.
#       MPI=/opt/openmpi-1.7/
# On the other hand, selecting the MPI flavor on alternatives-based
# systems is not supported, so the following does not work.
#       MPI=openmpi

if($ENV{MPI})
    set(findprog_flags)
    if("$ENV{MPI}" MATCHES "/")
        set(MPI_HOME $ENV{MPI})
        list(APPEND findprog_flags
            HINTS "$ENV{MPI}" "$ENV{MPI}/bin"
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)
    elseif("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
        # do nothing, we'll just rely on find_package below
    else()
        message(FATAL_ERROR "Selecting MPI alternative with MPI=$ENV{MPI} is no longer supported")
        # otherwise we make the .<variant> binary names higher
        # priority than others This is for finding things such as
        # mpicc.mpich2 which get installed by the alternatives
        # mechanism on debian-like systems.
        list(APPEND mpicc_names "mpicc.${MPI}")
        list(APPEND mpicxx_names "mpic++.${MPI}" "mpicxx.${MPI}" "mpiCC.${MPI}")
        list(APPEND mpiexec_names "mpiexec.${MPI}")
        # Well. Presently we're in fact *not* pushing the standard
        # names in the search list. Should we ?
    endif()

    find_package(MPI ${findprog_flags} REQUIRED)

    set(MPIEXEC ${MPIEXEC_EXECUTABLE})
    message(STATUS "Using MPI version ${MPI_VERSION}")
    message(STATUS "Using MPI C compiler ${MPI_C_COMPILER}")
    message(STATUS "Using MPI C++ compiler ${MPI_CXX_COMPILER}")
    message(STATUS "Using MPI driver ${MPIEXEC}")
    get_filename_component(HAVE_MPI ${MPIEXEC} PATH)
    # We're using this variable in the top-level substitution, so it needs
    # to escape its scope and go into the cache right now.
    set(WITH_MPI 1 CACHE INTERNAL "MPI is being used (for relevant code parts)")


    # for intel MPI, we used to add these to the command line.
    # set(MPI_C_COMPILER_CMDLINE_INSERTIONS "-cc=${CMAKE_C_COMPILER}")
    # set(MPI_CXX_COMPILER_CMDLINE_INSERTIONS "-cxx=${CMAKE_CXX_COMPILER}")

    # openmpi up until 1.10 should have the following.
    # message(STATUS "Enabling workaround for long-standing OpenMPI breakage (ompi/pull/1495)")
    # # tl;dr leave_pinned is just plain broken
    # # throughout most of the 1.7, 1.8. 1.9, and 1.10
    # # series of OpenMPI. The work to fix this is at
    # # https://github.com/open-mpi/ompi/pull/1495 ;
    # # see the attached commit logs (namely, commits
    # # 57035744 and 4b7cd1c0 in ompi-release carry the
    # # fix. Those are open-mpi/ompi@7aa03d66 and
    # # open-mpi/ompi@11e2d788 in the ompi repository).
    # set(MPIEXEC_EXTRA_STANZAS "--mca mpi_leave_pinned 0")
else()
    message(STATUS "MPI is not enabled")
    set(WITH_MPI CACHE INTERNAL "MPI is not used")
endif()

macro(mark_targets_as_mpi)
    if(WITH_MPI)
        foreach(t ${ARGN})
            set_property(TARGET ${t} PROPERTY CADO_USES_MPI ON)
            target_link_libraries(${t} PRIVATE MPI::MPI_C MPI::MPI_CXX)
            target_compile_definitions(${t} PRIVATE WITH_MPI)
        endforeach()
    endif()
endmacro()
