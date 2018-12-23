
#############################################################
# mpi
# Don't use the FindMPI module, it's buggy.

# I can't make if($ENV{MPI}) evaluate to true as I want. In particular, I
# want all of these yield true:
# MPI=1 , MPI=on, MPI=yes, etc.
# MPI=/opt/openmpi-1.7/
# MPI=openmpi
# The following excerpt from the doc seems to be just plain wrong:
# #         if(variable)
# #
# #       True if the variable's value is not empty, 0, N, NO, OFF, FALSE,
# #       NOTFOUND, or <variable>-NOTFOUND.
if("$ENV{MPI}" MATCHES "^(0|NO|no|OFF|off|)$")
    message(STATUS "MPI is not enabled")
    set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
    set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
else()
    set(findprog_flags)
    set(mpicc_names)
    set(mpicxx_names)
    set(mpiexec_names)
    if("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
        list(APPEND mpicc_names "mpicc")
        list(APPEND mpicxx_names "mpic++" "mpicxx" "mpiCC")
        list(APPEND mpiexec_names "mpiexec")
    else()
        if("$ENV{MPI}" MATCHES "/")
            # If MPI contains a /, then we assume it should be a path
            list(APPEND findprog_flags
                HINTS "$ENV{MPI}" "$ENV{MPI}/bin"
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)
            # These are the standard names.
            list(APPEND mpicc_names "mpicc")
            list(APPEND mpicxx_names "mpic++" "mpicxx" "mpiCC")
            list(APPEND mpiexec_names "mpiexec")
        else()
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
    endif()

    if(DEFINED ENV{MPI_C_COMPILER})
        set(MPI_C_COMPILER "$ENV{MPI_C_COMPILER}")
    else()
        find_program(MPI_C_COMPILER NAMES ${mpicc_names} NAMES_PER_DIR ${findprog_flags})
    endif()


    if(DEFINED ENV{MPI_CXX_COMPILER})
        set(MPI_CXX_COMPILER "$ENV{MPI_CXX_COMPILER}")
    else()
        find_program(MPI_CXX_COMPILER NAMES ${mpicxx_names} NAMES_PER_DIR ${findprog_flags})
    endif()

    if(DEFINED ENV{MPIEXEC})
        set(MPIEXEC "$ENV{MPIEXEC}")
    else()
        find_program(MPIEXEC NAMES ${mpiexec_names} NAMES_PER_DIR ${findprog_flags})
    endif()

    if (MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(STATUS "Using MPI C compiler ${MPI_C_COMPILER}")
        message(STATUS "Using MPI C++ compiler ${MPI_CXX_COMPILER}")
        message(STATUS "Using MPI driver ${MPIEXEC}")
        get_filename_component(HAVE_MPI ${MPIEXEC} PATH)
        # We're using this variable in the top-level substitution, so it needs
        # to escape its scope and go into the cache right now.
        set(WITH_MPI 1 CACHE INTERNAL "MPI is being used (for relevant code parts)")


        # Run mpicc -v to detect the MPI implementation. We need this at least
        # to add some command-line arguments to MPI builds for the Intel MPI
        # implementation.

        execute_process(COMMAND ${MPI_C_COMPILER} -v
                RESULT_VARIABLE test_return_code
                OUTPUT_VARIABLE test_stdout
                ERROR_VARIABLE test_stderr
                )
        string(STRIP "${test_stdout}" test_stdout)
        if (test_return_code)
        else()
            if(test_stdout MATCHES "mpi.*for MVAPICH2 version (.*)")
                message(STATUS "MPI C Compiler is mvapich2, version ${CMAKE_MATCH_1}")
                set(MPI_COMPILER_IS_MVAPICH2 1)
                set(MPI_MVAPICH2_COMPILER_VERSION ${CMAKE_MATCH_1})
            elseif(test_stdout MATCHES "mpi.*for MPICH.*version (.*)")
                message(STATUS "MPI C Compiler is mpich, version ${CMAKE_MATCH_1}")
                set(MPI_COMPILER_IS_MPICH 1)
                set(MPI_MPICH_COMPILER_VERSION ${CMAKE_MATCH_1})
            elseif(test_stdout MATCHES "^mpi.*for.*Intel.*MPI.*Library ([0-9].*) for")
                message(STATUS "MPI C Compiler is Intel MPI, version ${CMAKE_MATCH_1}")
                set(MPI_COMPILER_IS_INTEL_MPI 1)
                set(MPI_INTEL_COMPILER_VERSION ${CMAKE_MATCH_1})
		set(MPI_C_COMPILER_CMDLINE_INSERTIONS "-cc=${CMAKE_C_COMPILER}")
		set(MPI_CXX_COMPILER_CMDLINE_INSERTIONS "-cxx=${CMAKE_CXX_COMPILER}")
            else()
                # perhaps it's openmpi, but openmpi won't tell on mere
                # mpicc -v...
                execute_process(COMMAND ${MPI_C_COMPILER} "-showme:version"
                        RESULT_VARIABLE test_return_code
                        OUTPUT_VARIABLE test_stdout
                        ERROR_VARIABLE test_stderr
                        )
                string(STRIP "${test_stdout} ${test_stderr}" test_stdout)
                if(test_stdout MATCHES "Open MPI ([^ ]*)")
                    message(STATUS "MPI C Compiler is Open MPI, version ${CMAKE_MATCH_1}")
                    set(MPI_COMPILER_IS_OPEN_MPI 1)
                    set(MPI_OPEN_MPI_COMPILER_VERSION ${CMAKE_MATCH_1})
                    if(MPI_OPEN_MPI_COMPILER_VERSION VERSION_GREATER 1.6.5
                            AND
                            MPI_OPEN_MPI_COMPILER_VERSION VERSION_LESS 2.0.0)
                        message(STATUS "Enabling workaround for long-standing OpenMPI breakage (ompi/pull/1495)")
                        # tl;dr leave_pinned is just plain broken
                        # throughout most of the 1.7, 1.8. 1.9, and 1.10
                        # series of OpenMPI. The work to fix this is at
                        # https://github.com/open-mpi/ompi/pull/1495 ;
                        # see the attached commit logs (namely, commits
                        # 57035744 and 4b7cd1c0 in ompi-release carry the
                        # fix. Those are open-mpi/ompi@7aa03d66 and
                        # open-mpi/ompi@11e2d788 in the ompi repository).
                        set(MPIEXEC_EXTRA_STANZAS "--mca mpi_leave_pinned 0")
                        # This is solely to fix
                        # https://github.com/open-mpi/ompi/issues/299 but
                        # the problem is broader than that.
                        # set(MPI_C_COMPILER_CMDLINE_INSERTIONS "--openmpi:linkall")
                        # set(MPI_CXX_COMPILER_CMDLINE_INSERTIONS "--openmpi:linkall")
                    endif()
                else()
                    message(STATUS "MPI C Compiler front-end not recognized, proceeding anyway")
                endif()
            endif()
        endif()
        if (MPI_C_COMPILER_CMDLINE_INSERTIONS)
            message(STATUS "Adding ${MPI_C_COMPILER_CMDLINE_INSERTIONS} to mpicc command line")
        endif()
        if (MPI_CXX_COMPILER_CMDLINE_INSERTIONS)
            message(STATUS "Adding ${MPI_CXX_COMPILER_CMDLINE_INSERTIONS} to mpicxx command line")
        endif()
        # Now check for the MPI API version.
        macro(my_try_compile_mpicc SOURCE VAR)
            string(RANDOM LENGTH 8 ALPHABET "0123456789abcdef" uuid)
            string(TOLOWER "${VAR}" lcvar)
            set(checkfile
                "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/check_${lcvar}_${UUID}.c")

            FILE(WRITE "${checkfile}" "${SOURCE}\n")
            MESSAGE(STATUS "Performing Test ${VAR}")
            execute_process(
                COMMAND ${MPI_C_COMPILER} -c ${checkfile}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                RESULT_VARIABLE test_return_code
                OUTPUT_VARIABLE test_stdout
                ERROR_VARIABLE test_stderr
                )
            FILE(APPEND
                ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Summary for ${VAR} check:\n"
                "Source file was ${checkfile}\n"
                "Return code: ${test_return_code}\n"
                "stdout: ${test_stdout}\n"
                "stderr: ${test_stderr}\n"
                "\n"
                )
            if (test_return_code)
                SET(${VAR} 0 CACHE INTERNAL "Test ${VAR}")
                MESSAGE(STATUS "Performing Test ${VAR} -- Failed")
            else()
                SET(${VAR} 1 CACHE INTERNAL "Test ${VAR}")
                MESSAGE(STATUS "Performing Test ${VAR} -- Success")
            endif()
        endmacro()

        my_try_compile_mpicc("
        #include <mpi.h>
        #if !((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1))
        #error \"MPI version 2.1 not supported\"
        #endif
        int main(int argc, char** argv)
        {
        MPI_Init(&argc,&argv);
        MPI_Finalize();
        }
        " HAVE_MPI2_API)

        my_try_compile_mpicc("
        #include <mpi.h>
        #if !((MPI_VERSION > 3) || (MPI_VERSION == 3 && MPI_SUBVERSION >= 0))
        #error \"MPI version 3.0 not supported\"
        #endif
        int main(int argc, char** argv)
        {
        MPI_Init(&argc,&argv);
        MPI_Finalize();
        }
        " HAVE_MPI3_API)

    else()
        message(STATUS "Value found for MPI C compiler ${MPI_C_COMPILER}")
        message(STATUS "Value found for MPI C++ compiler ${MPI_CXX_COMPILER}")
        message(STATUS "Value found for MPI driver ${MPIEXEC}")
        message(FATAL_ERROR "Cannot find all of mpicc/mpic++/mpiexec with MPI=$ENV{MPI}")
    endif()
endif()

