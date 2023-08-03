
# The first argument of this macro is going to be used for the default
# compiler flag. But in fact, all strings that are formed from this can
# be customized with the optional arguments:
#  DISPLAYNAME <blah>   -> What gets printed is "Testing whether <blah> code can be used"
#  COMPILER_FLAG <blah> -> Optionally test with this flag (and not -m${archfeature})
#  HAVE <blah>          -> The resulting macro is HAVE_<blah>
#  TESTFILE <blah>      -> Use <blah> instead of config/${archfeature (sanitized)}.c

macro(check_cpu_feature archfeature)
    set(prerequisites 1)
    set(blockers)
    set(archfeaturename ${archfeature})
    set(compilerflag "-m${archfeature}")
    string(MAKE_C_IDENTIFIER "${archfeature}" archfeaturevar)
    string(REPLACE "_" "" archfeaturevar "${archfeaturevar}")
    string(TOUPPER "${archfeaturevar}" ARCHFEATUREVAR)
    set(testfile "${CMAKE_CURRENT_SOURCE_DIR}/config/${archfeaturevar}.c")
    set(current prerequisites)
    foreach(x ${ARGN})
        if (x STREQUAL "PREREQUISITES")
            SET(current prerequisites)
        elseif (x STREQUAL "BLOCKERS")
            SET(current blockers)
        elseif (x STREQUAL "DISPLAYNAME")
            SET(current archfeaturename)
        elseif (x STREQUAL "COMPILER_FLAG")
            SET(current compilerflag)
        elseif (x STREQUAL "HAVE")
            SET(current ARCHFEATUREVAR)
        elseif (x STREQUAL "TESTFILE")
            SET(current testfile)
        else()
            set(${current} ${x})
        endif()
    endforeach()

    message(STATUS "Testing whether ${archfeaturename} code can be used")

    if (${blockers})
        message(STATUS "Testing whether ${archfeaturename} code can be used -- skipped because ${blockers} is set")
    elseif (${prerequisites})
        # Try with the system as it is configured, with the default
        # compiler behaviour. This may include tweaks by the used, we
        # don't care. If it compiles, we're happy. And we do __NOT__ run
        # the code in that case, because it's the user's business to pass
        # a compiler that creates runnable code. Or maybe not runnable
        # code, there could be use cases for that (compiling for a
        # different microarchitecture for instance).


        # We have a problem with icc, which by default seems to compile
        # to the latest instruction set, without consideration of the
        # host machine. I know that this sounds totally bogus, but it's
        # apparently the way it is. Flags such as -xhost don't have any
        # impact. Apparently, any call to the intel intrinsics emits the
        # corresponding assembly instruction, regardless of the host
        # machine. So for ICC, short of knowing how to deal with this, we
        # have to effectively disable cross-microarchitecture
        # compilation and use try_run instead of try_compile

        if(CMAKE_C_COMPILER_ID MATCHES "Intel")
            try_run(${archfeaturevar}_runs ${archfeaturevar}_compiles
                ${PROJECT_BINARY_DIR}/config
                ${testfile})
            set(happy ${${archfeaturevar}_runs})
        else()
            try_compile(${archfeaturevar}_compiles
                ${PROJECT_BINARY_DIR}/config
                ${testfile})
            set(happy ${${archfeaturevar}_compiles})
        endif()
        if(${happy})
            message(STATUS "Testing whether ${archfeaturename} code can be used -- Yes")
            set (HAVE_${ARCHFEATUREVAR} 1)
        elseif(CMAKE_C_FLAGS MATCHES "-march")
            message(STATUS "Testing whether ${archfeaturename} code can be used -- No (not testing ${compilerflag} because -march is already present)")
        else()
            # We're going to _try_ to add the feature with an explicit
            # compiler flag, but noting that it makes sense to do this
            # _only_ if we can check that the resulting binary can run.
            try_run(${archfeaturevar}_runs ${archfeaturevar}_compiles
              ${PROJECT_BINARY_DIR}/config
              ${testfile}
              COMPILE_DEFINITIONS ${compilerflag})
            if(${archfeaturevar}_compiles)
                if (${archfeaturevar}_runs MATCHES FAILED_TO_RUN)
                    message(STATUS "Testing whether ${archfeaturename} code can be used -- No")
                    set (HAVE_${ARCHFEATUREVAR} 0)
                else()
                    message(STATUS "Testing whether ${archfeaturename} code can be used -- Yes, with ${compilerflag}")
                    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${compilerflag}")
                    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${compilerflag}")
                    set (HAVE_${ARCHFEATUREVAR} 1)
                endif()
            else()
                message(STATUS "Testing whether ${archfeaturename} code can be used -- No")
                set (HAVE_${ARCHFEATUREVAR} 0)
            endif()
        endif()
    else()
        message(STATUS "Testing whether ${archfeaturename} code can be used -- skipped because prerequisites are not met")
    endif()
endmacro()
