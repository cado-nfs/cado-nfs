
# The first argument of this macro is going to be used for the default
# compiler flag. But in fact, all strings that are formed from this can
# be customized with the optional arguments:
#  DISPLAYNAME <blah>   -> What gets printed is "Testing whether <blah> code can be used"
#  COMPILER_FLAG <blah> -> Optionally test with this flag (and not -m${archfeature})
#  HAVE <blah>          -> The resulting macro is HAVE_<blah>
#  TESTFILE <blah>      -> Use <blah> instead of config/${archfeature (sanitized)}.c

macro(check_cpu_feature archfeature)
    set(prerequisites 1)
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

    if (${prerequisites})
        try_run(${archfeaturevar}_runs ${archfeaturevar}_compiles
            ${PROJECT_BINARY_DIR}/config
            ${testfile}
            )
        if(${archfeaturevar}_compiles)
            if (${archfeaturevar}_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether ${archfeaturename} code can be used -- No")
                set (HAVE_${ARCHFEATUREVAR} 0)
            else()
                message(STATUS "Testing whether ${archfeaturename} code can be used -- Yes")
                set (HAVE_${ARCHFEATUREVAR} 1)
            endif()
        elseif(CMAKE_C_FLAGS MATCHES "-march")
            message(STATUS "Testing whether ${archfeaturename} code can be used -- No (not testing ${compilerflag} because -march is already present)")
        else()
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
        message(STATUS "Testing whether ${archfeaturename} code can be used -- skipped")
    endif()
endmacro()
