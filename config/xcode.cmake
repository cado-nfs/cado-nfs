# I'm not happy with _cmake_find_compiler_path

if(CMAKE_HOST_APPLE)

    macro(_query_xcrun compiler_name result_var_keyword result_var)
	if(NOT "x${result_var_keyword}" STREQUAL "xRESULT_VAR")
	    message(FATAL_ERROR "Bad arguments to macro")
	endif()
	execute_process(COMMAND xcrun --find ${compiler_name}
	    OUTPUT_VARIABLE _xcrun_out OUTPUT_STRIP_TRAILING_WHITESPACE
	    ERROR_VARIABLE _xcrun_err)
	set("${result_var}" "${_xcrun_out}")
    endmacro()

    # We'll do a rewrite, _BUT_ not change the cache variable. After all, this
    # wrapping is only plumbing of ours.

    macro(rewrite_xcrun lang)
	if (CMAKE_${lang}_COMPILER MATCHES "^.*/usr/bin/(.+)$")
	    set(compiler_name "${CMAKE_MATCH_1}")
	    set(xcrun_result)
	    _query_xcrun("${compiler_name}" RESULT_VAR xcrun_result)
	    if (CMAKE_${lang}_COMPILER STREQUAL xcrun_result)
		    find_program(xcrun_path xcrun)
		    execute_process(COMMAND
			    ${CMAKE_COMMAND} -E create_symlink ${xcrun_path} ${PROJECT_BINARY_DIR}/${compiler_name})
		    message(STATUS "rewriting CMAKE_${lang}_COMPILER=${CMAKE_${lang}_COMPILER} to CMAKE_${lang}_COMPILER=\"${PROJECT_BINARY_DIR}/${compiler_name}\" (pointing to ${xcrun_path})")
		    # set_property(CACHE CMAKE_${lang}_COMPILER PROPERTY VALUE "${PROJECT_BINARY_DIR}/${compiler_name}")
		    set(CMAKE_${lang}_COMPILER "${PROJECT_BINARY_DIR}/${compiler_name}")
	    endif()
	endif()
    endmacro()

    rewrite_xcrun(C)
    rewrite_xcrun(CXX)
endif()

