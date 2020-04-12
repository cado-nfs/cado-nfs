# I'm not happy with _cmake_find_compiler_path

if(CMAKE_HOST_APPLE)

	message(STATUS "World")
    macro(_query_xcrun compiler_name result_var_keyword result_var)
	if(NOT "x${result_var_keyword}" STREQUAL "xRESULT_VAR")
	    message(FATAL_ERROR "Bad arguments to macro")
	endif()
	execute_process(COMMAND xcrun --find ${compiler_name}
	    OUTPUT_VARIABLE _xcrun_out OUTPUT_STRIP_TRAILING_WHITESPACE
	    ERROR_VARIABLE _xcrun_err)
	set("${result_var}" "${_xcrun_out}")
    endmacro()

    macro(rewrite_xcrun lang)
	if (CMAKE_${lang}_COMPILER MATCHES "^.*/usr/bin/(.+)$")
	    set(compiler_name "${CMAKE_MATCH_1}")
	    set(xcrun_result)
	    _query_xcrun("${compiler_name}" RESULT_VAR xcrun_result)
	    if (CMAKE_${lang}_COMPILER STREQUAL xcrun_result)
		message(STATUS "rewriting CMAKE_${lang}_COMPILER=${CMAKE_${lang}_COMPILER} to CMAKE_${lang}_COMPILER=\"xcrun ${compiler_name}\"")
		set_property(CACHE CMAKE_${lang}_COMPILER PROPERTY VALUE "xcrun ${compiler_name}")
	    endif()
	endif()
    endmacro()

    rewrite_xcrun(C)
    rewrite_xcrun(CXX)
endif()

