if (NOT DEFINED ENV{NO_PYTHON_CHECK})
add_subdirectory (cadofactor)
endif()

if(HAVE_GDB_PYTHON_PRETTY_PRINT AND NOT VALGRIND)
add_subdirectory (gdb)
endif()
