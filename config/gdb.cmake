
find_program(GDB_EXECUTABLE NAMES gdb)

execute_process(COMMAND "${GDB_EXECUTABLE}" "--version"
    OUTPUT_VARIABLE GDB_OUT)

# message(STATUS "gdb version info: ${GDB_OUT}")
string(REGEX MATCH "^[^\n\r]+" GDB_OUT "${GDB_OUT}")
if(GDB_OUT)
    # message(STATUS "gdb version info: ${GDB_OUT}")
    # ubuntu has versions such as 15.0.50.20240403-git
    string(REGEX MATCH "([0-9a-z\.-]+)$" GDB_VERSION "${GDB_OUT}")
    # message(STATUS "gdb version info: ${GDB_VERSION}")
    if(GDB_VERSION)
        # message(STATUS "gdb version is ${CMAKE_MATCH_0}")
        # set(GDB_VERSION "${CMAKE_MATCH_0}")
        # to be honest, I have no idea what is the minimum version.
        if(GDB_VERSION VERSION_GREATER_EQUAL 15.0)
            set(HAVE_GDB_PYTHON_PRETTY_PRINT 1)
            message(STATUS "GNU gdb ${GDB_VERSION} detected, testing custom type printers")
        else()
            message(STATUS "GNU gdb ${GDB_VERSION} is probably too old for custom type printers")
        endif()
    endif()
endif()

