# - Look for an acceptable Python interpreter and test that importing the 
# module sqlite3 works

# Currently the cadofactor scripts use "python3" in the #! line, so we require
# a binary of that name to exist
find_program(PYTHON_EXECUTABLE NAMES python3)

if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Python interpreter not found")
endif()

# Test that the version is something sane
# First get the version string from "python3 --version"
execute_process(COMMAND "${PYTHON_EXECUTABLE}" --version OUTPUT_VARIABLE PYTHON_OUT ERROR_VARIABLE PYTHON_ERR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
if ("${PYTHON_OUT}" MATCHES "Python")
  set(_VERSION "${PYTHON_OUT}")
else()
  set(_VERSION "${PYTHON_ERR}")
endif()

message(STATUS "${PYTHON_EXECUTABLE} --version returned: ${_VERSION}")

string(REPLACE "Python " "" PYTHON_VERSION_STRING "${_VERSION}")

# Minumum acceptable Python version. Let's assume future versions are ok, too.
set(_Python_MINIMUM_ACCEPTED 3.6)

# Check that the interpreter is one of the accepted versions
if ("${PYTHON_VERSION_STRING}" VERSION_LESS "${_Python_MINIMUM_ACCEPTED}")
  message(FATAL_ERROR "Did not find a Python interpreter of version at least ${_Python_MINIMUM_ACCEPTED}. Please see README.Python")
endif()

# Check that importing the sqlite3 module works. Some distros don't have the
# sqlite3 library installed, and some omit the Python sqlite3 module... :-(

set(PYTHON_REQUIRED_MODS sqlite3 flask requests)
set(PYTHON_OPTIONAL_MODS mysql.connector)

foreach(mod ${PYTHON_REQUIRED_MODS} ${PYTHON_OPTIONAL_MODS})
    string(TOUPPER HAVE_PYTHON_${mod} var)
    string(MAKE_C_IDENTIFIER ${var} var)
    file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/python_${mod}_test.py" "import ${mod}\n")
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/python_${mod}_test.py" RESULT_VARIABLE _returncode OUTPUT_QUIET ERROR_QUIET)
    file(REMOVE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/python_${mod}_test.py")
    if (NOT _returncode EQUAL 0)
      set (${var} 0)
      if(mod IN_LIST PYTHON_REQUIRED_MODS)
      message(FATAL_ERROR "Importing the ${mod} Python module failed."
          " This may be caused by the ${mod} library package missing on"
          " your system. This package is usually called"
          " \"python3-${mod}\", \"python-${mod}\", or maybe just \"${mod}\";"
          " please ensure via your system's"
          " package manager that it is installed."
          " Alternatively, as a convenience means, we provide the"
          " script ./scripts/setup-venv.sh , which you can use to"
          " install the python requirements of cado-nfs in a venv."
          " Running this script is probably the quickest way to get"
          " you going.")
      else()
          message(STATUS "Importing optional module ${mod} in Python failed.")
      endif()
    else()
      if(mod IN_LIST PYTHON_REQUIRED_MODS)
          message(STATUS "Importing module ${mod} in Python succeeded.")
      else()
          message(STATUS "Importing optional module ${mod} in Python succeeded.")
      endif()
      set (${var} 1)
    endif()
endforeach()

mark_as_advanced(PYTHON_EXECUTABLE)
