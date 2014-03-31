# - Look for an acceptable Python interpreter and test that importing the 
# module sqlite3 works

# Currently the cadofactor scripts use "python3" in the #! line, so we require
# a binary of that name to exist
find_program(PYTHON_EXECUTABLE NAMES python3)

if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Python interpreter not found")
endif()

# Test that the version is something sane
# First get the version string from "python --version"
execute_process(COMMAND "${PYTHON_EXECUTABLE}" --version OUTPUT_VARIABLE PYTHON_OUT ERROR_VARIABLE PYTHON_ERR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
if ("${PYTHON_OUT}" MATCHES "Python")
  set(_VERSION "${PYTHON_OUT}")
else()
  set(_VERSION "${PYTHON_ERR}")
endif()

string(REPLACE "Python " "" PYTHON_VERSION_STRING "${_VERSION}")
string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" PYTHON_VERSION_MAJOR "${PYTHON_VERSION_STRING}")
string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+.*" "\\1" PYTHON_VERSION_MINOR "${PYTHON_VERSION_STRING}")
string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" PYTHON_VERSION_PATCH "${PYTHON_VERSION_STRING}")

# List acceptable Python versions. Although 3.4 is not released at this time,
# presumably will be ok, too
set(_Python_ACCEPTED 3.4 3.3 3.2)

# Check that the interpreter is one of the accepted versions
foreach(_TEST_VERSION ${_Python_ACCEPTED})
  if (${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EQUAL _TEST_VERSION)
    set (_PYTHON_OK 1)
  endif()
endforeach()

if (NOT _PYTHON_OK)
  message(FATAL_ERROR "Did not find a Python interpreter of one of the versions: ${_Python_ACCEPTED}. Please see README.Python")
else()
  message(STATUS "Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH}: ${PYTHON_EXECUTABLE}")
endif()

# Check that importing the sqlite3 module works. Some distros don't have the
# sqlite3 librarby installed, and some omit the Python sqlite3 module... :-(

file(WRITE "${CMAKE_BINARY_DIR}python_sqlite3_test.py" "import sqlite3\n")
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_BINARY_DIR}python_sqlite3_test.py" RESULT_VARIABLE _returncode OUTPUT_QUIET ERROR_QUIET)
file(REMOVE "${CMAKE_BINARY_DIR}python_sqlite3_test.py")
if (NOT _returncode EQUAL 0)
  message(FATAL_ERROR "Importing the sqlite3 Python module failed. "
"This may be caused by the sqlite3 library package missing on your system. "
"This package is usually called \"sqlite3\"; please ensure via your system's "
"package manager that it is installed. "
"\n"
"It may also be caused by the Python library module called \"sqlite3\" not "
"being installed as part of your Python package. One known Linux distribution "
"with this problem is Gentoo. Please read Readme.sqlite for information "
"how to fix this problem.")
else()
  message(STATUS "Importing module sqlite3 in Python succeeded.")
endif()

mark_as_advanced(PYTHON_EXECUTABLE)
