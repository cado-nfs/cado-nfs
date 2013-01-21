# Test for some POSIX headers and functions that may not exist on non-Unix systems

INCLUDE (CheckIncludeFiles)

# The fact that CMake prints "Looking for include files HAVE_RESOURCE_H"
# rather than "Looking for include files sys/resource.h" is a bug:
# http://public.kitware.com/Bug/view.php?id=13484

CHECK_INCLUDE_FILES (sys/resource.h HAVE_RESOURCE_H)
CHECK_INCLUDE_FILES (sys/utsname.h HAVE_UTSNAME_H)
CHECK_INCLUDE_FILES (sys/mman.h HAVE_MMAN_H)
CHECK_INCLUDE_FILES (sys/statvfs.h HAVE_STATVFS_H)
CHECK_INCLUDE_FILES (sys/wait.h HAVE_WAIT_H)

INCLUDE (CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SIGHUP "signal.h" HAVE_SIGHUP) 

# Unset the CMake variable that search_for_function() interprets
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

include(${CADO_NFS_SOURCE_DIR}/config/search_for_function.cmake)
search_for_function(getc_unlocked HAVE_GETC_UNLOCKED)
search_for_function(nanosleep HAVE_NANOSLEEP)
search_for_function(usleep HAVE_USLEEP)
search_for_function(popen HAVE_POPEN)
search_for_function(pclose HAVE_PCLOSE)
search_for_function(getrusage HAVE_GETRUSAGE)
search_for_function(lrand48 HAVE_LRAND48)
search_for_function(strdup HAVE_STRDUP)
search_for_function(strndup HAVE_STRNDUP)
search_for_function(sigaction HAVE_SIGACTION)
search_for_function(waitpid HAVE_WAITPID)
