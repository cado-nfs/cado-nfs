
# Define CURL=1 in local.sh to enable curl backend for bwc/u64_dispatch
# (if you do not know what this means, you do not need this
# functionality)
if($ENV{CURL})
    set(WANT_CURL 1)
else()
    string(REGEX MATCH "^/.*$" WANT_CURL "$ENV{CURL}")
    if(WANT_CURL STREQUAL "")
        set(WANT_CURL 0)
    else()
        set(WANT_CURL 1)
    endif()
endif()

if(${WANT_CURL})
# You can force a path to curl/curl.h using the environment variables CURL, or
# CURL_INCDIR and CURL_LIBDIR (note that curl/curl.h is also searched)
if("$ENV{CURL}" MATCHES "^(1|YES|yes|ON|on|)$")
    set(HAS_CURL_OVERRIDE 0)
else()
    set(HAS_CURL_OVERRIDE 1)
endif()

if (HAS_CURL_OVERRIDE)
    message(STATUS "Adding $ENV{CURL} to the search path for cURL")
    set(CURL_INCDIR_HINTS "$ENV{CURL}/include" ${CURL_INCDIR_HINTS})
    set(CURL_INCDIR_HINTS "$ENV{CURL}"         ${CURL_INCDIR_HINTS})
    set(CURL_LIBDIR_HINTS "$ENV{CURL}/lib"     ${CURL_LIBDIR_HINTS})
    set(CURL_LIBDIR_HINTS "$ENV{CURL}/.libs"   ${CURL_LIBDIR_HINTS})
endif()

string(COMPARE NOTEQUAL "$ENV{CURL_INCDIR}" "" HAS_CURL_INCDIR_OVERRIDE)
if (HAS_CURL_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{CURL_INCDIR} to the search path for cURL")
    set(CURL_INCDIR_HINTS "$ENV{CURL_INCDIR}" ${CURL_INCDIR_HINTS})
endif()

string(COMPARE NOTEQUAL "$ENV{CURL_LIBDIR}" "" HAS_CURL_LIBDIR_OVERRIDE)
if (HAS_CURL_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{CURL_LIBDIR} to the search path for cURL")
    set(CURL_LIBDIR_HINTS "$ENV{CURL_LIBDIR}"     ${CURL_LIBDIR_HINTS})
endif()

find_path   (CURL_INCDIR curl/curl.h HINTS ${CURL_INCDIR_HINTS} DOC "cURL headers")
find_library(CURL_LIB    curl   HINTS ${CURL_LIBDIR_HINTS} DOC "cURL library")
# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "CURL_INCDIR=${CURL_INCDIR}")
message(STATUS "CURL_LIB=${CURL_LIB}")
string(COMPARE NOTEQUAL "${CURL_INCDIR}" CURL_INCDIR-NOTFOUND CURL_INCDIR_OK)
string(COMPARE NOTEQUAL "${CURL_LIB}" CURL_LIB-NOTFOUND CURL_LIBDIR_OK)

get_filename_component(CURL_LIBDIR ${CURL_LIB} PATH)

if(CURL_INCDIR_OK AND CURL_LIBDIR_OK)
include_directories(${CURL_INCDIR})
link_directories(${CURL_LIBDIR})
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-L${CURL_LIBDIR}")
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES ${CURL_INCDIR})
set(CMAKE_REQUIRED_LIBRARIES "-lcurl")
CHECK_C_SOURCE_COMPILES("
    #include <curl/curl.h>
    int main(void)
    {
      CURL *curl;
      curl = curl_easy_init();
      curl_easy_cleanup(curl);
      return 0;
    }
" HAVE_CURL)
endif()
endif()
