# This checks for the presence of the fmt library. See https://fmt.dev/

# You can force a path to the fmt headers using the environment variables
# FMT, or FMT_INCDIR and FMT_LIBDIR

if (DEFINED ENV{FMT})
    message(STATUS "Adding $ENV{FMT} to the search path for fmt")
    set(FMT_INCDIR_HINTS ${FMT_INCDIR_HINTS} "$ENV{FMT}/include")
    set(FMT_INCDIR_HINTS ${FMT_INCDIR_HINTS} "$ENV{FMT}"        )
    set(FMT_LIBDIR_HINTS ${FMT_LIBDIR_HINTS} "$ENV{FMT}/lib"    )
endif()

if (DEFINED ENV{FMT_INCDIR})
    message(STATUS "Adding $ENV{FMT_INCDIR} to the search path for fmt")
    # prepend !
    set(FMT_INCDIR_HINTS "$ENV{FMT_INCDIR}" ${FMT_INCDIR_HINTS})
endif()

if (DEFINED ENV{FMT_LIBDIR})
    message(STATUS "Adding $ENV{FMT_LIBDIR} to the search path for fmt")
    # prepend !
    set(FMT_LIBDIR_HINTS "$ENV{FMT_LIBDIR}" ${FMT_LIBDIR_HINTS})
endif()

# Try in three passes, otherwise cmake gets in the way...
find_path   (FMT_INCDIR fmt/core.h HINTS ${FMT_INCDIR_HINTS} DOC "fmt headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
        )
if(NOT FMT_INCDIR)
find_path   (FMT_INCDIR fmt/core.h HINTS ${FMT_INCDIR_HINTS} DOC "fmt headers"
        NO_DEFAULT_PATH
        )
endif()
if(NOT FMT_INCDIR)
find_path   (FMT_INCDIR fmt/core.h HINTS ${FMT_INCDIR_HINTS} DOC "fmt headers")
endif()

find_library(FMT_LIB    fmt   HINTS ${FMT_LIBDIR_HINTS} DOC "fmt library"
    NO_DEFAULT_PATH
    )
if(NOT FMT_LIB)
find_library(FMT_LIB    fmt   HINTS ${FMT_LIBDIR_HINTS} DOC "fmt library")
endif()


# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "FMT_INCDIR=${FMT_INCDIR}")
message(STATUS "FMT_LIB=${FMT_LIB}")
string(COMPARE NOTEQUAL "${FMT_INCDIR}" FMT_INCDIR-NOTFOUND FMT_INCDIR_OK)
string(COMPARE NOTEQUAL "${FMT_LIB}" FMT_LIB-NOTFOUND FMT_LIBDIR_OK)

get_filename_component(FMT_LIBDIR ${FMT_LIB} PATH)

if(FMT_INCDIR_OK AND FMT_LIBDIR_OK)
include_directories(${FMT_INCDIR})
link_directories(${FMT_LIBDIR})
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-L${FMT_LIBDIR}")
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES ${FMT_INCDIR})
set(CMAKE_REQUIRED_LIBRARIES "-lfmt")
CHECK_CXX_SOURCE_COMPILES("
#include <iostream>
#include <fmt/core.h>
#ifdef FMT_VERSION
#if FMT_VERSION < 80000
#error \"too old, never checked\"
#endif
#else
#error \"MISSING FMT_VERSION\"
#endif
#include <fmt/format.h>

#if FMT_VERSION >= 90000
/* with fmt8, formatting a const reference to a type that defines a
 * conversion to a pointer seems to fail, even though we abide by the
 * recommended practice for the definition of the custom formatter.
 */
#include <sstream>
struct cxx_foo {
public:
    int x[1] { 42 };
    operator int const *() const { return x; }
};
template <> struct fmt::formatter<cxx_foo>: formatter<string_view> {
  // parse is inherited from formatter<string_view>.

  auto format(cxx_foo const & c, format_context& ctx) const
      -> format_context::iterator
      {
            std::ostringstream os;
            os << c.x[0];
            return formatter<string_view>::format( string_view(os.str()), ctx);
      }
};
#endif

int main(void)
{
#if FMT_VERSION >= 90000
    cxx_foo b;
    cxx_foo const & a(b);
    std::cout << fmt::format(FMT_STRING(\"{} {} {}\"), \"Catch\", 22, a) << std::endl;
#else
    std::cout << fmt::format(FMT_STRING(\"{} {}\"), \"Catch\", 22) << std::endl;
#endif
    return 0;
}

" HAVE_FMT)
if(HAVE_FMT)
    message(STATUS "Using the fmt library found on the system")
else()
    message(FATAL_ERROR "A version of the fmt library was found by cmake, but it apparently does not fit our criteria. This is a fatal error, as we cannot be sure that our embedded copy will work correctly in that case")
endif()
elseif(FMT_INCDIR_OK OR FMT_LIBDIR_OK)
    message(FATAL_ERROR "A partly installed version of the fmt library was found by cmake, but it lacks either the headers or the library. This is a fatal error, as we cannot be sure that our embedded copy will work correctly in that case")
else()
    message(STATUS "Using the embedded fmt library")
endif()
