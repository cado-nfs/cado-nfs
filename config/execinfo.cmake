message(STATUS "Checking for BSD/Linux execinfo interface")
include(CheckCSourceCompiles)
set(OLD_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET 1)
CHECK_C_SOURCE_COMPILES("
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
static void signal_handling (int signum)
{
   int sz = 100, i;
   void *buffer [sz];
   char** text;

   sz = backtrace (buffer, sz);
   text = backtrace_symbols (buffer, sz);
   for (i = 0; i < sz; i++)
       puts(text[i]);
   signal (signum, SIG_DFL);
   raise (signum);
}
int main()
{
    signal (SIGABRT, signal_handling);
    signal (SIGSEGV, signal_handling);
}
" HAVE_EXECINFO)
set(CMAKE_REQUIRED_QUIET "${OLD_CMAKE_REQUIRED_QUIET}")

if(HAVE_EXECINFO)
    message(STATUS "Checking for BSD/Linux execinfo interface -- Yes")
else()
    message(STATUS "Checking for BSD/Linux execinfo interface -- No")
endif()
