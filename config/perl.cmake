# - Look for an acceptable perl interpreter and test that importing the 
# module Digest::MD5 works. It should: it's been part of the standard
# modules for ages. But some distros prefer to take it out.

find_program(PERL_EXECUTABLE NAMES perl)

if(NOT PERL_EXECUTABLE)
    message(FATAL_ERROR "Perl interpreter not found")
endif()

# Test that the version is something sane
# First get the version string from $^V
execute_process(COMMAND "${PERL_EXECUTABLE}" "-e" "print \$^V;"
    OUTPUT_VARIABLE PERL_OUT)
string(REGEX REPLACE "^v" "" PERL_VERSION "${PERL_OUT}")
message(STATUS "Perl version is ${PERL_VERSION}")

# Minumum acceptable perl version. Well, to be honest we have so minimal
# requirements that probably any perl5 will do.
set(PERL_MINIMUM_ACCEPTED 5.10)

# Check that the interpreter is one of the accepted versions
if (PERL_VERSION VERSION_LESS PERL_MINIMUM_ACCEPTED)
    message(FATAL_ERROR "Could not find a Perl interpreter of version at least ${PERL_MINIMUM_ACCEPTED}.")
endif()

# Check that importing the Digest::MD5 module works.
execute_process(COMMAND "${PERL_EXECUTABLE}" -MDigest::MD5 -e 1 RESULT_VARIABLE _returncode OUTPUT_QUIET ERROR_QUIET)
if (NOT _returncode EQUAL 0)
    message(WARNING
        "Importing the standard Perl module Digest::MD5 has failed.  This is probably because the Digest::MD5 module, despite being among the standard set, is not installed alongside with Perl on your system. Some tests will be skipped."
        )
    set(HAVE_PERL_DIGEST_MD5 0)
else()
    set(HAVE_PERL_DIGEST_MD5 1)
endif()

mark_as_advanced(PERL_EXECUTABLE)
