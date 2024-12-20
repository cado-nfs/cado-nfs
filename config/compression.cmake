# These compression tools are defined in gzip.cpp, and therefore can be
# handled transparently by cado-nfs. Of course we're only willing to test
# them if the corresponding tools are available
find_program(HAVE_GZIP NAMES gzip)
find_program(HAVE_BZIP2 NAMES bzip2)
find_program(HAVE_ZSTD NAMES zstd)
find_program(HAVE_XZ NAMES xz)
find_program(HAVE_LZMA NAMES lzma)

if(NOT HAVE_GZIP)
    message(FATAL_ERROR "While not strictly necessary, the gzip tool is"
    " implicitly used in many of the cado-nfs scripts, and the test suite"
    " will fail without it. Please install it.")
endif()
