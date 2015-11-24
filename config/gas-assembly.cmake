# gas-syntax assembly sources.
message(STATUS "Trying to compile gas-syntax assembly sources")
if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Intel")
    # -x assembler below is gcc syntax, understood by icc as well. Do not
    # attempt to do the same with other compilers.
    #
    # Note that renaming asm_sample.c to asm_sample.S does not work,
    # because cmake's try_compile insists on knowing which compiler to
    # use directly from the file name extension. So far we haven't
    # decided to go the enable_languages(ASM-ATT) route, but if we do,
    # then this would make the -x assembler below unnecessary ; at the
    # expense of doing recognition for gas on top of all the rest.
    # https://cmake.org/Wiki/CMake/Assembler
try_compile(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/asm_sample.c
    COMPILE_DEFINITIONS "-x assembler")
if(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
    message(STATUS "Trying to compile gas-syntax assembly sources -- Success")
else()
    message(STATUS "Trying to compile gas-syntax assembly sources -- Failed")
endif()
else()
    message(STATUS "Trying to compile gas-syntax assembly sources -- disabled (need GCC)")
endif()

