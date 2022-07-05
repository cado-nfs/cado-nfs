
message(STATUS "Testing whether we know about the std::vector internals")
try_compile(HAVE_KNOWN_CXX_STD_VECTOR_INTERNALS
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/cxx_std_vector_internals.cpp

    CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${PROJECT_SOURCE_DIR}/utils
    )
if(HAVE_KNOWN_CXX_STD_VECTOR_INTERNALS)
    message(STATUS "Testing whether we know about the std::vector internals -- yes")
else()
    message(STATUS "Testing whether we know about the std::vector internals -- no")
endif()

