message(STATUS "Testing C++ access to FP environment (fenv)")
try_compile(HAVE_CXX_FENV
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/fenv.cpp)
if(HAVE_CXX_FENV)
    message(STATUS "Testing C++ access to FP environment (fenv) -- Success")
else()
    message(STATUS "Testing C++ access to FP environment (fenv) -- Failed")
endif()

