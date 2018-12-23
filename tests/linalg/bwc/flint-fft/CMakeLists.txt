
# The flint-fft library is only used if GF(p) lowlevel arithmetic is
# enabled, and we never worked on that in the 32-bit context. So we need
# to disable the flint-fft-only checks if BWC_GFP_ARITHMETIC_BACKENDS is
# empty.
if(BWC_GFP_ARITHMETIC_BACKENDS)
    include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/flint-fft)
    include_directories(${PROJECT_SOURCE_DIR}/utils)
    cado_define_test(test_transform_length.cpp LIBRARIES flint-fft utils)
    set_target_properties(test_transform_length PROPERTIES
        COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
        LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()
