build_tree="/tmp"
PREFIX=/usr/local
# It's a bit tricky. What should we use as a default? And in fact, why
# would we have to care? Could it make sense to have images with multiple
# arch settings available?
CFLAGS="-O3 -DNDEBUG -march=haswell"
CXXFLAGS="$CFLAGS"
BWC_GF2_ARITHMETIC_BACKENDS="m128;u64k1"
BWC_GFP_ARITHMETIC_BACKENDS=""
BWC_GF2_MATMUL_BACKENDS="bucket"
BWC_GFP_MATMUL_BACKENDS="zone"
ENABLE_SHARED=1
