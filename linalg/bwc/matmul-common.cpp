#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstring>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <memory>

#include "fmt/base.h"

#include "matmul.hpp"   // for matmul_public_s
#include "matmul-common.hpp"
#include "verbose.h"
#include "utils_cxx.hpp"        // for unique_ptr<FILE, delete_FILE>

#define MM_COMMON_MAGIC 0xb0010003UL

const char * const rowcol[2] = { "row", "col", };

/* Factor out some stuff which turns out to appear fairly often */

std::unique_ptr<FILE, delete_FILE> matmul_common_reload_cache_fopen(size_t stride, matmul_public & mm, uint32_t magic)
{
    std::unique_ptr<FILE, delete_FILE> f;
    if (mm.cachefile_name.empty()) return f;
    f.reset(fopen(mm.cachefile_name.c_str(), "rb"));
    if (!f) return f;

    // mm.cachefile_name is a cache file for mm.locfile (which in
    // general never exists)

    uint32_t magic_check;
    MATMUL_COMMON_READ_ONE32(magic_check, f.get());
    
    if (magic_check != magic) {
        fmt::print(stderr, "Wrong magic in cached matrix file\n");
        f.reset();
        return f;
    }   
    
    MATMUL_COMMON_READ_ONE32(magic_check, f.get());
    if (magic_check != MM_COMMON_MAGIC) {
        fmt::print(stderr, "Wrong magic in cached matrix file\n");
        f.reset();
        return f;
    }   

    /* Four reserved bytes for alignment */
    MATMUL_COMMON_READ_ONE32(magic_check, f.get());

    uint32_t nbytes_check;
    MATMUL_COMMON_READ_ONE32(nbytes_check, f.get());
    /* It's not fatal. It only deserves a warning */
    if (nbytes_check != stride) {
        fmt::print(stderr, "Warning: cached matrix file fits data with different stride\n");
    }

    MATMUL_COMMON_READ_ONE32(mm.dim[0], f.get());
    MATMUL_COMMON_READ_ONE32(mm.dim[1], f.get());
    MATMUL_COMMON_READ_ONE64(mm.ncoeffs, f.get());

    return f;
}

std::unique_ptr<FILE, delete_FILE> matmul_common_save_cache_fopen(size_t stride, matmul_public const & mm, uint32_t magic)
{
    std::unique_ptr<FILE, delete_FILE> f;
    if (mm.cachefile_name.empty()) return f;
    f.reset(fopen(mm.cachefile_name.c_str(), "wb"));
    if (!f) {
        fmt::print(stderr, "Cannot open {} for writing: {}\n",
                mm.cachefile_name,
                strerror(errno));
        abort();
    }

    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
        fmt::print("Saving {} to cache file {}\n",
                mm.locfile,
                mm.cachefile_name);

    MATMUL_COMMON_WRITE_ONE32(magic, f.get());
    MATMUL_COMMON_WRITE_ONE32(MM_COMMON_MAGIC, f.get());
    MATMUL_COMMON_WRITE_ONE32(0, f.get());
    MATMUL_COMMON_WRITE_ONE32(stride, f.get());
    MATMUL_COMMON_WRITE_ONE32(mm.dim[0], f.get());
    MATMUL_COMMON_WRITE_ONE32(mm.dim[1], f.get());
    MATMUL_COMMON_WRITE_ONE64(mm.ncoeffs, f.get());

    return f;
}
