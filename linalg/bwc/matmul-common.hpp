#ifndef CADO_MATMUL_COMMON_HPP
#define CADO_MATMUL_COMMON_HPP

#include <cstdio>
#include <cstdint>

#include <memory>

#include <sys/stat.h>

#include "utils_cxx.hpp"        // for unique_ptr<FILE, delete_FILE>
#include "macros.h"  // for DIE_ERRNO_DIAG, FATAL_ERROR_CHECK

struct matmul_public;

std::unique_ptr<FILE, delete_FILE> matmul_common_reload_cache_fopen(size_t, matmul_public & mm, uint32_t magic);
std::unique_ptr<FILE, delete_FILE> matmul_common_save_cache_fopen(size_t, matmul_public const & mm, uint32_t magic);

extern const char * const rowcol[2];  // [0] = "row" [1] = "col"

/* All matmul implementation are peppered with such markers */
#ifndef ASM_COMMENT
#ifdef  __GNUC__
#define ASM_COMMENT(x)  __asm__("#\t" x "\n")
#else
#define ASM_COMMENT(x)  /**/
#endif
#endif

/* Use this instead fo just vector::resize before using
 * MATMUL_COMMON_READ_MANYxxx ; this at least makes sure that we're not
 * allocating more than the file size!
 *
 * This must be called only when T is a type with rigid allocation size
 * (immediate POD types, or maybe pairs/tuples)
 */
template<typename T>
void resize_and_check_meaningful(std::vector<T> & a, size_t n, FILE * f)
{
    struct stat sbuf[1];
    int const rc = fstat(fileno(f), sbuf);
    ASSERT_ALWAYS(rc == 0);
    long here = ftell(f);
    ASSERT_ALWAYS(here >= 0);
    size_t there = (size_t) here + n * sizeof(T);
    ASSERT_ALWAYS(there <= (size_t) sbuf->st_size);
    a.resize(n);
}

/* I/O with cache files is made easier with these macros ; rather than
 * having to check for errors over and over again... */
#define MATMUL_COMMON_READ_ONE64(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint64_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    (final_v__) = storage_v__;                                          \
} while (0)
#define MATMUL_COMMON_READ_ONE32(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint32_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    (final_v__) = storage_v__;                                          \
} while (0)
#define MATMUL_COMMON_READ_ONE8(final_v__, file__)  do {                \
    size_t rc;                                                          \
    uint8_t storage_v__;                                                \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    (final_v__) = storage_v__;                                          \
} while (0)
#define MATMUL_COMMON_READ_ONE16(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint16_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    (final_v__) = storage_v__;                                          \
} while (0)
#define MATMUL_COMMON_READ_MANY32(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fread(ptr__, sizeof(uint32_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short read from cached matrix file");	\
} while (0)
#define MATMUL_COMMON_READ_MANY16(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fread(ptr__, sizeof(uint16_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short read from cached matrix file");	\
} while (0)
#define MATMUL_COMMON_READ_MANY8(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fread(ptr__, sizeof(uint8_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short read from cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_ONE64(final_v__, file__)  do {              \
    size_t rc;								\
    uint64_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    DIE_ERRNO_DIAG(rc < 1, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_ONE32(final_v__, file__)  do {              \
    size_t rc;								\
    uint32_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    DIE_ERRNO_DIAG(rc < 1, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_ONE8(final_v__, file__)  do {              \
    size_t rc;								\
    uint8_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    DIE_ERRNO_DIAG(rc < 1, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_ONE16(final_v__, file__)  do {              \
    size_t rc;								\
    uint16_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    DIE_ERRNO_DIAG(rc < 1, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_MANY32(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fwrite(ptr__, sizeof(uint32_t), n__, f__);			\
    DIE_ERRNO_DIAG(rc < n__, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_MANY16(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fwrite(ptr__, sizeof(uint16_t), n__, f__);			\
    DIE_ERRNO_DIAG(rc < n__, "write(%s)", "cached matrix file");	        \
} while (0)
#define MATMUL_COMMON_WRITE_MANY8(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fwrite(ptr__, sizeof(uint8_t), n__, f__);			\
    DIE_ERRNO_DIAG(rc < n__, "write(%s)", "cached matrix file");	        \
} while (0)


#endif	/* MATMUL_COMMON_HPP_ */
