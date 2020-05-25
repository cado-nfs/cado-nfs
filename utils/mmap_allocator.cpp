#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <string>       // std::string
#include <exception>    // std::terminate

#include <sys/mman.h>  // for PROT_READ, mmap, munmap, MAP_SHARED, PROT_WRITE
#include <sys/stat.h>  // for fstat, stat
#include <fcntl.h>      // O_RDONLY etc // IWYU pragma: keep
#include <unistd.h>    // for sysconf, _SC_PAGE_SIZE, close

#include "mmap_allocator.hpp"


/* This is inspired from https://github.com/johannesthoma/mmap_allocator
 * License is LGPL.
 * The version here has been trimmed down significantly, look for uptream
 * version for more features */

#define ALIGN_TO_PAGE(x) ((x) & ~(sysconf(_SC_PAGE_SIZE) - 1))
#define UPPER_ALIGN_TO_PAGE(x) ALIGN_TO_PAGE((x)+(sysconf(_SC_PAGE_SIZE)-1))
#define OFFSET_INTO_PAGE(x) ((x) & (sysconf(_SC_PAGE_SIZE) - 1))

namespace mmap_allocator_details {
    mmapped_file::mapping::mapping(const char * filename, access_mode amode, offset_type offset, size_type length) {
        if (!filename || *filename == '\0') {
            throw mmap_allocator_exception("mmapped_file not correctly initialized: filename is empty.");
        }
        int mode;
        int prot;
        int mmap_mode = 0;

        switch (amode) {
            case READ_ONLY:
                mode = O_RDONLY;
                prot = PROT_READ;
                mmap_mode |= MAP_SHARED;
                break;
            case READ_WRITE_SHARED:
                mode = O_RDWR;
                prot = PROT_READ | PROT_WRITE;
                mmap_mode |= MAP_SHARED;
                break;
            case READ_WRITE_PRIVATE:
                mode = O_RDONLY;
                prot = PROT_READ | PROT_WRITE;
                mmap_mode |= MAP_PRIVATE;
                break;
            default:
                throw mmap_allocator_exception("Internal error");
                break;
        }

        fd = open(filename, mode);
        if (fd < 0)
            throw mmap_allocator_exception(std::string("Error opening file ") + filename);

        if (length == std::numeric_limits<size_type>::max()) {
            /* well, we really want the file length, not more ! */
            struct stat sbuf[1];
            if (fstat(fd, sbuf) < 0)
                throw mmap_allocator_exception("stat() error");
            length = sbuf->st_size;
        }
        offset_mapped = ALIGN_TO_PAGE(offset);
        length_mapped = UPPER_ALIGN_TO_PAGE(length + offset - offset_mapped);
        area = mmap(NULL, length_mapped, prot, mmap_mode, fd, offset_mapped);
    }

    mmapped_file::mapping::~mapping() {
        if (munmap(area, length_mapped) < 0)
            std::terminate(); /* munmap() error in dtor, fatal */
        if (close(fd) < 0)
            std::terminate(); /* close() error in dtor, fatal */
    }

    void * mmapped_file::mapping::get(offset_type offset, size_type length) {
        if (offset >= offset_mapped && length + offset <= offset_mapped + length_mapped) {
            return ((char*)area)+offset-offset_mapped;
        } else {
            throw mmap_allocator_exception("Cannot get range outside mapping bounds");
        }
    }

    void mmapped_file::mapping::put(void *, offset_type, size_type)
    {
        /* in fact, we do nothing. We _could_ imagine keeping track
         * of things, but what for, really ?
         */
    }
}
