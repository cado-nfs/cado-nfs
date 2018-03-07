#ifndef MMAP_ALLOCATOR_H
#define MMAP_ALLOCATOR_H

#include "macros.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <memory>
#include <limits>
#include <string>
#include <vector>
#include <stdexcept>

/* This is inspired from https://github.com/johannesthoma/mmap_allocator
 * License is LGPL.
 * The version here has been trimmed down significantly, look for uptream
 * version for more features */

namespace mmap_allocator_namespace
{
    enum access_mode {
        DEFAULT_STL_ALLOCATOR, /* Default STL allocator (malloc based). Reason is to have containers that do both and are compatible */
        READ_ONLY,  /* Readonly modus. Segfaults when vector content is written to */
        READ_WRITE_PRIVATE, /* Read/write access, writes are not propagated to disk */
        READ_WRITE_SHARED  /* Read/write access, writes are propagated to disk (file is modified) */
    };

    struct mmap_allocator_exception: public std::runtime_error {
        mmap_allocator_exception(const char *msg_param):
            std::runtime_error(msg_param) { }

        mmap_allocator_exception(std::string msg_param):
            std::runtime_error(msg_param) { }

        virtual ~mmap_allocator_exception(void) { }
    };

    /* Only the mmapped_file is something we can use to grab mapping
     * segments. Once the mmapped_file object is destroyed, the segments
     * survive because of shared_ptr's.
     */
    class mmapped_file {
        public:
            typedef off_t offset_type;
            typedef size_t size_type;
            private:
            class mapping {
                int fd = -1;
                void * area = NULL;
                offset_type offset_mapped;     /* page-aligned */
                size_type length_mapped;       /* page-aligned */
                public:
                mapping(const char * fname, enum access_mode access_mode, offset_type offset, size_type length);
                ~mapping();
                void * get(offset_type, size_type);
                void put(void *, offset_type, size_type);
            };
            std::shared_ptr<mapping> m;
            public:
            class segment {
                public:
                std::shared_ptr<mapping> file;
                offset_type offset;
                size_type length;
                operator bool() const { return file != nullptr; }
                void * get(size_type s) {
                    if (s == 0) return NULL;
                    ASSERT_ALWAYS(s == length);
                    return file->get(offset, length);
                }
                void put(void * p, size_type s) {
                    if (s == 0) return;
                    ASSERT_ALWAYS(s == length);
                    file->put(p, offset, length);
                }
            };
        public:
            segment get_segment(offset_type offset, size_type length) {
                return segment { m, offset, length };
            }
            mmapped_file(const char * fname, access_mode mode = READ_ONLY, offset_type offset = 0, size_type length = std::numeric_limits<size_type>::max()) : m(std::make_shared<mapping>(fname, mode, offset, length)) {}
    };

    template <typename T, typename A> class mmappable_vector;

    template <typename T> class mmap_allocator: public std::allocator<T>
    {
        public:
            using typename std::allocator<T>::size_type;
            using typename std::allocator<T>::pointer;
            using typename std::allocator<T>::const_pointer;
            typedef mmapped_file::offset_type offset_type;
        private:
            mmapped_file::segment s;
            friend class mmappable_vector<T, mmap_allocator<T> >;
        public:
            template<typename _Tp1>
            struct rebind { typedef mmap_allocator<_Tp1> other; };

            pointer allocate(size_type n, const void *hint=0)
            {
                if (!s)
                    return std::allocator<T>::allocate(n, hint);

                if (n == 0) return NULL;

                return (pointer) s.get(n*sizeof(T));
            }

            void deallocate(pointer p, size_type n)
            {
                if (!s) {
                    std::allocator<T>::deallocate(p, n);
                    return;
                }

                if (n == 0) return;

                s.put(p, n*sizeof(T));
            }

            mmap_allocator() = default;
            mmap_allocator(const std::allocator<T> &a): std::allocator<T>(a) {}
            mmap_allocator(const mmap_allocator &a) = default;
            ~mmap_allocator() = default;

            /* These are the only ctors that enable a behaviour that
             * differs from the stl container */
            mmap_allocator(mmapped_file::segment const & s) : s(s) {}
            mmap_allocator(mmapped_file m, offset_type offset, size_type length) : s(m.get_segment(offset, length * sizeof(T))) {}
            mmap_allocator(const char * f, access_mode mode, offset_type offset, size_type length) : mmap_allocator(mmapped_file(f, mode, offset, length*sizeof(T)), offset, length) {}

    };
}

#endif /* MMAP_ALLOCATOR_H */
