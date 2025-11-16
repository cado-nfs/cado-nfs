#ifndef CADO_UTILS_MMAP_ALLOCATOR_HPP
#define CADO_UTILS_MMAP_ALLOCATOR_HPP

#include <cstdio>

#include <memory>
#include <limits>
#include <string>
#include <stdexcept>
#include <utility>

#include <sys/types.h>  // off_t // IWYU pragma: keep

#include "macros.h"

/* This is inspired from https://github.com/johannesthoma/mmap_allocator
 * License is LGPL.
 * The version here has been trimmed down significantly, look for uptream
 * version for more features */

namespace mmap_allocator_details
{
    // NOLINTNEXTLINE(misc-include-cleaner)
    using offset_type = off_t;
    using size_type = size_t;

    enum access_mode {
        DEFAULT_STL_ALLOCATOR, /* Default STL allocator (malloc based). Reason is to have containers that do both and are compatible */
        READ_ONLY,  /* Readonly modus. Segfaults when vector content is written to */
        READ_WRITE_PRIVATE, /* Read/write access, writes are not propagated to disk */
        READ_WRITE_SHARED  /* Read/write access, writes are propagated to disk (file is modified) */
    };

    struct mmap_allocator_exception: public std::runtime_error {
        explicit mmap_allocator_exception(const char *msg_param):
            std::runtime_error(msg_param) { }

        explicit mmap_allocator_exception(std::string const & msg_param):
            std::runtime_error(msg_param) { }
    };

    /* Only the mmapped_file is something we can use to grab mapping
     * segments. Once the mmapped_file object is destroyed, the segments
     * survive because of shared_ptr's.
     */
    class mmapped_file {
        public:
            private:
            class mapping {
                int fd = -1;
                void * area = nullptr;
                offset_type offset_mapped;     /* page-aligned */
                size_type length_mapped;       /* page-aligned */
                public:
                mapping(const char * filename, enum access_mode amode, offset_type offset, size_type length);
                mapping(std::string const & filename, enum access_mode amode, offset_type offset, size_type length)
                    : mapping(filename.c_str(), amode, offset, length)
                {}
                ~mapping();
                void * get(offset_type, size_type);
                void put(void *, offset_type, size_type);
                mapping(mapping const &) = delete;
                mapping(mapping &&) = delete;
                mapping& operator=(mapping const &) = delete;
                mapping& operator=(mapping &&) = delete;
            };
#if 0
            /* use this to trace shared_ptr games. */
            class shared_mapping : private std::shared_ptr<mapping> {
                public:
                template<typename... Args> shared_mapping(Args... args) :
                    std::shared_ptr<mapping>(args...) {
                        fprintf(stderr, "Creating a shared_ptr [%p]\n", (void*) get());
                    }
                shared_mapping(std::shared_ptr<mapping> const & s) : std::shared_ptr<mapping>(s) { 
                    fprintf(stderr, "Copying a shared_ptr [%p]\n", (void*) get());
                }
                shared_mapping(shared_mapping const & s) : std::shared_ptr<mapping>(s) {
                    fprintf(stderr, "Copying a shared_ptr [%p]\n", (void*) get());
                }
                ~shared_mapping() {
                    fprintf(stderr, "Destroying a shared_ptr [%p]\n", (void*) get());
                }
                mapping * operator->() { return std::shared_ptr<mapping>::operator->(); }
                mapping const * operator->() const { return std::shared_ptr<mapping>::operator->(); }
                operator bool() const { return (std::shared_ptr<mapping> const&)(*this) != nullptr; }
            };
#endif
            std::shared_ptr<mapping> m;
            public:
            class segment {
                public:
                std::shared_ptr<mapping> file;
                offset_type offset;
                size_type length;
                /* with shared_mapping above that defines an operator
                 * bool, we just say "return file" */
                explicit operator bool() const { return file != nullptr; }
                void * get(size_type s) {
                    if (s == 0) return nullptr;
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
                return { m, offset, length };
            }
            explicit mmapped_file(const char * fname, access_mode mode = READ_ONLY, offset_type offset = 0, size_type length = std::numeric_limits<size_type>::max()) : m(std::make_shared<mapping>(fname, mode, offset, length)) {}
            explicit mmapped_file(std::string const & fname, access_mode mode = READ_ONLY, offset_type offset = 0, size_type length = std::numeric_limits<size_type>::max()) : mmapped_file(fname.c_str(), mode, offset, length) {}
    };

    template <typename T> class mmap_allocator: public std::allocator<T>
    {
        private:
            mmapped_file::segment s;
        public:
            bool has_defined_mapping() const { return bool(s); }
            template<typename Tp1>
            struct rebind { using other = mmap_allocator<Tp1>; };

            T * allocate(size_t n)
            {
                if (!s) {
                    return std::allocator<T>::allocate(n);
                }

                if (n == 0) return nullptr;

                return (T *) s.get(n*sizeof(T));
            }

            void deallocate(T * p, size_t n)
            {
                if (!s) {
                    std::allocator<T>::deallocate(p, n);
                    return;
                }

                if (n == 0) return;

                s.put(p, n*sizeof(T));
            }

            explicit mmap_allocator(const std::allocator<T> &a): std::allocator<T>(a) {}
            mmap_allocator() = default;

            /* These are the only ctors that enable a behaviour that
             * differs from the stl container */
            explicit mmap_allocator(mmapped_file::segment s)
                : s(std::move(s))
            {}
            mmap_allocator(mmapped_file m, offset_type offset, size_t length) : s(m.get_segment(offset, length * sizeof(T))) {}
            mmap_allocator(const char * f, access_mode mode, offset_type offset, size_t length) : mmap_allocator(mmapped_file(f, mode, offset, length*sizeof(T)), offset, length) {}

    };
}

#endif /* CADO_UTILS_MMAP_ALLOCATOR_HPP */
