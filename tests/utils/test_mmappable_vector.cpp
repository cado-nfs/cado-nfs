#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>

#include "macros.h"
#include "mmap_allocator.hpp"
#include "mmappable_vector.hpp"
#include "portability.h" // strdup // IWYU pragma: keep

/* This is inspired from https://github.com/johannesthoma/mmap_allocator
 * License is LGPL.
 * The version here has been trimmed down significantly, look for uptream
 * version for more features */

using namespace std;
using namespace mmap_allocator_details;

struct testcase {
    std::string TESTFILE;
    std::string TESTFILE2;

    static void generate_test_file(int count, std::string const & fname)
    {
        FILE * f = fopen(fname.c_str(), "w+");
        DIE_ERRNO_DIAG(f == nullptr, "fopen(%s)", fname.c_str());
        for (int i = 0; i < count; i++) {
            fwrite(&i, 1, sizeof(i), f);
        }
        fclose(f);
    }

    static void test_test_file(int count, std::string const & fname, bool expect_zeros)
    {
        FILE * f = fopen(fname.c_str(), "r");
        DIE_ERRNO_DIAG(f == nullptr, "fopen(%s)", fname.c_str());
        for (int i = 0; i < count; i++) {
            int j;
            size_t const n = fread(&j, 1, sizeof(j), f);
            ASSERT_ALWAYS(n == sizeof(j));
            ASSERT_ALWAYS(j == (expect_zeros ? 0 : i));
        }
        fclose(f);
    }

    void test_mmap() const
    {
        generate_test_file(1024, TESTFILE);
        generate_test_file(1024 * 1024, TESTFILE2);

        {
            fprintf(stderr, "Testing R/O mapping\n");
            mmapped_file const M(TESTFILE, READ_ONLY);
            mmappable_vector<uint32_t> int_vec_default(
                mmap_allocator<uint32_t>(M, 0, 1024));
            int_vec_default.mmap(1024);
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (uint32_t i = 0; i < 1024; i++) {
                // this will segfault because the mapping is read-only.
                // int_vec_default[i] = i;
                ASSERT_ALWAYS(int_vec_default[i] == i); /* just to be sure */
            }
            test_test_file(1024, TESTFILE, false);
        }

        /* Now do the same test read-write */
        {
            fprintf(stderr, "Testing private R/W mapping\n");
            mmapped_file const M(TESTFILE, READ_WRITE_PRIVATE);
            mmappable_vector<uint32_t> int_vec_default(
                mmap_allocator<uint32_t>(M, 0, 1024));
            int_vec_default.mmap(1024);
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (int i = 0; i < 1024; i++) {
                int_vec_default[i] = 0; /* should not segfault */
            }
            /* because the mapping is private, we should still have the
             * normal stuff */
            test_test_file(1024, TESTFILE, false);
        }

        /* If we test that with a shared read-write, this will be
         * different */
        {
            fprintf(stderr, "Testing shared R/W mapping\n");
            mmapped_file const M(TESTFILE, READ_WRITE_SHARED);
            mmappable_vector<uint32_t> int_vec_default(
                mmap_allocator<uint32_t>(M, 0, 1024));
            int_vec_default.mmap(1024);
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (int i = 0; i < 1024; i++) {
                int_vec_default[i] = 0; /* should not segfault */
            }
            /* because the mapping is shared, we should see zeroes now.  */
            test_test_file(1024, TESTFILE, true);

            /* clean up our mess */
            generate_test_file(1024, TESTFILE);
        }

        /* Now how does it go if we map only part of a file */
        {
            fprintf(stderr, "Testing fragment mapping\n");
            mmapped_file const M(TESTFILE2, READ_WRITE_SHARED, 8000, 1040576);
            mmappable_vector<uint32_t> int_vec_default(
                mmap_allocator<uint32_t>(M, 8000, 1024));
            int_vec_default.mmap(1024);
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (int i = 0; i < 1024; i++) {
                ASSERT_ALWAYS(int_vec_default[i] == (uint32_t) (i + 2000)); /* just to be sure */
            }
        }

        /* explicitly zero-initialize elements, but on a private mapping,
         * so that we don't see the result in the file.
         */
        {
            fprintf(stderr, "Testing value-initialized + private mapping\n");
            mmapped_file const M(TESTFILE, READ_WRITE_PRIVATE);
            mmappable_vector<uint32_t> int_vec_default(
                1024, 0, mmap_allocator<uint32_t>(M, 0, 1024));
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (int i = 0; i < 1024; i++) {
                ASSERT_ALWAYS(int_vec_default[i] == 0); /* just to be sure */
            }
            test_test_file(1024, TESTFILE, false);
        }

        /* on a shared mapping, we're supposed to get the zeroes */
        {
            fprintf(stderr, "Testing value-initialized + shared mapping\n");
            mmapped_file const M(TESTFILE, READ_WRITE_SHARED);
            mmappable_vector<uint32_t> int_vec_default(
                1024, 0, mmap_allocator<uint32_t>(M, 0, 1024));
            ASSERT_ALWAYS(int_vec_default.size() == 1024);
            for (int i = 0; i < 1024; i++) {
                ASSERT_ALWAYS(int_vec_default[i] == 0); /* just to be sure */
            }
            test_test_file(1024, TESTFILE, true);

            /* clean up our mess */
            generate_test_file(1024, TESTFILE);
        }
    }

    void test_conversion() const
    {
        fprintf(stderr,
                "Testing conversion between STL vector and mmap vector.\n");
        generate_test_file(1024, TESTFILE);

        mmappable_vector<uint32_t> mmap_vector;
        mmap_vector.mmap_file(TESTFILE, READ_ONLY, 0, 1024);
        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(mmap_vector[i] == i);
        }

        vector<uint32_t> std_vector(mmap_vector.begin(), mmap_vector.end());
        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(std_vector[i] == i);
        }
        for (uint32_t i = 0; i < 1024; i++) {
            std_vector[i] *= 2;
        }
        mmappable_vector<uint32_t> mmap_vector2(std_vector.begin(),
                                           std_vector.end());
        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(mmap_vector2[i] == i * 2);
        }
    }

    void test_shortcut_interface() const
    {
        fprintf(stderr, "Testing shortcut interface\n");

        generate_test_file(1024, TESTFILE);

        mmappable_vector<uint32_t> vec;
        vec.mmap_file(TESTFILE, READ_ONLY, 0, 1024);

        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(vec[i] == i);
        }
        try {
            /* This is expected to fail, because the vector is already mapped */
            vec.mmap_file(TESTFILE, READ_ONLY, 0, 1024);
            ASSERT_ALWAYS(0);
        } catch (mmap_allocator_exception const & e) {
            fprintf(stderr, "Exception message (expected): %s\n", e.what());
        }
        vec.munmap_file();

        generate_test_file(2048, TESTFILE);
        vec.mmap_file(TESTFILE, READ_ONLY, 4096, 1024);
        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(vec[i] == i + 1024);
        }
    }

    void test_cache_bug() const
    {
        mmappable_vector<uint32_t> vec;

        fprintf(stderr, "Testing if wrong offset bug in pool is fixed.\n");
        generate_test_file(2048, TESTFILE);
        vec.mmap_file(TESTFILE, READ_ONLY, 4096, 1024);

        for (uint32_t i = 0; i < 1024; i++) {
            ASSERT_ALWAYS(vec[i] == i + 1024);
        }
    }

    static constexpr size_t FILESIZE = 1024 * 1024 * 16;

    void read_large_file(enum access_mode mode) const
    {
        // struct timeval t, t2;
        mmappable_vector<uint32_t> vec;

        // gettimeofday(&t, NULL);

        vec.mmap_file(TESTFILE, mode, 0, FILESIZE);
        for (size_t i = 0; i < FILESIZE; i++) {
            ASSERT_ALWAYS(vec[i] == (size_t) i);
        }
        // gettimeofday(&t2, NULL);
        // fprintf(stderr, "Mode: %d Time: %lu.%06lu\n", mode, (t2.tv_sec -
        // t.tv_sec)-(t2.tv_usec < t.tv_usec), (t2.tv_usec < t.tv_usec)*1000000
        // + (t2.tv_usec - t.tv_usec));
    }

    void test_large_file() const
    {
        fprintf(stderr, "Testing large file.\n");
        generate_test_file(FILESIZE, TESTFILE); /* 1G */

        read_large_file(READ_ONLY);
        read_large_file(READ_WRITE_PRIVATE);
        read_large_file(READ_WRITE_SHARED);
    }

    void test_multiple_open() const
    {
        generate_test_file(1024, TESTFILE);
        generate_test_file(1024, TESTFILE2);

        /* first we create a file mapping for each vector, which causes many
         * different mmap() calls to be issued */
        {
            fprintf(stderr,
                    "Testing multiple open (you need to strace this).\n");
            mmappable_vector<uint32_t> vec1, vec2, vec3, vec4;
            vec1.mmap_file(TESTFILE, READ_ONLY, 0, 1024);
            vec2.mmap_file(TESTFILE, READ_ONLY, 0, 1024);
            vec3.mmap_file(TESTFILE2, READ_ONLY, 0, 1024);
            vec4.mmap_file(TESTFILE2, READ_ONLY, 0, 1024);
        }

        /* It is better to specifically create a mapping for the file, and
         * use it. That way, we have only one mmap() per file */
        {
            mmapped_file const M(TESTFILE, READ_ONLY);
            mmapped_file const M2(TESTFILE2, READ_ONLY);
            fprintf(stderr,
                    "Testing multiple open (you need to strace this).\n");
            mmappable_vector<uint32_t> vec1(mmap_allocator<uint32_t>(M, 0, 1024));
            mmappable_vector<uint32_t> vec2(mmap_allocator<uint32_t>(M, 0, 1024));
            mmappable_vector<uint32_t> vec3(mmap_allocator<uint32_t>(M, 0, 1024));
            mmappable_vector<uint32_t> vec4(mmap_allocator<uint32_t>(M, 0, 1024));
            vec1.mmap(1024);
            vec2.mmap(1024);
            vec3.mmap(1024);
            vec4.mmap(1024);
        }
    }

    void test_allocate_0_bytes() const /* shouldn't segfault */
    {
        fprintf(stderr, "Testing vectors of mmappable_vectors.\n");

        vector<mmappable_vector<uint32_t>> vecs;
        vecs.resize(2);
        for (int i = 0; i < 2; i++) {
            vecs[i].mmap_file(TESTFILE, READ_ONLY, 0, 1024);
            for (uint32_t j = 0; j < 1024; j++) {
                ASSERT_ALWAYS(vecs[i][j] == j);
            }
        }
    }

    void test_all() const
    {
        test_mmap();
        test_conversion();
        test_cache_bug();

        test_shortcut_interface();
        test_large_file();
        test_multiple_open();
        test_allocate_0_bytes();
    }

    explicit testcase(std::string const & tmpdir)
        : TESTFILE(tmpdir + "/testfile")
        , TESTFILE2(tmpdir + "/testfile2")
    {
    }
};

// coverity[root_function]
int main(int argc MAYBE_UNUSED, char const * argv[] MAYBE_UNUSED)
{
    char const * tmpdir = "/tmp";

    char const * env_wdir = getenv("wdir");
    if (env_wdir)
        tmpdir = env_wdir;

    testcase(tmpdir).test_all();
}
