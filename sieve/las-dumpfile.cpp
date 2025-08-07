#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>

#include <gmp.h>

#include "las-dumpfile.hpp"
#include "special-q.hpp"
#include "macros.h"

dumpfile_t::~dumpfile_t() {
    if (f) fclose(f);
}

void dumpfile_t::close() {
    if (f) fclose(f);
}

void dumpfile_t::open(const char *filename_stem, special_q const & doing, int side)
{
    ASSERT_ALWAYS(!f);
    if (filename_stem != NULL) {
        char *filename;
        int const rc = gmp_asprintf(&filename, "%s.%d.sq%Zd.rho%Zd.side%d.dump",
            filename_stem,
            doing.side,
            (mpz_srcptr) doing.p, 
            (mpz_srcptr) doing.r, 
            side);
        ASSERT_ALWAYS(rc > 0);
        f = fopen(filename, "w");
        if (f == NULL) {
            perror("Error opening dumpfile");
        }
        free(filename);
    }
}

size_t dumpfile_t::write(const unsigned char * const data, const size_t size) const {
    if (!f) return 0;

    size_t const rc = fwrite(data, sizeof(unsigned char), size, f);
    ASSERT_ALWAYS(rc == size);
    return rc;
}
