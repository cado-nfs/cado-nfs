#ifndef CADO_LAS_DUMPFILE_HPP
#define CADO_LAS_DUMPFILE_HPP

#include <cstdio>  // for size_t, FILE, NULL
struct special_q;

/* If -dumpregion is given, dump sieve region to a file to be able to
   compare new sieving code with a known good reference. Beware of
   resulting large files. */
class dumpfile_t {
    FILE *f = NULL;
public:
    ~dumpfile_t();
    void close();
    void open(const char *filename_stem, special_q const & doing, int side);
    size_t write(const unsigned char *, size_t) const;
};


#endif	/* CADO_LAS_DUMPFILE_HPP */
