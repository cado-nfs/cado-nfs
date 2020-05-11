#ifndef LAS_DUMPFILE_HPP_
#define LAS_DUMPFILE_HPP_

#include <cstdio>  // for size_t, FILE, NULL
struct las_todo_entry;

/* If -dumpregion is given, dump sieve region to a file to be able to
   compare new sieving code with a known good reference. Beware of
   resulting large files. */
class dumpfile_t {
    FILE *f = NULL;
public:
    ~dumpfile_t();
    void close();
    void open(const char *filename_stem, las_todo_entry const & doing, int side);
    size_t write(const unsigned char *, size_t) const;
};


#endif	/* LAS_DUMPFILE_HPP_ */
