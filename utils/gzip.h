#ifndef CADO_GZIP_H
#define CADO_GZIP_H

#include "cado_config.h"  // for HAVE_GETRUSAGE
#include <stdio.h>

#ifdef __cplusplus
#include <vector>
#include <string>
#endif

#ifdef  HAVE_GETRUSAGE
#include <sys/time.h>   // IWYU pragma: keep
#include <sys/resource.h>       // IWYU pragma: keep
#endif

#ifdef __cplusplus
/* Return a unix commands list with antebuffer. Example:
 * antebuffer X file_relation1 | cat -
 * antebuffer X file_relation2.gz file_relation3.gz | gzip -dc -
 * antebuffer X file_relation4.bz2 file_relation5.bz2 | bzip2 -dc -
 * [empty string]

 * antebuffer_cmd is /path/to/antebuffer <X>, where <X> is an integer
 * denoting the size of the antebuffer (24 means 2^24 bytes).
 */
std::vector<std::string> prepare_grouped_command_lines(std::vector<std::string> const & list_of_files);
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Search the executable in PATH and, if found, return in real_path the
   complete path WITH the executable in the end */
char * search_real_exec_in_path(const char *executable, char *real_path);

/* Search the path of antebuffer and put the complete path + name of the
 * executable "antebuffer" in a static variable used by the gzip.c layer.
 * Arguments are the $0 variable (path from cwd to the current executable
 * file), and path_antebuffer, which is obtained from the command line if
 * it happens to exist. Either, or both, may be NULL. The rule for
 * deriving the path to the antebuffer binary is as follows (first match
 * wins):
 *
 *  - if path_antebuffer is a complete path to an executable program, use
 *  it.
 *  - if executable_filename is non-NULL, use `dirname
 *  $0`/../utils/antebuffer, if that happens to point to a valid
 *  executable filename.
 *  - otherwise, antebuffer is disabled.
 *
 * Configuring to *not* use antebuffer can be done by calling
 * set_antebuffer_path(NULL, NULL), or not calling it at all.
 *
 * The return value is 1 if an executable has been found.
 *
 * This is a configuration function which must be called at most once, and in
 * a monothreaded context.
 */
int set_antebuffer_path (const char *executable_filename, const char *path_antebuffer);

/* There are of course scores of existing basename() codes accessible,
 * starting with POSIX basename. However we fear posible inconsistencies
 * here, so we stick to a simple-and-stupid version, whose specification
 * meets our needs. Here, the returned string is always a substring of
 * the input string, and the latter never undergoes any modification.
 */
extern const char * path_basename(const char * path);

extern int is_supported_compression_format(const char * s);
extern int filename_matches_one_compression_format(const char * path);

/* Put in sfx the suffix in s (can be "" or NULL) */
extern void get_suffix_from_filename (const char *s, char const **sfx);

/* Takes a filename, possibly ending with any recognized compression
 * extension, and returns the associated file stream. The stream may have
 * been opened by either fopen() of popen(), depending on whether an
 * external decompression program has to be called. If the caller is
 * intersted in knowing, the integer *p_pipeflag is filled with 1 (for
 * popen) or 0 (for fopen). In fact, the caller should care, because this
 * can be used to decide whether to close the stream with pclose or
 * fclose (even if fclose works, it's almost guaranteed to create
 * zombies).
 * If non-NULL, suf is the location of a pointer which is position to the
 * recognized suffix, which has been used to decide on which compression
 * method.
 */
extern FILE * fopen_maybe_compressed2(const char * name, const char * mode, int* p_pipeflag, char const ** suf);
extern FILE * fopen_maybe_compressed(const char * name, const char * mode);

/* This one just looks at the file name, and guesses again whether popen() or
 * fopen() was used. The file stream is then closed with pclose() or
 * fclose() accordingly.  */
extern int fclose_maybe_compressed(FILE *, const char * name);

#ifdef  HAVE_GETRUSAGE
/* Same, but recovers the time taken by the underlying process */
extern int fclose_maybe_compressed2 (FILE * f, const char * name, struct rusage * r);
#endif

#ifdef __cplusplus
namespace cado::details {
struct suffix_handler {
    char const * suffix;
    char const * pfmt_in;
    char const * pfmt_out;
};

static constexpr suffix_handler supported_compression_formats[] = {
    {
        .suffix = ".gz",
        .pfmt_in = "gzip -dc {}",
        .pfmt_out = "gzip -c --fast > {}",
    },
    {
        .suffix = ".bz2",
        .pfmt_in = "bzip2 -dc {}",
        .pfmt_out = "bzip2 -c -1 > {}",
    },
    {
        /* zstd seems to be uniformly better than any other alternative */
        .suffix = ".zstd",
        .pfmt_in = "zstd -dcf {}",
        .pfmt_out = "zstd --fast > {}",
    },
    {
        /* xz is really slow */
        .suffix = ".xz",
        .pfmt_in = "xz -dc {}",
        .pfmt_out = "xz --fast > {}",
    },
    {
        .suffix = ".lzma",
        .pfmt_in = "lzma -dc {}",
        .pfmt_out = "lzma -c -0 > {}",
    },

    {
        /* must be present: we must have "" among the suffixes because
         * that is the always-true match */
        .suffix = "",
        .pfmt_in = "",
        .pfmt_out = ""
    },
};

} /* namespace cado::details */
#endif

#ifdef __cplusplus
}
#endif


#endif	/* CADO_GZIP_H */
