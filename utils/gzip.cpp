#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio> // FILE // IWYU pragma: keep
#include <cstring>
#include <cerrno>

#include <array>
#include <fstream>  // filebuf
#include <ios>  // std::ios_base::openmode // IWYU pragma: keep
#include <iosfwd>   // filebuf too?
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#include <sys/wait.h>  // WIFEXITED WEXITSTATUS (on freebsd at least)
#include <unistd.h>     // close getpid
#include <sys/stat.h> // stat // IWYU pragma: keep
#ifdef HAVE_GETRUSAGE
#include <sys/time.h> // IWYU pragma: keep
#include <sys/resource.h> // IWYU pragma: keep
#endif

#include "fmt/base.h"
#include "fmt/format.h"

#include "macros.h"
#include "gzip.h"
#include "misc.h"
#include "cado_popen.h"
#include "cado_pipe_streambuf.hpp"
#include "portability.h" // realpath

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
static char antebuffer[PATH_MAX];       /* "directory/antebuffer" or "cat" */
static int antebuffer_buffer_size = 24; /* default value 2^24 = 16 Mo */
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

const std::array<suffix_handler, 6> supported_compression_formats {{
    { ".gz", "gzip -dc {}", "gzip -c --fast > {}", },
    { ".bz2", "bzip2 -dc {}", "bzip2 -c -1 > {}", },
    /* zstd seems to be uniformly better than any other alternative */
    { ".zstd", "zstd -dcf {}", "zstd --fast > {}", },
    /* xz is really slow */
    { ".xz", "xz -dc {}", "xz --fast > {}", },
    { ".lzma", "lzma -dc {}", "lzma -c -0 > {}", },
    /* These two have to be present */
    { "", "", "" },
}};

const char * path_basename(const char * path)
{
    const char *p = strrchr(path, '/');
    if (!p) {
        p = path;
    } else {
        p = p + 1;
    }
    return p;
}

int is_supported_compression_format(const char * s)
{
    for(auto const & r : supported_compression_formats) {
        if (r.suffix == s)
            return 1;
    }
    return 0;
}

int filename_matches_one_compression_format(const char * path)
{
    for(auto const & r : supported_compression_formats) {
        if (!*r.suffix) continue;
        if (has_suffix(path, r.suffix)) return 1;
    }
    return 0;
}

void get_suffix_from_filename (const char *s, char const **sfx)
{
  for(auto const & r : supported_compression_formats) {
    if (has_suffix(s, r.suffix))
    {
      *sfx = r.suffix;
      return;
    }
  }

  /* If we arrive here, it's because "" is not among the suffixes */
  abort();
}

static int try_antebuffer_path()
{
    int const rc = access(antebuffer, X_OK);
    if (rc >= 0) {
        fprintf(stderr, "antebuffer set to %s\n", antebuffer);
        return 1;
    }
    fprintf(stderr, "access to %s: %s\n", antebuffer, strerror(errno));
    *antebuffer = 0;
    return 0;
}

int set_antebuffer_path (const char *executable_filename, const char *path_antebuffer)
{
  *antebuffer = 0;
  antebuffer[PATH_MAX-1]='\0';
  /* First, if we have path_antebuffer, we must have antebuffer or error */
  if (path_antebuffer) {
      struct stat sbuf[1];
      int const rc = stat(path_antebuffer, sbuf);
      if (rc < 0) {
          fprintf(stderr, "%s: path_antebuffer=\"%s\" access error: %s\n",
                  __func__, path_antebuffer, strerror(errno));
      } else {
          /* Older versions had path_antebuffer be a directory. We still
           * support this, but only as a compatibility measure. */
          if (S_ISDIR(sbuf->st_mode)) {
#ifdef EXECUTABLE_SUFFIX
              snprintf(antebuffer, PATH_MAX-1, "%s/antebuffer" EXECUTABLE_SUFFIX, path_antebuffer);
#else
              snprintf(antebuffer, PATH_MAX-1, "%s/antebuffer", path_antebuffer);
#endif
          } else {
              strncpy(antebuffer, path_antebuffer, PATH_MAX-1);
          }
          if (try_antebuffer_path()) return 1;
      }
  }
  /* Second option: if we failed for any reason, and if $0 was given to
   * us, use that as a potential fallback */
  if (executable_filename) {
      char dummy[PATH_MAX];
      char dummy2[PATH_MAX + 64];
      const char * slash = strrchr(executable_filename, '/');
      if (slash) {
          int const len = MIN(PATH_MAX - 1, slash - executable_filename);
          strncpy(dummy, executable_filename, len);
          dummy[len]='\0';
      } else {
          dummy[0]='.';
          dummy[1]='\0';
      }
#ifdef EXECUTABLE_SUFFIX
      snprintf(dummy2, sizeof(dummy2), "%s/../utils/antebuffer" EXECUTABLE_SUFFIX, dummy);
#else
      snprintf(dummy2, sizeof(dummy2), "%s/../utils/antebuffer", dummy);
#endif
      if (realpath(dummy2, antebuffer) && try_antebuffer_path())
          return 1;
  }
  /* Third option: walk $PATH */
  if (path_resolve("antebuffer", antebuffer) && try_antebuffer_path()) {
      return 1;
  }
  *antebuffer = 0;
  fprintf(stderr, "No antebuffer configured\n");
  *antebuffer='\0';
  return 0;
}

/* Return a list of unix commands to _read_ a set of files. Consecutive
 * files sharing the same decompression mechanism are grouped into a
 * single command line.
 *
 * antebuffer file1.gz file2.gz file3.gz | gzip -dc
 * antebuffer file4.bz2 | gzip -dc
 * antebuffer file5.gz | gzip -dc
 * antebuffer file6.gz | cat    // useless use of cat, should be fixed.
 *
 * Note that antebuffer may also not be defined. In that case, the
 * simpler command formats like "gzip -dc file1.gz file2.gz file3.gz" are
 * used.
 *
 * All strings returned are meant to be passed to popen(). The return
 * value is a malloc()-ed array of malloc()-ed strings, and the caller is
 * in charge of freeing it (with filelist_clear, for instance).
 */
std::vector<std::string> prepare_grouped_command_lines(std::vector<std::string> const & list_of_files)
{
    std::vector<std::string> new_commands;
    
    /* Allow a few bytes extra for popen's "/bin/sh" "-c" prefix */
    ASSERT_ALWAYS(get_arg_max() >= 20);
    size_t const arg_max = get_arg_max() - 20;
    
    for(auto it = list_of_files.begin() ; it != list_of_files.end() ; ) {
        std::string cmd_prefix, cmd_postfix;
        const struct suffix_handler * this_suffix = nullptr;
        for (auto const & r : supported_compression_formats) {
            if (has_suffix(it->c_str(), r.suffix)) {
                this_suffix = &r;
                break;
            }
        }
        ASSERT_ALWAYS(this_suffix);

        if (*antebuffer) {
            if (*this_suffix->pfmt_in) {
                /* antebuffer 24 file1.gz file2.gz file3.gz | gzip -dc - */
                cmd_prefix  = fmt::format("{} {}", antebuffer, antebuffer_buffer_size);
                cmd_postfix = fmt::format(" | {}", fmt::format(fmt::runtime(this_suffix->pfmt_in), "-"));
            } else {
                /* antebuffer 24 file1.txt file2.txt file3.txt */
                /* avoid piping through cat */
                cmd_prefix  = fmt::format("{} {}", antebuffer, antebuffer_buffer_size);
            }
        } else {
            if (*this_suffix->pfmt_in) {
                /* gzip -dc file1.gz file2.gz file3.gz */
                cmd_prefix = fmt::format(fmt::runtime(this_suffix->pfmt_in), "");
            } else {
                /* cat file1.txt file2.txt file3.txt */
                /* There's potential for this to qualify as a useless use
                 * of cat, but anyway we don't expect to meet this case
                 * often.
                 */
                cmd_prefix = "cat ";
            }
        }
        
        std::string cmd = cmd_prefix;

        for( ; it != list_of_files.end() ; ++it) {
            const struct suffix_handler * other_suffix = nullptr;
            for (auto const & r : supported_compression_formats) {
                if (has_suffix(it->c_str(), r.suffix)) {
                    other_suffix = &r;
                    break;
                }
            }
            ASSERT_ALWAYS(other_suffix);
            if (other_suffix != this_suffix)
                break;

            cmd += " ";
            cmd += *it;

            if (cmd.size() + cmd_postfix.size() > arg_max)
                break;
        }

        /* Now all file names referenced by pointers in the interval
         * [grouphead..grouptail[ have the same suffix.
         */

        cmd += cmd_postfix;

        new_commands.push_back(cmd);
    }
    return new_commands;
}

FILE*
fopen_maybe_compressed2 (const char * orig_name, const char * mode, int* p_pipeflag, char const ** suf)
{
    FILE * f;

    std::string name = orig_name;

    // coverity[fs_check_call]
    if (strchr(mode, 'r') && access(name.c_str(), R_OK) != 0)
        return nullptr;

    for(auto const & r : supported_compression_formats) {
        if (!has_suffix(name.c_str(), r.suffix)) continue;
        if (suf) *suf = r.suffix;
        std::string command, tempname;

        /* Just *any* file that we write to will get a .tmp.$PID suffix
         */
        if (strchr(mode, 'w'))
            name = tempname = fmt::format("{}.tmp.{}", name, getpid());

        if (strchr(mode, 'r') && *r.pfmt_in)
            command = fmt::format(fmt::runtime(r.pfmt_in), name);
        else if (strchr(mode, 'w') && *r.pfmt_out)
            command = fmt::format(fmt::runtime(r.pfmt_out), name);

        if (!command.empty()) {
          /* apparently popen() under Linux does not accept the 'b' modifier */
            char pmode[2] = "x";
            pmode[0] = mode[0];
            f = cado_popen(command.c_str(), pmode);
            if (p_pipeflag) *p_pipeflag = 1;
#ifdef F_SETPIPE_SZxxx
            /* The pipe capacity is 2^16 by default; we can increase it,
             * but it does not seem to make a difference, thus we don't
             * change it by default (patch from Alain Filbois). */
            fcntl (fileno (f), F_SETPIPE_SZ, 1UL << 20);
#endif
        } else {
            f = fopen(name.c_str(), mode);
            if (p_pipeflag) *p_pipeflag = 0;
        }
        return f;
    }
    /* If we arrive here, it's because "" is not among the suffixes */
    abort();
    return nullptr;
}


FILE*
fopen_maybe_compressed (const char * name, const char * mode)
{
    return fopen_maybe_compressed2(name, mode, nullptr, nullptr);
}

#ifdef  HAVE_GETRUSAGE
int
fclose_maybe_compressed2 (FILE * f, const char * orig_name, struct rusage * rr)
#else
/* if we don't even have getrusage, then no fclose_maybe_compressed2 is
 * exposed. Yet, we use one as a code shortcut
 */
static int
fclose_maybe_compressed2 (FILE * f, const char * orig_name, void * rr MAYBE_UNUSED)
#endif
{
    std::string name = orig_name;

    for(auto const & r : supported_compression_formats) {
        if (!has_suffix(name.c_str(), r.suffix)) continue;
        /* It doesn't really make sense to imagine that one of these two
         * may exist and not the other */
        ASSERT_ALWAYS((!*r.pfmt_out) == (!*r.pfmt_in));

        const std::string tempname = fmt::format("{}.tmp.{}", name, getpid());

        int ret;

        if (*r.pfmt_in || *r.pfmt_out) {
#ifdef  HAVE_GETRUSAGE
            if (rr)
                ret = cado_pclose2(f, rr);
            else
#endif
                ret = cado_pclose(f);

#if defined(WIFEXITED) && defined(WEXITSTATUS)
            /* Unless child process finished normally and with exit ret 0,
               we return an error */
            if (ret == -1 || !WIFEXITED(ret) || WEXITSTATUS(ret) != 0) {
                return EOF;
            }
#else
            /* What do under MinGW? -1 definitely means an error, but how do
               we parse the other possible ret codes? */
            if (ret == -1)
                return EOF;
#endif

        } else {
#ifdef  HAVE_GETRUSAGE
            if (rr) memset(rr, 0, sizeof(*rr));
#endif
            ret = fclose(f);
            if (ret != 0)
                return ret;
        }

        /* do the rename only if the child completed successfully */

        // coverity[fs_check_call]
        struct stat sbuf[1];
        if (stat(tempname.c_str(), sbuf) == 0) {
            ret = rename(tempname.c_str(), name.c_str());
            if (ret != 0) return EOF;
        }

        return 0;
    }
    /* If we arrive here, it's because "" is not among the suffixes */
    abort();
    return EOF;
}

int
fclose_maybe_compressed (FILE * f, const char * name)
{
    return fclose_maybe_compressed2(f, name, nullptr);
}

streambase_maybe_compressed::streambase_maybe_compressed(std::string const & name, std::ios_base::openmode mode)
{
    open(name, mode);
    init(buf);
}

void streambase_maybe_compressed::open(std::string const & name_arg, std::ios_base::openmode mode)
{
    std::string name = name_arg;
    orig_name = name;
    if (mode & std::ios_base::out) {
        // fmtlib's fmt::format oddly mentions that it can throw a format
        // error, while its constexpr nature should be able to mark it as
        // impossible.
        // coverity[exception_thrown]
        tempname = fmt::format("{}.tmp.{}", name, getpid());
        name = tempname;
    }

    if (mode & std::ios_base::in && access(name.c_str(), R_OK) != 0)
        throw std::runtime_error("cannot open file for reading");
    /* creating is ok, of course
    if (mode & std::ios_base::out && access(name, W_OK) != 0)
        throw std::runtime_error("cannot open file for writing");
     */
    for(auto const & r : supported_compression_formats) {
        if (!has_suffix(orig_name.c_str(), r.suffix)) continue;
        std::string command;
        if (mode & std::ios_base::in && *r.pfmt_in)
            command = fmt::format(fmt::runtime(r.pfmt_in), name);
        else if (mode & std::ios_base::out && *r.pfmt_out)
            command = fmt::format(fmt::runtime(r.pfmt_out), name);

        if (!command.empty()) {
            /* apparently popen() under Linux does not accept the 'b' modifier */
            pbuf = std::make_unique<cado_pipe_streambuf>(command.c_str(), mode);
            buf = pbuf.get();
            pipe = true;
        } else {
            fbuf = std::make_unique<std::filebuf>();
            fbuf->open(name, mode);
            buf = fbuf.get();
            pipe = false;
        }
        /* hmmm */
        return;
    }
};

void streambase_maybe_compressed::close()
{
    if (pipe) pbuf->close();
    else fbuf->close();
    if (!tempname.empty()) {
        int const rc = rename(tempname.c_str(), orig_name.c_str());
        ASSERT_ALWAYS(rc == 0);
        tempname.clear();
    }
}

// we're in a dtor, exceptions can turn your computer into a coconut.
// yet we have an ASSERT_ALWAYS in close()
// coverity[exn_spec_violation]
streambase_maybe_compressed::~streambase_maybe_compressed()
{
    close();
}
