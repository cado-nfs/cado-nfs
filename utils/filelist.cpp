#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <string>

#include "macros.h"
#include "gzip.h"
#include "filelist.hpp"
#include "utils_cxx.hpp"
#include "filelist.hpp"

#include "portability.h" // asprintf // IWYU pragma: keep

static void chomp(char *s) {
    char *p;
    if (s && (p = strrchr(s, '\n')) != NULL)
        *p = '\0';
}

/* Return a NULL-terminated list of file names read from filename.
   Empty lines and comment lines (starting with '#') are skipped.
   If basepath != NULL, it is used as path before each read filename
*/
char const ** filelist_from_file(const char * basepath, const char * filename,
                           int typ)
{
    char const ** files = NULL;
    int nfiles_alloc = 0;
    int nfiles = 0;
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
      if (typ == 0)
        perror ("Problem opening filelist");
      else
        perror ("Problem opening subdirlist");
      exit (1);
    }
    char relfile[FILENAME_MAX + 10];
    while (fgets(relfile, FILENAME_MAX + 10, f) != NULL) {

        // skip leading blanks
        char *rfile = relfile;
        while (isspace((int)(unsigned char)rfile[0]))
            rfile++;
        // if empty line or comment line, continue
        if ((rfile[0] == '#') || (rfile[0] == '\0') || (rfile[0] == '\n'))
            continue;
        chomp(rfile);

        if (nfiles == nfiles_alloc) {
            nfiles_alloc += nfiles_alloc / 2 + 16;
            checked_realloc(files, nfiles_alloc);
        }
        if (basepath) {
            char * name;
            int ret = asprintf(&name, "%s/%s", basepath, rfile);
            ASSERT_ALWAYS(ret >= 0);
            files[nfiles] = name;
        } else {
            files[nfiles] = strdup(rfile);
        }
        nfiles++;
    }
    fclose(f);

    if (nfiles == nfiles_alloc) {
        nfiles_alloc += nfiles_alloc / 2 + 16;
        checked_realloc(files, nfiles_alloc);
    }
    files[nfiles++] = NULL;
    return files;
}

void filelist_clear(char const ** filelist)
{
    if (!filelist) return;
    for(char const ** p = filelist ; *p ; p++)
        free((char *) *p);
    free(filelist);
}

/* Build the file list (ugly). It is the concatenation of all
 *  b s p
 * where:
 *    b is the basepath (empty if not given)
 *    s ranges over all subdirs listed in the subdirlist (empty if no
 *    such list)
 *    p ranges over all paths listed in the filelist.
 */
static char const ** filelist_from_file_with_subdirlist(const char *basepath,
					  const char *filelist,
					  const char *subdirlist)
{
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char const ** fl = filelist_from_file(NULL, filelist, 0);
    for (char const ** p = fl; *p; p++, nfiles++);

    char const ** sl = filelist_from_file(basepath, subdirlist, 1);
    for (char const ** p = sl; *p; p++, nsubdirs++);

    char ** fic = (char ** ) malloc((nsubdirs * nfiles + 1) * sizeof(char *));
    ASSERT_ALWAYS(fic != NULL);

    char ** full = fic;
    for (char const ** f = fl; *f; f++) {
	for (char const ** s = sl; *s; s++, full++) {
	    int ret = asprintf(full, "%s/%s", *s, *f);
	    ASSERT_ALWAYS(ret >= 0);
	}
    }
    *full = NULL;
    filelist_clear(fl);
    filelist_clear(sl);
    return (char const **) fic;
}


std::vector<std::string> filelist_from_file(std::string const & basepath, std::string const & filename, int typ)
{
    auto * files = filelist_from_file(basepath.c_str(), filename.c_str(), typ);
    std::vector<std::string> res;
    for (const auto * p = files; *p; p++) {
        res.emplace_back(*p);
    }
    filelist_clear(files);
    return res;
}

filelist::filelist(cxx_param_list & pl, /* {{{ */
            int argc,
            char const **argv)
        : argc(argc)
        , argv(argv)
        , my_files(pl)
        , basepath(pl)
        , subdirlist(pl)
        , path_antebuffer_loc(pl)
{
    if (basepath.is_provided() && my_files.is_default())
        throw parameter_error("-basepath only valid with -filelist");
    if (subdirlist.is_provided() && my_files.is_default())
        throw parameter_error("-basepath only valid with -filelist");
    if (my_files.is_provided() + (argc != 0) != 1)
        throw parameter_error("provide either -filelist or freeform file names");
    if (path_antebuffer_loc.is_provided())
        set_antebuffer_path(pl.binary_name().c_str(),
                path_antebuffer_loc.parameter_value().c_str());
} /* }}} */

/*{{{ build the list of input files from the given args
 * If no filelist is given, files are on the command-line.
 * If no subdirlist is given, files are easily construct from basepath and
 * filelist.
 * If subdirlist is given, it is a little bit trickier, see
 * filelist_from_file_with_subdirlist for more details. */
std::vector<std::string> filelist::create_file_list() const
{
    if (my_files.is_default())
        return { argv, argv + argc };

    std::vector<std::string> res;

    char const ** ifiles;
    if (subdirlist.is_default())
        ifiles = filelist_from_file(basepath.parameter_value().c_str(),
                my_files.parameter_value().c_str(), 0);
    else
        ifiles = filelist_from_file_with_subdirlist(
                basepath.parameter_value().c_str(),
                my_files.parameter_value().c_str(),
                subdirlist.parameter_value().c_str());

    for (const auto * p = ifiles; *p; p++)
        res.emplace_back(*p);

    /* Free allocated stuff */
    filelist_clear(ifiles);
    return res;
} /*}}}*/
