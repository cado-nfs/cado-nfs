#ifndef UTILS_FILELIST_HPP_
#define UTILS_FILELIST_HPP_

#include <vector>
#include <string>
#include <utility>

#include "params.hpp"

/* These are legacy calls */
extern char const ** filelist_from_file(const char * basepath, const char * filename,
                                  int typ);
extern void filelist_clear(char const ** filelist);


/* the interface still needs to be polished, obviously */
extern std::vector<std::string> filelist_from_file(std::string const & basepath, std::string const & filename, int typ);



struct filelist { /* {{{ */
    int argc;
    char const ** argv;

    /* The three parameters below specify the set of input files, of the
     * form <base path>/<one of the possible subdirs>/<one of the
     * possible file names>
     *
     * possible subdirs are lister in the file passed as subdirlist.
     * Ditto for possible file names.
     *
     * file names need not be basenames, i.e. they may contain directory
     * components. subdirlist and basepath are optional.
     */

    parameter<std::string,
        "filelist",
        "file containing a list of input files">
            my_files;
    parameter<std::string,
        "basepath",
        "path added to all file in filelist">
            basepath;
    parameter<std::string,
        "subdirlist",
        "file containing a list of subdirectories">
            subdirlist;

    /* path_antebuffer is actually stored elsewhere. Our goal is that
     * eventually, no code beyond this filelist interface remains tainted
     * by this antebuffer thing.
     */
    parameter<std::string,
        "path_antebuffer",
        "path to antebuffer program">
            path_antebuffer_loc;


    static void configure(cado::params::cxx_param_list & pl)
    {
        pl.declare_usage_section("input specification");
        decltype(my_files)::configure(pl);
        decltype(basepath)::configure(pl);
        decltype(subdirlist)::configure(pl);
        decltype(path_antebuffer_loc)::configure(pl);
    }

    filelist(cado::params::cxx_param_list & pl,
            int argc,
            char const **argv);

    /* build the list of input files from the given args
     * If no filelist is given, files are on the command-line.
     * If no subdirlist is given, files are easily construct from basepath and
     * filelist.
     * If subdirlist is given, it is a little bit trickier, see
     * filelist_from_file_with_subdirlist for more details. */
    std::vector<std::string> create_file_list() const;

    template<typename F>
    std::pair<std::vector<std::string>, std::vector<std::string>> separate_file_list(F const & f) const
    {
        std::vector<std::string> u, v;
        for(auto const & p : create_file_list())
            (f(p) ? u : v).push_back(p);
        return { u, v };
    }
}; /* }}} */

#endif	/* UTILS_FILELIST_HPP_ */
