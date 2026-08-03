#include "cado.h"       // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"
#include "nlohmann/json.hpp"

#include "cxx_mpz.hpp"
#include "utils_cxx.hpp"
#include "macros.h"
#include "relation_cache.hpp"
#include "params.hpp"

using json = nlohmann::json;

std::string relation_cache::subdir_name(std::vector<unsigned long> const & split_q) const/*{{{*/
{
    std::string d;
    /* find the file */
    for(unsigned int i = 0 ; i + 1 < split_q.size() ; i++) {
        int l = 0;
        for(unsigned long s = 1 ; splits[i] > s ; s*=10, l++);
        d += fmt::format("/{:0{}}", split_q[i], l);
    }
    return d;
}/*}}}*/

static std::string find_filepath_inner(std::string const & d, unsigned long qq)/*{{{*/
{
    std::string filepath;
    DIR * dir = opendir(d.c_str());
    DIE_ERRNO_DIAG(dir == nullptr, "opendir(%s)", d.c_str());
    for(struct dirent * ent ; (ent = readdir(dir)) != nullptr ; ) {
        unsigned long q0, q1;
        if (sscanf(ent->d_name, "%lu-%lu", &q0, &q1) != 2) continue;
        if (qq < q0 || qq >= q1) continue;
        filepath = d + "/" + ent->d_name;
        break;
    }
    closedir(dir);

    return filepath;
}/*}}}*/

std::string relation_cache::find_filepath(cxx_mpz const & q0) const/*{{{*/
{
    cxx_mpz q = q0;
    std::vector<std::string> searched;

    /* write q in the variable basis given by the splits */
    std::vector<unsigned long> split_q = splits;
    for(unsigned int i = splits.size() ; i-- ; ) {
        split_q[i] = mpz_fdiv_ui(q, splits[i]);
        mpz_fdiv_q_ui(q, q, splits[i]);
    }
    if (mpz_cmp_ui(q, 0) != 0) {
        fmt::print(stderr, "# q is too large for relation cache\n", q0);
        exit(EXIT_FAILURE);
    }

    std::string d = path() + subdir_name(split_q);

    std::string filepath = find_filepath_inner(d, split_q.back());

    if (filepath.empty() && split_q.size() > 1) {
        searched.push_back(d);

        /* Try the previous directory, if qranges cross the
         * boundaries at powers of ten */
        split_q[split_q.size() - 2] -= 1;
        split_q[split_q.size() - 1] += splits[splits.size() - 1];
        d = path() + subdir_name(split_q);
        filepath = find_filepath_inner(d, split_q.back());
    }

    if (filepath.empty()) {
        searched.push_back(d);
        fmt::print(stderr, "# no file found in relation cache for q={} (searched directories: {})\n", q, join(searched, " "));
        exit(EXIT_FAILURE);
    }

    return filepath;
}
/*}}}*/


relation_cache::relation_cache(cxx_param_list & pl)
    : path(pl)
{
    if (!active())
        return;

    std::vector<unsigned long> splits;

    try {
        /* recover the list of splits from the config file */
        json dirinfo;
        if (!(std::ifstream(path() + "/dirinfo.json") >> dirinfo))
            throw std::exception();
        for(long s : dirinfo["splits"])
            splits.push_back(s);
    } catch (std::exception const & e) {
        fmt::print(stderr, "# Cannot read relation cache, or dirinfo.json in relation cache\n");
        exit(EXIT_FAILURE);
    }
}
