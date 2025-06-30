#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <sstream>
#include <algorithm>
#include <string>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fmt/base.h"
#include "fmt/format.h"

#include "rolling.hpp"
#include "bw-common.h"
#include "macros.h"
#include "portability.h"
#include "utils_cxx.hpp"


void keep_rolling_checkpoints(std::string const & stem, unsigned int v)
{
    if (bw->keep_rolling_checkpoints == 0)
        return;

    std::vector<unsigned int> vs;

    DIR * d = opendir(".");
    for(struct dirent * de; (de = readdir(d)) != NULL ; ) {
        unsigned int k;
        std::istringstream is(de->d_name);
        if (!(is >> expect(stem) >> expect(".") >> k))
            continue;
        if (v && k > v)
            continue;
        vs.push_back(k);
    }
    closedir(d);

    if (vs.size() <= (size_t) bw->keep_rolling_checkpoints) {
        return;
    }
    std::sort(vs.begin(), vs.end());
    vs.erase(vs.end() - bw->keep_rolling_checkpoints, vs.end());
    for(unsigned int const k : vs) {
        if (bw->checkpoint_precious && (k % bw->checkpoint_precious == 0))
            continue;
        if (k == 0)
            continue;
        std::string v = fmt::format("{}.{}", stem, v);
        struct stat sbuf[1];
        int rc = stat(v.c_str(), sbuf);
        if (rc < 0) {
            if (errno == ENOENT) {
                fmt::print("Old checkpoint {} is gone already\n", v);
            } else {
                fmt::print("Old checkpoint {}: {}\n", v, strerror(errno));
            }
        } else {
            ASSERT_ALWAYS(rc == 0);
            time_t const now = time(nullptr);
            int age = now - sbuf->st_mtime;
            if (age < bw->keep_checkpoints_younger_than) {
                fmt::print("Not discarding old checkpoint {}, too recent ({} s < {})\n", v, age, bw->keep_checkpoints_younger_than);
            } else {
                fmt::print("Discarding old checkpoint {}\n", v);
                rc = unlink(v.c_str());
                if (rc < 0) perror(v.c_str());
            }
        }
    }
}

