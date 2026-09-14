#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#include <string>
#include <vector>

#include "fmt/format.h"

#include "nlohmann/json.hpp"

#include "matmul-heatmap.hpp"
#include "fstream_maybe_compressed.hpp"
#include "gzip.h"
#include "verbose.hpp"

/* Insert tag before the extension of filename, keeping any compression
 * suffix last: "foo.json.gz" and "h0.v0" give "foo.h0.v0.json.gz"
 */
static std::string insert_tag(std::string const & filename,
        std::string const & tag)
{
    if (tag.empty())
        return filename;

    std::string const compression = get_suffix(filename, std::string());
    std::string stem = filename.substr(0, filename.size() - compression.size());

    auto const dot = stem.rfind('.');
    auto const slash = stem.rfind('/');

    if (dot != std::string::npos &&
            (slash == std::string::npos || dot > slash))
        return fmt::format("{}.{}{}{}",
                stem.substr(0, dot), tag, stem.substr(dot), compression);

    return fmt::format("{}.{}{}", stem, tag, compression);
}

void heatmap_dump(std::string const & filename,
        std::string const & tag,
        heatmap_info const & info,
        std::vector<heatmap_block> const & blocks)
{
    std::string const name = insert_tag(filename, tag);

    uint64_t const niter = uint64_t(info.iterations[0]) + info.iterations[1];

    if (!niter) {
        verbose_fmt_print(0, 0,
                "# Not writing {}: no iteration was timed\n", name);
        return;
    }

    using json = nlohmann::json;

    json J;
    J["format"] = 20260913;
    J["matrix"] = info.matrix;
    J["nrows"] = info.nrows;
    J["ncols"] = info.ncols;
    J["ncoeffs"] = info.ncoeffs;
    J["iterations"] = info.iterations;
    J["total"] = info.total / double(niter);
    J["total_max"] = info.total_max / double(niter);
    J["types"] = info.type_names;
    if (info.nh && info.nv) {
        J["grid"] = { info.nh, info.nv };
        J["submatrix"] = { info.submatrix_nrows, info.submatrix_ncols };
    }
    J["blocks"] = json::array();

    for(auto const & B : blocks) {
        json E;
        E.push_back(B.i0);
        E.push_back(B.i1);
        E.push_back(B.j0);
        E.push_back(B.j1);
        E.push_back(B.ncoeffs);
        E.push_back(B.type);
        E.push_back(B.dispatch / double(niter));
        E.push_back(B.combine / double(niter));
        J["blocks"].push_back(E);
    }

    ofstream_maybe_compressed out(name);
    out << J;

    verbose_fmt_print(0, 0,
            "# Heat map for {} ({} blocks) written to {}\n"
            "# It can be viewed with scripts/bwc-heatmap/heatmap.html\n",
            info.matrix, blocks.size(), name);
}
