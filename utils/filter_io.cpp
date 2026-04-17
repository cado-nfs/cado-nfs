#include "cado.h" // IWYU pragma: keep

#include <string>
#include <vector>
#include <istream>

#include "ringbuf.hpp"
#include "timing.h"
#include "filter_io.hpp"
#include "utils_cxx.hpp"
#include "fstream_maybe_compressed.hpp"

void cado::filter_io_details::filter_rels_producer_thread(
    ringbuf & r,
    std::istream& f,
    timingstats_dict_ptr stats)
{
    r.feed_stream(f);
    if (!f.good() && !f.eof())
        throw cado::error("read error on input stream");

    r.mark_done();
    if (stats) timingstats_dict_add_mythread(stats, "producer");
}

void cado::filter_io_details::filter_rels_producer_thread(
    ringbuf & r,
    std::vector<std::string> const & input_files,
    timingstats_dict_ptr stats)
{
    for(auto const & filename : input_files) {
        ifstream_maybe_compressed f(filename);
        if (!f)
            throw cado::error("cannot open {}", filename);
        r.feed_stream(f);
        if (!f.good() && !f.eof())
            throw cado::error("read error on {}", filename);
    }
    r.mark_done();
    if (stats) timingstats_dict_add_mythread(stats, "producer");
}
