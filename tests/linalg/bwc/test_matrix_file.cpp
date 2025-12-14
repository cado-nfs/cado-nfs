#include "cado.h" // IWYU pragma: keep

#include <cstring>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

#include "matrix_file.hpp"
#include "sha1.h"
#include "params.h"

static std::ostream * bind(cxx_param_list & pl, const char * option_name, std::shared_ptr<std::ostream>& b)
{
    const char * filename = param_list_lookup_string(pl, option_name);

    if (!filename)
        return nullptr;

    if (strcmp(filename, "-") == 0)
        return &std::cout;

    auto * p = new std::ofstream(filename, std::ios::out | std::ios::binary);
    if (!*p) throw std::runtime_error(std::string(filename) + ": cannot open");
    b.reset(p);

    return p;
}


int main(int argc, char const * argv[])
{
    cxx_param_list pl;
    argv++, argc--;
    std::vector<std::string> wild;

    int direction = 0;

    param_list_configure_switch(pl, "--columns", &direction);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild.size() < 3) {
            wild.emplace_back(argv[0]);
            argv++, argc--;
            continue;
        }
    }

    if (wild.empty())
        throw std::runtime_error("missing verb");

    matrix_file M(wild[1]);
    M.lookup();

    std::shared_ptr<std::ostream> bm;
    std::shared_ptr<std::ostream> bc;
    std::shared_ptr<std::ostream> bo;
    std::shared_ptr<std::ostream> bd;

    std::ostream *pbm = bind(pl, "binary-mixed", bm);
    std::ostream *pbc = bind(pl, "binary-coeffs", bc);
    std::ostream *pbo = bind(pl, "binary-offsets", bo);
    std::ostream *pbd = bind(pl, "binary-data", bd);


    if (wild[0] == "read") {
        if (wild.size() == 2) {
            M.read(direction);
            if (pbm) M.dump_mixed(*pbm);
            if (pbd) M.dump_data(*pbd);
            if (pbo) M.dump_offsets(*pbo);
            if (pbc && M.withcoeffs) M.dump_coeffs(*pbc);
            if (!(pbm || pbo || pbd || pbc)) {
                {
                    sha1_checksumming_stream sbuf;
                    M.dump_mixed(sbuf);
                    std::cout << "bin " << sbuf.digest() << "\n";
                }
                {
                    sha1_checksumming_stream sbuf;
                    M.dump_data(sbuf);
                    std::cout << M.rowcol[direction] << "s " << sbuf.digest() << "\n";
                }
                {
                    sha1_checksumming_stream sbuf;
                    M.dump_offsets(sbuf);
                    std::cout << M.rowcol[direction] << "_offsets " << sbuf.digest() << "\n";
                }
                if (M.withcoeffs) {
                    sha1_checksumming_stream sbuf;
                    M.dump_coeffs(sbuf);
                    std::cout << M.rowcol[direction] << "_coeffs " << sbuf.digest() << "\n";
                }
            }
        } else if (wild.size() == 3) {
            throw std::runtime_error("read with bfile is collective, test is not ready yet");
        } else {
            throw std::runtime_error("read requires matrix name");
        }
    }
}
