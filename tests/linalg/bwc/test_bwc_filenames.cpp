#include "cado.h" // IWYU pragma: keep

#include "fmt/base.h"

#include "bwc_filenames.hpp"

int main()
{
    for(const char * name: {
        "/tmp/cado.4nt73is9/c30.bwc/S.sols0-64.0-64",
        "/tmp/cado.4nt73is9/c30.bwc/F.sols0-64.0-64",
        "/tmp/cado.4nt73is9/c30.bwc/V0-64.64",
        "/tmp/cado.4nt73is9/c30.bwc/A0-64.0-64",
        "/tmp/cado.4nt73is9/c30.bwc/Cd0-64.80",
        "/tmp/cado.4nt73is9/c30.bwc/Cv0-64.80",
        "/tmp/cado.4nt73is9/c30.bwc/Cd0-64.64",
        "/tmp/cado.4nt73is9/c30.bwc/Cv0-64.64",
        "/tmp/cado.4nt73is9/c30.bwc/Cd0-64.16",
        "/tmp/cado.4nt73is9/c30.bwc/Cv0-64.16",
        "/tmp/cado.4nt73is9/c30.bwc/Cv0-64.0",
        "/tmp/cado.4nt73is9/c30.bwc/Ct0-64.0-64",
        "/tmp/cado.4nt73is9/c30.bwc/Cr0-64.0-64",
        "/tmp/cado.4nt73is9/c30.bwc/X",
        "/tmp/cado.4nt73is9/c30.bwc/V0-64.0",
        "/tmp/cado.4nt73is9/c30.bwc/c30.sparse.1x2",
        "/tmp/cado.4nt73is9/c30.bwc/c30.sparse.1x2/c30.sparse.1x2.2707664384.h0.v1-bucketT.bin",
        "/tmp/cado.4nt73is9/c30.bwc/c30.sparse.1x2/c30.sparse.1x2.2707664384.h0.v0-bucketT.bin",
        "/tmp/cado.4nt73is9/c30.bwc/c30.sparse.1x2/c30.sparse.1x2.bin",
        /* Add a few extras, including near-matches */
        "A0-64.0-64.txt",
        "VV0-64.0",
        "/VV0-64.0",
        "W0-64.0",
        })
    {
        bwc_V_file v;
        bwc_S_file s;
        bwc_Cv_file cv;
        bwc_Cd_file cd;
        bwc_Cr_file cr;
        bwc_Ct_file ct;
        bwc_A_file a;
        bwc_F_file f;

        if (bwc_V_file::match(v, name)) {
            fmt::print("{}: recognized V file,"
                    " columns {} to {}, iteration {}\n",
                    name, v.j0, v.j1, v.n);
        } else if (bwc_S_file::match(s, name)) {
            fmt::print("{}: recognized S file,"
                    " solutions {} to {}, iterations {} to {}\n",
                    name, s.s0, s.s1, s.n0, s.n1);
        } else if (bwc_Cv_file::match(cv, name)) {
            fmt::print("{}: recognized Cv file,"
                    " columns {} to {}, stretch {}\n",
                    name, cv.j0, cv.j1, cv.stretch);
        } else if (bwc_Cd_file::match(cd, name)) {
            fmt::print("{}: recognized Cd file,"
                    " columns {} to {}, stretch {}\n",
                    name, cd.j0, cd.j1, cd.stretch);
        } else if (bwc_Cr_file::match(cr, name)) {
            fmt::print("{}: recognized Cr file,"
                    " {} checks\n",
                    name, cr.nchecks);
        } else if (bwc_Ct_file::match(ct, name)) {
            fmt::print("{}: recognized Ct file,"
                    " {} checks, m={}\n",
                    name, ct.nchecks, ct.m);
        } else if (bwc_A_file::match(a, name)) {
            fmt::print("{}: recognized A file,"
                    " columns {} to {}, iterations {} to {}\n",
                    name, a.j0, a.j1, a.n0, a.n1);
        } else if (bwc_F_file::match(f, name)) {
            fmt::print("{}: recognized F file,"
                    " solutions {} to {}, columns {} to {}\n",
                    name, f.s0, f.s1, f.j0, f.j1);
        } else {
            fmt::print("{}: no file name recognized\n", name);
        }
    }
    return 0;
}

