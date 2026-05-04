#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <string>
#include <vector>
#include <ostream>
#include <fstream>
#include <iostream>

#include <gmp.h>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "params.hpp"
#include "utils_cxx.hpp"
#include "gmp_aux.h"

const std::vector<std::string> purged_sample = {
	"-27b241f,74c:7,9,c2,a3ee,12045,152c2,5,16,22,dd7,55d7,8c89,bc8d,be8c,1454a\n",
	"-228698d,6bc:7,12,33,37,1eba,23a0c,ad14e,10,1a,11b,44c2,8c89,16945,25fcd,2e768\n",
	"1d9e8d7,51f:4,196,af7,641b,13dc0,43082,0,1b,1f,47,93,fb0,8c89,17296,1ab65,30c92\n",
	"22935f5,538:9,1e,3c,6d,ca,1ee,65dc,1c6ca,e2,140,554,8c89,1731e,27462,105ac3\n",
	"-29a4211,73:4,37,138,1913,277d,358d,992f,a,19,121,12a,5f0,14b4,8c89,97c4,6b48c\n",
	"-ee0523,1bf:12,3c,48,4a5c,cb09,2c3b8,da,19f,176a,73ca,8c89,c441,5d2bd,1\n",
	"-1432aa7,ec:531,63f,2be6,3c16,eab88,5,14,6c,86,bfc,6684,8c89,d555,4140b\n",
	"1b94c99,b:b2,ed,16b,365,4d7a,d9cdc,13,19,5a3,11a0,1578,1750,227c,2bb4,8c89\n",
	"f49f9,19c:9,c,1c,45,6d,1128,1cfa2,5f3b0,5,10,16,30,109,2f7,1f6f,55af,8c89,4a8fb\n",
	"-4600fb,28:7,c,84,2af,3d2,5082d,51b55,581,8af,e1b,8c89,9478,9bffa\n",
	"2ddaf63,2a6:1d,21,58,bf,16c,1379b,6a612,5,2,f,16,18,27,8b,233,5f0,10d5,8c89,50344,1,3\n",
};

const std::vector<std::string> sieve_sample  {
	"181011149,3191:2,2,7,133,4513,1c433,267e3,49410b:3,3,5,5,5,b,b,d,7f,2a5,b5d,2889b,64205,148b3b,230a7f\n",
	"156003139,4096:17,89,337,5c9,1079,1d011,66f271:2,2,3,3,5,5,b,d,d,d,11,2b,e4d3,50de9,5730b,64205,311b87\n",
	"172515198,2611:3,3,3,43,67,2ab,26029,3f289,245c17:2,11,35,4f,4a89,2330b,2b373,4d87b,64205,52a951\n",
	"147063796,1295:59,199,6f7,39bf,127dd1,7b0019:2,5,b,b,b,d,1d,ef,139,599,6e3b,8b07,64205,5d927d\n",
	"104349541,2463:2,2,119,128f,3fc6d,10da41,52efa7:3,3,3,5,5,d,d,29,26b3,3bcb,8047,64205,11ea3f,1e6be3\n",
	"106006947,2435:2,3,3,7,b,d,11,34f1,91ed,2545f,56f2f:5,29,a7,fa3,64205,7e151,c294b,e5bdf,f954b\n",
	"83543671,2281:2,2,2,5,b,83,4a3,529,1fe1,12a5b,1e8cb:3,3,29,3b,125,20b,2ec3,12811,64205,1f3dab,4ef3a1\n",
	"56909413,1855:2,2,17,25,407,318dd,17364d,485a89:3,3,5,5,b,d,11,35,5fcf,64205,c4837,39dd83,49afe9\n",
	"42911531,1363:2,25,25,71,bf,60d,a49,5453,b3ef9:5,a7,10f,aed,dd5f,3ad65,64205,8e225,9ad81\n",
	"30942432,1585:3,69d,2405,a493,4aaaf,56e96d:2,5,b,d,d,11,4f,3f3d,64205,a85b3,1269a3,6d2305\n",
	"29111011,1029:2,2,2,2,47f,1cb7,32f9,45c23,522959:3,3,3,5,5,d,25,101,1bb,209,a93,64205,249b99,694487\n",
	"14697633,1127:2,2,2,2,2,2,3,1d,a39,c6eb,5f251,6c3055:5,5,d,a7,b9b,2d9d,18d31,64205,c0c77,1739ef\n",
	"-5507682,1405:3,3,3,3,5615,14b49,169c13,567f4f:2,5,5,11,265,1567,224f,173d7,36a25,64205,16cc21\n",
	"-31712453,608:7,7,29,89,26b9d,3e815d,57efe7:2,2,5,5,1d,61,4d5,c11,3611f,5b219,64205,2f0725\n",
	"16862664,163:3,3,11,47,175,1949,3a51,a09d,108af:2,83,209,67df,1ef29,64205,7da5b,625567\n",
	"-23247059,270:13,35,3d,49,ef,a99,7713d,a3d57:2,3,3,5,b,d,2b,125,14b,1505,2fcf,3e05,64205,a0ae9\n",
	"14820548,149:287,3d7b,9efb,729ef,4bb0f3:2,3,3,5,d,11,11,47,199,355,3605,1190b,22f5d,64205\n",
	"6840977,265:2,11,59,6c77,203db,28745,583eb:5,b,29,35,9d,1607,191b,64205,ee12f,4a2d3b\n",
	"-14240393,490:2b,df,10d,aa9,1a17,4505,40d47:2,3,3,5,b,d,11,83,8ef,4ac47,64205,7ea93,173143\n",
};

struct generator {
    parameter_with_default<size_t,
        "target-size",
        "Approximate size of the file to generate, in bytes",
        "65536">
            target_size;
    parameter<std::string,
        "o",
        "output file name (default stdout)">
            outname;
    parameter_with_default<std::string,
        "mode",
        R"(type of relations to generate ("purged" or "sieve"))",
        "purged">
            type;

    size_t target_nrels;

    static void configure(cxx_param_list & pl)
    {
        decltype(target_size)::configure(pl);
        decltype(outname)::configure(pl);
        decltype(type)::configure(pl);
    }

    explicit generator(cxx_param_list & pl)
        : target_size(pl)
        , outname(pl)
        , type(pl)
    {
        if (!psource())
            pl.fail(R"(--mode should be either "purged" or "sieve")");

        auto const & S = source();

        size_t s = 0;
        for(auto const & t : S) s += t.size();

        target_nrels = static_cast<size_t>(double_ratio(
                    target_size.parameter_value() * S.size(),
                    s));

        /*
        fmt::print(stderr,
                "# for a file of approximately {} bytes,"
                " we need roughly {} relations from our sample\n",
                target_size.parameter_value(), target_nrels);
                */
    }

    std::vector<std::string> const * psource() const {
        if (type.parameter_value() == "purged")
            return &purged_sample;
        else if (type.parameter_value() == "sieve")
            return &sieve_sample;
        return nullptr;
    }
    std::vector<std::string> const & source() const {
        return *psource();
    }

    void generate(std::ostream& o)
    {
        if (type.parameter_value() == "purged") {
            /* 0x105ac4 is one past the last encountered column index in
             * our sample.
             * 143 is very abnormal, since our sample comprises only
             * repetitions of the same relations over and over again.
             */
            fmt::print(o, "# {} {} {}\n",
                    target_nrels,
                    0x105ac4, 143);
        }
        auto const & S = source();
        const size_t K = S.size();
        cxx_gmp_randstate R;
        for(size_t i = 0 ; i < target_nrels ; i++)
            o << S[gmp_urandomm_ui(R, K)];
    }

    void generate()
    {
        if (outname.is_provided()) {
            std::ofstream f(outname);
            generate(f);
        } else {
            generate(std::cout);
        }
    }

};

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    generator::configure(pl);

    pl.process_command_line(argc, argv);
    
    generator G(pl);

    if (pl.warn_unused())
        pl.fail("Unused arguments were given");

    
    G.generate();
}

