#include "cado.h"       // IWYU pragma: keep

#include <cstdint>

#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "nlohmann/json.hpp"
#include "fstream_maybe_compressed.hpp"
#include "params.hpp"

struct clamp_tool {
    parameter_with_default<std::string, "o", "output file", "-"> out;
    parameter_with_default<std::string, "i", "input file", "-"> in;
    parameter_mandatory<double, "t0", "begin time"> t0;
    parameter_mandatory<double, "t1", "begin time"> t1;

    static void configure(cxx_param_list & pl)
    {
        pl.declare_usage_section("General operational flags");

        decltype(out)::configure(pl);
        decltype(in)::configure(pl);
        decltype(t0)::configure(pl);
        decltype(t1)::configure(pl);
    }

    explicit clamp_tool(cxx_param_list & pl)
        : out(pl)
        , in(pl)
        , t0(pl)
        , t1(pl)
    {
    }

    void doit(std::istream & input, std::ostream& output) const
    {
        using json = nlohmann::json;
        json J;
        input >> J;

        if (J.contains("format")) {
            const uint64_t nt0 = t0() * 1.0e9;
            const uint64_t nt1 = t1() * 1.0e9;

            if (J["format"] != 20260727)
                throw std::runtime_error("Unexpected chronogram format");
            for(auto & [ k, v ] : J["events"].items()) {
                std::vector<json> clamped;
                for(auto const & w : v) {
                    const auto event_t0 = w[0].get<uint64_t>();
                    const auto event_t1 = event_t0 + w[1].get<uint64_t>();
                    if (event_t0 < nt1 && event_t1 > nt0)
                        clamped.push_back(w);
                }
                v = clamped;
            }
        } else {
            const double T0 = t0();
            const double T1 = t1();
            /* legacy (in seconds) */
            for(auto & [ k, v ] : J["events"].items()) {
                std::vector<json> clamped;
                for(auto const & w : v) {
                    const auto event_t0 = w[0].get<double>();
                    const auto event_t1 = event_t0 + w[1].get<double>();
                    if (event_t0 < T1 && event_t1 > T0)
                        clamped.push_back(w);
                }
                v = clamped;
            }
        }

        output << J;
    }
    void doit() const
    {
        std::istream * ptr_in = &std::cin;
        std::unique_ptr<std::istream> p_in;
        if (in() != "-") {
            p_in = std::unique_ptr<std::istream>(new ifstream_maybe_compressed(in()));
            ptr_in = p_in.get();
        }

        std::ostream * ptr_out = &std::cout;
        std::unique_ptr<std::ostream> p_out;
        if (out() != "-") {
            p_out = std::unique_ptr<std::ostream>(new ofstream_maybe_compressed(out));
            ptr_out = p_out.get();
        }

        doit(*ptr_in, *ptr_out);
    }
    
};

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    pl.declare_usage_header("This program is a help to deal with chronogram files\n");

    clamp_tool::configure(pl);

    pl.process_command_line(argc, argv);

    const clamp_tool p(pl);

    if (pl.warn_unused())
        pl.fail("not all parameters were parsed");


    p.doit();
}
