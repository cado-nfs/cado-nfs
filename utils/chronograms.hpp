#ifndef CHRONOGRAMS_HPP_
#define CHRONOGRAMS_HPP_

#include <cstdint>

#include <array>
#include <string>
#include <type_traits>
#include <vector>

#include "params.hpp"
#include "timing.h"

namespace chronograms /* {{{ */ {
    extern void configure_switches(cxx_param_list & pl);
    extern void declare_usage(cxx_param_list & pl);
    extern void interpret_parameters(cxx_param_list & pl);
    extern bool is_enabled();

    // Empty tag structs for parameterless events
    struct NONE {};
    struct INIT {};
    struct QLATTICE {};
    struct SLICING {};
    struct ALLOC {};
    struct AB {};
    struct ECM {};
    struct DUPCHECK {};
    struct BOTCHED {};

    struct SSS   { int side; int level; };
    struct FIB   { int side; int level; uint32_t B; size_t slice; };
    struct DS    { int side; int level; uint32_t B; };
    struct PCLAT { int side; int level; size_t slice; };
    struct PBR   { int M; size_t B; };

    // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
    struct bubble_info {
        enum class kind_t : uint8_t {
            NONE,
            INIT,
            QLATTICE,
            SLICING,
            ALLOC,
            SSS,
            FIB,
            AB,
            DS,
            PCLAT,
            PBR,
            ECM,
            DUPCHECK,
            BOTCHED
        };

        static constexpr std::array<const char*, 14> kind_names = {
            "NONE",
            "INIT",
            "QLATTICE",
            "SLICING",
            "ALLOC",
            "SSS",
            "FIB",
            "AB",
            "DS",
            "PCLAT",
            "PBR",
            "ECM",
            "DUPCHECK",
            "BOTCHED", };

        kind_t kind;

        union payload {
            SSS   sss;
            FIB   fib;
            DS    ds;
            PCLAT pclat;
            PBR   pbr;
            DUPCHECK   dupcheck;

            payload() = default;
            explicit payload(SSS e)   : sss(e) {}
            explicit payload(FIB e)   : fib(e) { }
            explicit payload(DS e)    : ds(e) { }
            explicit payload(PCLAT e) : pclat(e) { }
            explicit payload(PBR e)   : pbr(e) { }
        } data = {};

        // NOLINTBEGIN(google-explicit-constructor,hicpp-explicit-conversions)
        bubble_info(NONE)       : kind(kind_t::NONE) {}
        bubble_info(INIT)       : kind(kind_t::INIT) {}
        bubble_info(SSS e)      : kind(kind_t::SSS), data(e) { }
        bubble_info(FIB e)      : kind(kind_t::FIB), data(e) { }
        bubble_info(DS e)       : kind(kind_t::DS), data(e) { }
        bubble_info(PCLAT e)    : kind(kind_t::PCLAT), data(e) { }
        bubble_info(PBR e)      : kind(kind_t::PBR), data(e) { }
	bubble_info(QLATTICE)   : kind(kind_t::QLATTICE) {}
	bubble_info(SLICING)    : kind(kind_t::SLICING) {}
	bubble_info(ALLOC)      : kind(kind_t::ALLOC) {}
	bubble_info(AB)         : kind(kind_t::AB) {}
	bubble_info(ECM)        : kind(kind_t::ECM) {}
	bubble_info(DUPCHECK)   : kind(kind_t::DUPCHECK) {}
	bubble_info(BOTCHED)    : kind(kind_t::BOTCHED) {}
        // NOLINTEND(google-explicit-constructor,hicpp-explicit-conversions)
    };
    // NOLINTEND(cppcoreguidelines-pro-type-union-access)

    struct bubble {
        uint64_t t0 = 0;
        uint64_t t1 = 0;
        uint64_t on_cpu = 0;
        bubble_info info = NONE {};

        static_assert(std::is_trivially_copyable_v<bubble_info>);

        bubble() = default;
        bubble(bubble_info info)
            : info(info)
        {}

        void start() {
            t0 = wct_nanoseconds();
            on_cpu = - microseconds_thread() * 1000;
        }
        void stop() {
            t1 = wct_nanoseconds();
            on_cpu += microseconds_thread() * 1000;
        }
    };

    class [[nodiscard]] bubble_guard {
        std::vector<chronograms::bubble> * destination = nullptr;
        bubble bubble_;

        public:
        static_assert(std::is_trivially_copyable_v<bubble_info>);
        bubble_guard(std::vector<chronograms::bubble> & destination, bubble_info info)
            : destination(&destination)
            , bubble_(info)
        {
            if (!chronograms::is_enabled()) return;
            bubble_.start();
        }

        ~bubble_guard() {
            if (!chronograms::is_enabled() || !destination) return;
            bubble_.stop();
            static_assert(std::is_trivially_copyable_v<chronograms::bubble_info>);
            destination->push_back(bubble_);
        }

        bubble_guard(const bubble_guard&) = delete;
        bubble_guard(bubble_guard&&) = delete;
        bubble_guard& operator=(const bubble_guard&) = delete;
        bubble_guard& operator=(bubble_guard&&) = delete;
        private:
        /* we intentionally make the default ctor private */
        bubble_guard() = default;
        public:
        static bubble_guard dummy() { return {}; }
    };

    std::string format_as(const bubble_info& info);
    void display(std::map<size_t, std::vector<chronograms::bubble>> const & M);
} /* namespace chronograms */ /* }}} */

#endif	/* CHRONOGRAMS_HPP_ */
