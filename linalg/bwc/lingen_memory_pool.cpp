#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <memory>
#include <utility>

#ifdef HAVE_EXECINFO
#include <execinfo.h>
#ifdef HAVE_CXXABI_H
/* We use that to demangle C++ names */
#include <cxxabi.h>
#endif
#endif

#include "lingen_memory_pool.hpp"
#include "misc.h"
#include "select_mpi.h"
#include "utils_cxx.hpp"

void memory_pool_details::inaccuracy_handler<true>::handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed)
{
    if (already_allocated + asked > previously_allowed) {
        size_t const d = (already_allocated + asked) - previously_allowed;
        {
            char buf[20];
            /*
               if (d < 0 && (size_t) -d > max_inaccuracy) {
               fprintf(stderr, "# Over-estimating the amount of reserved RAM by %asked\n",
               size_disp((size_t) -d, buf));
               max_inaccuracy = -d;
               } else if (d > 0 && (size_t) d > max_inaccuracy) {
               fprintf(stderr, "# Under-estimating the amount of reserved RAM by %asked\n",
               size_disp((size_t) d, buf));
               max_inaccuracy = d;
               }
               */
            cumulated_inaccuracy += d;
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            fprintf(stderr, "# Under-estimating the amount of reserved RAM by %s\n",
                    size_disp(cumulated_inaccuracy, buf));
        }

        previously_allowed += d;
    }
}

void memory_pool_details::alloc_check(const char * text, bool condition)
{
    if (condition) return;
    throw memory_pool_exception(text);
}

/* copied from las_debug.cpp
 * TODO: refactor? */

#ifdef HAVE_EXECINFO
static std::string remove_trailing_address_suffix(std::string const& a, std::string& suffix)
{
    size_t const pos = a.find('+');
    if (pos == std::string::npos) {
        suffix.clear();
        return a;
    }
    suffix = a.substr(pos);
    return a.substr(0, pos);
}

static std::string get_parenthesized_arg(std::string const& a, std::string& prefix, std::string& suffix)
{
    size_t const pos = a.find('(');
    if (pos == std::string::npos) {
        prefix=a;
        suffix.clear();
        return {};
    }
    size_t const pos2 = a.find(')', pos + 1);
    if (pos2 == std::string::npos) {
        prefix=a;
        suffix.clear();
        return {};
    }
    prefix = a.substr(0, pos);
    suffix = a.substr(pos2 + 1);
    return a.substr(pos+1, pos2-pos-1);
}
#endif

memory_pool_exception::memory_pool_exception(std::string message)
    : message(std::move(message))
{
#ifdef HAVE_EXECINFO
    int sz = 100;

    void * callers_addresses[sz];
    sz = backtrace(callers_addresses, sz);
    std::unique_ptr<char *[], free_delete<char *>> callers;
    callers.reset(backtrace_symbols(callers_addresses, sz));

    message += "\n";
    message += "======= Backtrace: =========\n";

    for (int i = 0; i < sz; i++) {
        std::string caller = callers[i];
        std::string xx,yy,zz;
        yy = get_parenthesized_arg(caller, xx, zz);
        if (!yy.empty()) caller = yy;

        if (caller.empty()) {
            caller="<no symbol (static?)>";
        } else {
#ifdef HAVE_CXXABI_H
            std::string address_suffix;
            caller = remove_trailing_address_suffix(caller, address_suffix);
            {
                int demangle_status;
                const std::unique_ptr<char, free_delete<char>> freeme(
                        abi::__cxa_demangle(caller.c_str(), nullptr, nullptr, &demangle_status));
                if (freeme)
                    caller = freeme.get();
            }

            /* Get rid of the type signature, it rather useless */
            yy = get_parenthesized_arg(caller, xx, zz);
            caller = std::move(xx);

            /* could it be that we have the return type in the name
             * as well ? */

            caller+=address_suffix;
        }
#endif
        message += caller + "\n";
    }
#endif
}
