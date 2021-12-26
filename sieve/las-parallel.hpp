#ifndef LAS_PARALLEL_HPP_
#define LAS_PARALLEL_HPP_
// IWYU pragma: no_include <bits/exception.h>
#include "cado_config.h"     // for HAVE_HWLOC
#include <exception> // IWYU pragma: keep  // for exception
#include <memory>            // for shared_ptr
#include <string>            // for string
#include <sstream>           // IWYU pragma: keep
#include "macros.h"          // for GNUC_VERSION_ATMOST, MAYBE_UNUSED
#include "params.h"           // for param_list_decl_usage, cxx_param_list
#ifdef HAVE_HWLOC
#include "hwloc-aux.h"  // cxx_hwloc_nodeset
#endif

class las_parallel_desc {
    std::string description_string;
    double jobram = -1;


    int nmemory_binding_zones = 1;
    int ncpu_binding_zones_per_memory_binding_zone = 1;
    int nsubjobs_per_cpu_binding_zone = 1;
    int nthreads_per_subjob = 1;

#ifdef HAVE_HWLOC
    int memory_binding_size = 0;
    int cpu_binding_size = 0;
#endif

    struct helper; // IWYU pragma: keep
    friend struct helper;
    std::shared_ptr<helper> help;
public:
#if 0
    /* The "loose" cpu binding, is often used to process
     * or prepare data that is of interest to all cores, and we wish to
     * do that collectively. We wonder, however, how many threads this
     * means. The following rules are followed:
     *
     *  - no hwloc: same integer as the number of threads of the unique
     *  subjob.
     *  - without job replication: total number of PUs in the binding context.
     *  - with job replication: total number of PUs on the machine.
     */
    int number_of_threads_loose() const;
#else
    /* changed my mind. I don't see why number_of_threads_loose() and
     * number_of_threads_total() would be different things.
     */
    inline int number_of_threads_loose() const { return number_of_threads_total();
    }
#endif

    int number_of_memory_binding_zones() const { return nmemory_binding_zones; }
    int number_of_subjobs_per_cpu_binding_zone() const { return nsubjobs_per_cpu_binding_zone; }
    int number_of_threads_per_subjob() const { return nthreads_per_subjob; }

    int number_of_subjobs_per_memory_binding_zone() const {
        return nsubjobs_per_cpu_binding_zone *
            ncpu_binding_zones_per_memory_binding_zone;
    }
    int number_of_subjobs_total() const {
        return nmemory_binding_zones *
            number_of_subjobs_per_memory_binding_zone();
    }
    int number_of_threads_total() const {
        return number_of_subjobs_total() * nthreads_per_subjob;
    }
    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "t", "Number of threads and subjobs. Use -t help for extended documentation");
#ifdef HAVE_HWLOC
        param_list_decl_usage(pl, "job-memory", "Estimated memory per subjobs, used for job placement (see -t help) (in gigabytes, floating point values allowed)");
        param_list_decl_usage(pl, "memory-margin", "Keep this amount of RAM free for the rest of the world (see -t help) (in gigabytes, floating point values allowed)");
#else
        param_list_decl_usage(pl, "job-memory", "(unused, needs hwloc)\n");
#endif
    }
    las_parallel_desc();
    las_parallel_desc(cxx_param_list & pl, double jobram = -1);
    void display_binding_info() const;
    int set_loose_binding() const;
    int set_subjob_binding(int k) const;
    int set_subjob_cpu_binding(int k) const;
    int set_subjob_mem_binding(int k) const;

#ifdef HAVE_HWLOC
    cxx_hwloc_nodeset current_memory_binding() const;
#endif

    struct needs_job_ram : public std::exception {
        virtual const char * what() const noexcept {
            return "The \"fit\" specifier requires an estimate of the available RAM\n";
        }
    };

    class bad_specification : public std::exception {
        std::string message;
#if ! GNUC_VERSION_ATMOST(4,7,99)
        void build_message(std::ostream &) { }
        template <typename Car, typename... Cdr>
        void build_message(std::ostream & os, Car car, Cdr... cdr) {
            build_message(os << car, cdr...);
        }
#endif
        public:
        template<typename ...Args>
        bad_specification(Args&&... args MAYBE_UNUSED) {
#if GNUC_VERSION_ATMOST(4,7,99)
            message = "no message for exception, missing compiler feature";
            /* ostream:600 in gcc-4.7.2 has this code:
             *
    template<typename _CharT, typename _Traits, typename _Tp>
    inline basic_ostream<_CharT, _Traits>&
    operator<<(basic_ostream<_CharT, _Traits>&& __os, const _Tp& __x)
    { return (__os << __x); }
             *
             * which has been replaced by a much more subtle
             * implementation in recent g++. Seems to guard against
             * compilation ambiguities: in particular, the use of
             * variadic templates triggers difficulties.
             */
#else
            std::ostringstream os;
            build_message(os, args...);
            message = os.str();
#endif
        }
        virtual const char * what() const noexcept { return message.c_str(); }
    };
};

#endif	/* LAS_PARALLEL_HPP_ */
