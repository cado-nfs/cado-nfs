#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <hwloc/bitmap.h>
// IWYU pragma: no_include "hwloc/bitmap.h"

#include <cerrno>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

// #include <mutex>
#include <string>
#include <tuple>
#include <vector>
#include <memory>
#include <sstream>

#include <strings.h>
#include <regex.h>

#ifdef HAVE_HWLOC
#include <hwloc.h>
#include "hwloc-aux.h"
#endif
#include "las-parallel.hpp"

#include "misc.h"
#include "utils_cxx.hpp"
#include "verbose.h"
#include "macros.h"
#include "params.h"

static const char * const default_placement_with_auto = "node,fit*4,fit,pu,loose";

static bool parse_number(std::string const & s, int & x, std::string::size_type pos = 0) /*{{{*/
{
    const char * digits = "0123456789";
    if (s.empty() || s.find_first_not_of(digits, pos) != std::string::npos)
        return false;
    std::istringstream is(s.substr(pos));
    return static_cast<bool>(is >> x);
}/*}}}*/

/* used to help with some of the job achieved by the las_parallel ctor.
 * Otherwise no reason to expose within the public class, as it's of no
 * real use (not even expository).
 */
struct las_parallel_desc::helper {
    std::string memory_binding_specifier_string;
    std::string cpu_binding_specifier_string;
    std::string jobs_within_cpu_binding_string;
    std::string threads_per_job_string;

    /* The banner is always printed once at least. The full description
     * only with -v */
    std::string banner;
    std::string full_diagnostics;

#ifdef HAVE_HWLOC
    hwloc_topology_t topology;

    std::string synthetic_topology_string;
    int depth = 0;
    std::vector<int> depth_per_level;

    /* for the "fit" keyword */
    mutable int computed_min_pu_fit = -1;

    /* size() == number of memory binding zones */
    std::vector<cxx_hwloc_nodeset> memory_binding_nodesets;

    /* size() == total number of subjobs */
    std::vector<cxx_hwloc_cpuset> subjob_binding_cpusets;

    double total_ram_margin = 0;

    int memory_binding_size;
    int cpu_binding_size;
#endif
    int nsubjobs_per_cpu_binding_zone = 1;
    int nthreads_per_subjob = 1;

    bool loose = false;
    bool replicate = true;

    bool is_loose() const { return loose; }
    bool is_strict() const { return !loose; }
    bool want_replicate() const { return replicate; }
    /*
    bool has_option(std::string & const opt) {
        return binding_options.find(opt) != binding_options.end();
    }
    */

    void set_banner(std::string const & description_string) {
        std::ostringstream os;
#ifdef HAVE_HWLOC
        os  << "# Applying binding "
            << description_string
            << " on a machine with topology "
            << synthetic_topology_string
            << " (" << (total_ram() >> 30) << " GB RAM)\n",
#else
            os  << "# Applying binding " << description_string << " with hwloc disabled\n",
#endif
            banner += os.str();
    }
#ifdef HAVE_HWLOC
    int mod(int&i) const {/*{{{*/
        i = i % depth; if (i < 0) i += depth;
        return i;
    }/*}}}*/
    int number_at(int i) const {/*{{{*/
        return hwloc_get_nbobjs_by_depth(topology, mod(i));
    }/*}}}*/
    int number_of(int child_depth, int parent_depth) const/*{{{*/
    {
        mod(child_depth);
        mod(parent_depth);
        if (child_depth < parent_depth) return 0;
        int n = 1;
        for(int i = parent_depth ; i < child_depth ; ++i)
            n *= depth_per_level[i+1];
        return n;
    }/*}}}*/
#endif
    helper() {/*{{{*/
        /* Here we need hwloc, of course */
#ifdef HAVE_HWLOC
        hwloc_topology_init(&topology);
        {/*{{{ load the topology `*/
            unsigned long flags = 0;
#if HWLOC_API_VERSION >= 0x010700
            flags = hwloc_topology_get_flags(topology);
#endif  /* HWLOC_API_VERSION >= 0x010700 */
            /* we must make sure to remove these flags, but it's likely that
             * they're off by default anyway */
#if HWLOC_API_VERSION < 0x020000
            flags &= ~(HWLOC_TOPOLOGY_FLAG_IO_DEVICES | HWLOC_TOPOLOGY_FLAG_IO_BRIDGES);
#endif
            hwloc_topology_set_flags(topology, flags);

            hwloc_topology_load(topology);

        }/*}}}*/
        depth = hwloc_topology_get_depth(topology);
        depth_per_level.reserve(depth);
        int n = 1;
        std::ostringstream const os;
        for(int i = 0 ; i < depth ; i++) {
            int const x = hwloc_get_nbobjs_by_depth(topology, i);
            depth_per_level.push_back(x / n);
            n = x;
        }

        hwloc_obj_t root = hwloc_get_root_obj(topology);
        if (!root->symmetric_subtree) {
            fprintf(stderr, "# Topology is not symmetric,"
                    " cannot proceed with replication"
                    " of the las process with the current code."
                    " No cpu/memory binding will be set.\n");
            /* simply stick to the default. We'll use depth=0 as an
             * indicator of this situation, and we must make sure
             * that we return sane defaults in all cases. It's quite
             * ugly and should be reworked.
             */
            return;
        }

        /* Form a sensible set of defaults, for one unbound thread */
        memory_binding_size = number_of(-1,0);
        cpu_binding_size = memory_binding_size;
        compute_binding_bitmaps();

        char buf[1024];
        hwloc_topology_export_synthetic(topology, buf, sizeof(buf), HWLOC_TOPOLOGY_EXPORT_SYNTHETIC_FLAG_NO_ATTRS );
        synthetic_topology_string = buf;

#endif
    }/*}}}*/

    ~helper() {/*{{{*/
#ifdef HAVE_HWLOC
        hwloc_topology_destroy(topology);
#endif
    }/*}}}*/
    std::vector<std::string> tokenize(std::string const & s) const {/*{{{*/
        std::vector<std::string> tokens;
        for(std::string::size_type x = 0, y; x != std::string::npos ; x = y) {
            y = s.find(',',x);
            if (y == std::string::npos) {
                tokens.push_back(s.substr(x));
            } else {
                tokens.push_back(s.substr(x, y-x));
                y++;
            }
        }
        return tokens;
    }/*}}}*/
    void replace_aliases(std::string & desc) const {/*{{{*/
        int k = 0;
        if (parse_number(desc, k)) {
            std::ostringstream os;
            os << "machine,1," << k; 
            desc = os.str();
            return;
        }
        if (desc == "single") { desc = "machine,1,pu"; return; }
        if (desc == "auto") { desc = default_placement_with_auto; return; }
        if (desc == "auto,no-replicate") { desc = default_placement_with_auto; desc += ",no-replicate"; return; }
        if (desc.substr(0,7) == "single-") {
            std::ostringstream os;
            os << desc.substr(7) << ",1,pu";
            desc = os.str();
            return;
        }
    }/*}}}*/
    void parse(std::string & desc) { /* {{{ */
        replace_aliases(desc);

        std::vector<std::string> tokens = tokenize(desc);

        bool loose_opt = false;
        bool strict_opt = false;
        bool norepl_opt = false;
        for( ; tokens.size() >= 3 ; ) {
            std::string const & s(tokens.back());
            if (s == "strict") { strict_opt = true; }
            else if (s == "loose") { loose_opt = true; }
            else if (s == "no-replicate") { norepl_opt = true; }
            else break;
            tokens.erase(--tokens.end());
        }
        if (loose_opt && strict_opt)
            throw bad_specification("loose and strict are incompatible");
        if (loose_opt) loose = true;
        if (strict_opt) loose = false;
        if (norepl_opt) replicate = false;

        if (tokens.size() == 3) {
            memory_binding_specifier_string = tokens[0];
            jobs_within_cpu_binding_string = tokens[1];
            threads_per_job_string = tokens[2];
        } else if (tokens.size() == 4) {
            memory_binding_specifier_string = tokens[0];
            cpu_binding_specifier_string = tokens[1];
            jobs_within_cpu_binding_string = tokens[2];
            threads_per_job_string = tokens[3];
        } else {
            throw bad_specification("not the right number of tokens, or wrong options");
        }
    }/*}}}*/
#ifdef HAVE_HWLOC
    uint64_t total_ram() const {/*{{{*/
        hwloc_obj_t root = hwloc_get_root_obj(topology);
        /* I'm not sure about how I want to get the info on the amount of
         * available ram. hwloc-info, as such, does not seem to document
         * a way (beyond the fact that this bit of info is output with
         * "machine" and "Numanode" specifiers). There is an
         * hwloc_obj_memory_s type in the hwloc C api, and that seems to
         * be attached to hwloc_obj structures as well. But at any rate,
         * this gives some notion of available memory. This is well and
         * good, except that it might be misleading if a job that takes
         * "almost" 4GB is scheduled to be placedonly 4 times on a 16GB
         * machine when we specify --job-memory 4 (with some confidence
         * that this is a convenient enough upper bound).
         */

#if HWLOC_API_VERSION < 0x020000
        uint64_t ram = root->memory.total_memory;
#else
        uint64_t ram = root->total_memory;
#endif
        /* Round this up to at most 1/16-th. This will do rubbish on a
         * machine where the 1 bits in the binary expansion of the
         * hardware ram size spread more than 4 positions, but we find
         * that unlikely.
         */
        for(uint64_t x = ram >> 4; x ; x >>= 1) ram |= x;
        ram = ram + 1;
        ram = ram - ((uint64_t) (total_ram_margin * (1<<30)));
        return ram;
    }/*}}}*/
   int enclosing_depth(int n) const {/*{{{*/
       /* return deepest level k such that an object at depth k contains
        * more than, or exactly N PUs.
        */
       int k;
       int m = 1;
       ASSERT_ALWAYS(n >= 1);
       for(k = depth - 1; n > m && k >= 0 ; k--) {
           ASSERT(number_of(-1, k) == m);
           m *= depth_per_level[k];
           ASSERT(k == 0 || number_of(-1, k-1) == m);
       }

       if (k >= 0) {
           ASSERT(number_of(-1, k) >= n);
           if (k == depth - 1) {
               ASSERT(n == 1);
           } else {
               ASSERT(number_of(-1, k + 1) < n);
           }
       }

       return k;
   }/*}}}*/
   std::tuple<int, int, int> flat_to_hierarchical(int n) const {/*{{{*/
       int const k = enclosing_depth(n);
       /* a binding scope that has to meet the constraint of containing
        * at least n PUs can do so with an object at depth
        * enclosing_depth(n). However, a finer grain might be possible,
        * and that may be an important optimization when we have many
        * edges in the topology tree at this level: it is perhaps
        * possible to divide the object at depth k into equal parts made
        * of several objects at depth k+1, all parts meeting the
        * constraint. We count how large these parts must be (counted in
        * terms of number of objects at depth k+1, called "children"
        * here.).
        */
       if (k == depth - 1) {
           ASSERT(n == 1);
           return std::make_tuple(k, 1, 1);
       }
       int child_depth = k + 1;
       int child_size = number_of(-1, k + 1);
       int part_size = iceildiv(n, child_size);
       int const v = number_of(k+1, k);
       ASSERT(part_size <= v);
       if (replicate)
           for( ; v % part_size ; part_size++);
       if (part_size == v) {
           child_depth--;
           child_size *= v;
           part_size = 1;
       }
       return std::make_tuple(child_depth, child_size, part_size);
   }/*}}}*/
   int acceptable_binding(int n) const {/*{{{*/
       auto x = flat_to_hierarchical(n);
       return std::get<1>(x) * std::get<2>(x);
   }/*}}}*/
   std::string textual_description_for_binding(int n) const {/*{{{*/
       int k, child_size, part_size;
       ASSERT(n <= number_of(-1, 0));
       std::tie(k, child_size, part_size) = flat_to_hierarchical(n);
       for( ; k < depth - 1 && number_of(k+1, k) == 1 ; k++);
       hwloc_obj_t obj = hwloc_get_obj_by_depth(topology, k, 0);
       char s[256];
       hwloc_obj_type_snprintf(s, sizeof(s), obj, 1);
       if (part_size == 1) {
           return s;
       } else {
           std::ostringstream os;
           os << std::string(s) << ":0-" << part_size - 1;
           return os.str();
       }
   }/*}}}*/
   std::vector<std::string> all_textual_descriptions_for_binding(int n) const {/*{{{*/
       /* n is a binding granularity. Return all the subsets of the
        * machine, all identical, that have some hardware relevance (=
        * correspond to an integer fraction of the subtree at some depth,
        * and made of entire subtrees at the level below), and all
        * contain at least n PUs.
        */
       int k, child_size, part_size;
       ASSERT(n <= number_of(-1, 0));
       std::tie(k, child_size, part_size) = flat_to_hierarchical(n);
       for( ; k < depth - 1 && number_of(k+1, k) == 1 ; k++);
       std::vector<std::string> res;
       int const nk = number_of(k, 0);
       for(int i = 0 ; i < nk ; i+=part_size) {
           char s[256];
           hwloc_obj_t obj = hwloc_get_obj_by_depth(topology, k, i);
           hwloc_obj_type_snprintf(s, sizeof(s), obj, 0);
           std::ostringstream os;
           int const a = obj->logical_index;
           int const b = a + part_size - 1;
           if (part_size == 1) {
               os << s << ":" << a;
           } else {
               os << s << ":" << a << '-' << b;
           }
           res.push_back(os.str());
       }
       return res;
   }/*}}}*/
    std::vector<cxx_hwloc_cpuset> all_cpu_bitmaps(int width, int multiply = 1) const /*{{{*/
    {
        int const n = width;
        int k, child_size, part_size;
        ASSERT(n <= number_of(-1, 0));
        std::tie(k, child_size, part_size) = flat_to_hierarchical(n);
        std::vector<cxx_hwloc_cpuset> res;
        int const nk = number_of(k, 0);
        for(int i = 0 ; i < nk ; i+=part_size) {
            cxx_hwloc_cpuset c;
            for(int j = 0 ; j < part_size ; j++)
                c = c | hwloc_get_obj_by_depth(topology, k, i + j)->cpuset;
            for(int k = 0 ; k < multiply ; k++)
                res.push_back(c);
        }
        return res;
    }/*}}}*/
    std::vector<cxx_hwloc_nodeset> all_mem_bitmaps(int width, int multiply = 1) const /*{{{*/
    {
        int const n = width;
        int k, child_size, part_size;
        ASSERT(n <= number_of(-1, 0));
        std::tie(k, child_size, part_size) = flat_to_hierarchical(n);
        std::vector<cxx_hwloc_cpuset> res;
        int const nk = number_of(k, 0);
        for(int i = 0 ; i < nk ; i+=part_size) {
            cxx_hwloc_nodeset c;
            for(int j = 0 ; j < part_size ; j++)
                c = c | hwloc_get_obj_by_depth(topology, k, i + j)->nodeset;
            for(int k = 0 ; k < multiply ; k++)
                res.push_back(c);
        }
        return res;
    }/*}}}*/
   int min_pu_fit(double jobram) {/*{{{*/
       if (computed_min_pu_fit >= 0) return computed_min_pu_fit;
       /* This computes the minimal number n of PUs such that
        * floor(total_PUs/n) jobs, each using jobram GB, fit on the
        * machine and do not exceed the total amount of __physical__ RAM.
        *
        * hwloc does give us a useful bit of information regarding
        * available, and not physical ram. We strive to get the physical
        * ram because that amount is generally better known to the user,
        * and one expects the automatic placement logic to find that
        * *four* jobs of say 12GB would fit on a 48GB machine.
        *
        * Note that this integer is not yet of interest to the
        * hierarchical topology.
        */
       /* how many jobs can fit ? */
       int const how_many = total_ram() / (jobram * (1<<30));
       if (how_many == 0) {
           fprintf(stderr, "This machine does not have enough memory"
                   " (%" PRIu64 " GB) to fit jobs that need %.2f GB\n",
                   total_ram() >> 30, jobram);
           exit(EXIT_FAILURE);
       }
       /* This means subjobs can't be smaller than this: */
       computed_min_pu_fit = iceildiv(number_of(-1, 0), how_many);
       std::string const s = textual_description_for_binding(computed_min_pu_fit);
       std::ostringstream os;
       os  << "# Given our estimate of "
           << size_disp(jobram * (1<<30)) << " RAM per subjob,"
           << " the \"fit\" binding\n"
           << "# corresponds to a minimum size of "
           << computed_min_pu_fit << " covered PUs in the hwloc hierarchy,\n"
           << "# i.e. no less than \"" << s << "\"\n";
       full_diagnostics += os.str();
       return computed_min_pu_fit;
   }/*}}}*/
    void compute_binding_bitmaps() {/*{{{*/
        subjob_binding_cpusets = all_cpu_bitmaps(cpu_binding_size, nsubjobs_per_cpu_binding_zone);
        memory_binding_nodesets = all_mem_bitmaps(memory_binding_size);
    }/*}}}*/
    cxx_hwloc_nodeset current_memory_binding() const {/*{{{*/
        cxx_hwloc_nodeset nn;
        hwloc_membind_policy_t pol;
#if HWLOC_API_VERSION < 0x010b03
        /* this legacy called remained valid throughout hwloc 1.x */
        int rc = hwloc_get_membind_nodeset(topology,  nn, &pol,
                HWLOC_MEMBIND_THREAD);
#else
        /* newer call */
        int const rc = hwloc_get_membind(topology,  nn, &pol,
                HWLOC_MEMBIND_THREAD | HWLOC_MEMBIND_BYNODESET);
#endif
        if (rc < 0) {
            /*
            static std::mutex mm;
            std::lock_guard<std::mutex> const dummy(mm);
            static int got_message = 0;
            if (!got_message++)
                fprintf(stderr, "Error while attempting to get memory binding\n");
                */
            hwloc_bitmap_zero(nn);
        }
        return nn;
    }/*}}}*/
#endif
   int interpret_memory_binding_specifier() {
#ifdef HAVE_HWLOC
       if (depth) {
           memory_binding_size = interpret_generic_binding_specifier(memory_binding_specifier_string);
           std::ostringstream os;
           os
               << "# Memory binding specifier: request: "
               << memory_binding_specifier_string
               << " ; resolves to: "
               << textual_description_for_binding(memory_binding_size)
               << "\n";
           full_diagnostics += os.str();
           return memory_binding_size;
       } else
#endif
       {
           /* just check for errors */
           return interpret_generic_binding_specifier(memory_binding_specifier_string);
       }
   }
   int interpret_cpu_binding_specifier() {
#ifdef HAVE_HWLOC
       if (depth) {
           if (cpu_binding_specifier_string.empty()) {
               cpu_binding_size = memory_binding_size;
           } else {
               /* cap to memory binding size */
               cpu_binding_size = interpret_generic_binding_specifier(cpu_binding_specifier_string, memory_binding_size);
               if (memory_binding_size % cpu_binding_size)
                   throw bad_specification(
                           "\n",
                           banner,
                           full_diagnostics,
                           "# ERROR: cpu binding "
                           , cpu_binding_specifier_string,
                           " yields a size of ", cpu_binding_size, " PUs,"
                           " but this must be an integer divisor of",
                           " the memory binding size of ", memory_binding_size,
                           " PUs, which follows from the specifier ",
                           memory_binding_specifier_string);
           }
           std::ostringstream os;
           os
               << "# Cpu binding specifier: request: "
               << cpu_binding_specifier_string
               << " ; resolves to: "
               << textual_description_for_binding(cpu_binding_size)
               << "\n";
           full_diagnostics += os.str();
           return cpu_binding_size;
       } else
#endif
       {
           /* just check for errors */
           if (cpu_binding_specifier_string.empty()) {
               return 0;
           } else {
               return interpret_generic_binding_specifier(cpu_binding_specifier_string);
           }
       }
   }
   int interpret_generic_binding_specifier(std::string const & specifier, int cap MAYBE_UNUSED = -1) const {/*{{{*/
#ifdef HAVE_HWLOC
       if (!depth) {
           if (strcasecmp(specifier.c_str(), "machine") != 0)
               throw bad_specification("hwloc detected asymmetric topology, the only accepted memory binding specifier is \"machine\"");
           return 0;
       }
       int binding_size = 0;
       /* and apply our different calculation rules to the provided
        * string (either memory_binding_specifier_string or
        * cpu_binding_specifier_string if there is one). */
       if (parse_number(specifier, binding_size))
           return binding_size;
       const char * binding_specifier_regexp = 
           "^("                         // 1 : non-limited part.
               "([^*/[:space:]]*)"      // 2: object, or "fit"
               "(\\*([0-9]+))?"          // 3,4: coarsening
               "(/([0-9]+))?"           // 5,6: restricting
           ")"
           "(/([^[:space:]]*))?"        // 7,8: limiting
           "$";
       regex_t R;
       regcomp(&R, binding_specifier_regexp, REG_ICASE|REG_EXTENDED);
       auto dummy = call_dtor([&](){regfree(&R);});
       regmatch_t m[9];
       int const r = regexec(&R, specifier.c_str(), 9, m, 0);
       if (r != 0)
           throw bad_specification("binding specifier ",
                   specifier, " is invalid");
       int objsize;
       auto sub = [&](int i, bool want = false) {
           std::string res;
           if (want) ASSERT_ALWAYS(m[i].rm_so >= 0);
           if (m[i].rm_so >= 0)
               res = specifier.substr(m[i].rm_so, m[i].rm_eo - m[i].rm_so);
           return res;
       };
       std::string const base_object = sub(2, true);
       bool const is_fit = strcasecmp(base_object.c_str(), "fit") == 0;
       if (is_fit) {
           if (computed_min_pu_fit < 0) throw needs_job_ram();
           objsize = computed_min_pu_fit;
       } else {
           int argdepth = hwloc_aux_get_depth_from_string(topology, base_object.c_str());
#if HWLOC_API_VERSION >= 0x020000
           if (argdepth == HWLOC_TYPE_DEPTH_NUMANODE) {
               argdepth = hwloc_get_memory_parents_depth(topology);
           }
#endif
           if (argdepth < 0)
               throw bad_specification(base_object, " is invalid");
           objsize = number_of(-1, argdepth);
       }
       std::string const multiplier_string = sub(4);
       if (!multiplier_string.empty()) {
           int x = 0;
           bool const t = parse_number(multiplier_string, x);
           ASSERT_ALWAYS(t);
           objsize *= x;
       }
       if (cap < 0) cap = number_of(-1, 0);
       objsize = std::min(objsize, cap);
       if (is_fit)
           objsize = acceptable_binding(objsize);
       std::string const divisor_string = sub(6);
       if (!divisor_string.empty()) {
           int x = 0;
           bool const t = parse_number(divisor_string, x);
           ASSERT_ALWAYS(t);
           if (!x || (objsize % x)) {
               std::ostringstream os;
               os << "specifier binding specifier is invalid. Cannot divide "
                   << base_object << sub(3) << " (=" << objsize << ") into "
                   << x << " parts";
               if (is_strict()) {
                   throw bad_specification(os.str());
               } else {
                   fprintf(stderr, "ERROR: %s\n", os.str().c_str());
               }
           }
           if (x) objsize /= x;
       }
       std::string const limiting_string = sub(8);
       if (!limiting_string.empty()) {
           int const argdepth = hwloc_aux_get_depth_from_string(topology, limiting_string.c_str());
           if (argdepth < 0)
               throw bad_specification(limiting_string, " is invalid");
           int const compare = number_of(-1, argdepth);
           if (objsize > compare) {
               std::ostringstream os;
               os << "automated binding "
                   << textual_description_for_binding(objsize)
                   << " deduced from \"" << sub(2) << sub(4) << sub(6) << "\""
                   << " (= " << objsize << ")"
                   << " escapes the advised binding given by " << sub(7)
                   << " (which is "
                   << textual_description_for_binding(compare) 
                   << ")";
               if (is_strict()) {
                   throw bad_specification(os.str());
               } else {
                   fprintf(stderr, "Warning: %s\n", os.str().c_str());
               }
           }
       }
       return objsize;
#else
       if (strcasecmp(specifier.c_str(), "machine") != 0)
           throw bad_specification("hwloc being disabled, the only accepted binding specifier is \"machine\"");
       return 0;
#endif
   }/*}}}*/
   int interpret_njobs_specifier() {/*{{{*/
       if (parse_number(jobs_within_cpu_binding_string, nsubjobs_per_cpu_binding_zone))
           return nsubjobs_per_cpu_binding_zone;
#ifdef HAVE_HWLOC
       int objsize;
       if (strcasecmp(jobs_within_cpu_binding_string.c_str(), "fit") == 0) {
           if (computed_min_pu_fit < 0) throw needs_job_ram();
           objsize = acceptable_binding(computed_min_pu_fit);
       } else if (jobs_within_cpu_binding_string.find_first_of("*/") != std::string::npos) {
           throw bad_specification("Only binding specifiers allow multipliers");
       } else {
           int const argdepth = hwloc_aux_get_depth_from_string(topology, jobs_within_cpu_binding_string.c_str());
           if (argdepth < 0)
               throw bad_specification(jobs_within_cpu_binding_string, " is invalid");
           objsize = number_of(-1, argdepth);
       }
       if (cpu_binding_size % objsize)
           throw bad_specification(
                   "\n",
                   banner,
                   full_diagnostics,
                   "ERROR: cannot place jobs according to",
                   " the ", jobs_within_cpu_binding_string, " rule",
                   " (which resolves to: ",
                   textual_description_for_binding(objsize), "),"
                   " as there is not an integer number of these",
                   " in the imposed binding ",
                   textual_description_for_binding(cpu_binding_size),
                   " [HINT: are you using -t auto and running out of RAM on NUMA nodes?]");
       nsubjobs_per_cpu_binding_zone = cpu_binding_size / objsize;
       std::ostringstream os;
       os
           << "# Njobs specifier: request: "
           << jobs_within_cpu_binding_string
           << " ; resolves to: "
           << textual_description_for_binding(objsize)
           << " (" << nsubjobs_per_cpu_binding_zone << " jobs)"
           << "\n";
       full_diagnostics += os.str();
       return nsubjobs_per_cpu_binding_zone;
#else
       throw bad_specification("hwloc being disabled, the only accepted specifiers for the number of jobs are integers");
#endif
   }/*}}}*/
   int interpret_nthreads_specifier() {/*{{{*/
       if (parse_number(threads_per_job_string, nthreads_per_subjob))
           return nthreads_per_subjob;
#ifdef HAVE_HWLOC
       int objsize;
       int const argdepth = hwloc_aux_get_depth_from_string(topology, threads_per_job_string.c_str());
       if (argdepth < 0)
           throw bad_specification(threads_per_job_string, " is invalid");
       objsize = number_of(-1, argdepth);
       int const loose_per_job_scale = cpu_binding_size / nsubjobs_per_cpu_binding_zone;
       if (loose_per_job_scale % objsize)
           throw bad_specification("cannot place threads according to",
                   " the ", threads_per_job_string, " rule,",
                   " as there is not an integer number of these",
                   " in the per-job fraction (1/", nsubjobs_per_cpu_binding_zone, ")",
                   " of the binding ",
                   textual_description_for_binding(cpu_binding_size));
       nthreads_per_subjob = loose_per_job_scale / objsize;
       std::ostringstream os;
       os
           << "# Nthreads specifier: request: "
           << threads_per_job_string
           << " ; resolves to: "
           << textual_description_for_binding(objsize)
           << " (" << nthreads_per_subjob << " threads)"
           << "\n";
       full_diagnostics += os.str();
       return nthreads_per_subjob;
#else
       throw bad_specification("hwloc being disabled, the only accepted specifiers for the number of threads are integers");
#endif
   }/*}}}*/
};

static void extended_usage()/*{{{*/
{
    std::ostringstream os;
    os << R"(
Documentation for the --job-binding-policy (-t) option.
=======================================================

Option takes either an alias, or a longer specificier.  Most of the
functionality depends on the hwloc library being usable. Examples are
provided after the syntax description


Aliases, first:
===============

    [an integer N]: equivalent to Machine,1,N
        (run one N-threaded job, not bound to anything).
    single: equivalent to single-machine
    single-XXX: equivalent to XXX,1,PU
)"
    << "    auto: equivalent to " << default_placement_with_auto << " [see below]\n"
    << "    auto,no-replicate: equivalent to " << default_placement_with_auto << ",no-replicate [see below]\n"
    << R"(
Syntax of the specifiers
========================

[memory binding level],[cpu binding level],[number of jobs with this binding],[threads per job][,modifiers]

The cpu binding level may be omitted, and defaults to the same value as
the machine binding level.

All three (or four) main items may be set as integers. However this is
probably not the most useful, and certainly not the most portable way to
set them.

Binding level:
--------------

Both the memory binding level and the cpu binding level follow the exact
same syntax.

First observe the output ow hwloc-info (may vary depending on the
machine).
        localhost ~ $ hwloc-info --no-io
        depth 0:	1 Machine (type #1)
         depth 1:	2 NUMANode (type #2)
          depth 2:	2 Package (type #3)
           depth 3:	2 L3Cache (type #4)
            depth 4:	16 L2Cache (type #4)
             depth 5:	16 L1Cache (type #4)
              depth 6:	16 Core (type #5)
               depth 7:	32 PU (type #6)
This means 2 NUMANode per machine, 8 Core per NUMANode, 2 PU per
Core.

The syntax of the binding level is
    [integer]
    OR [hwloc object name](\*[integer]|/[integer])?
    OR fit(\*[integer])?(/[hwloc object name])?

Here's how it goes to specify the binding level as an integer. As
mentioned, this is neither easy nor portable, look further for the
recommended way to set via level aliases.  An integer N defines a binding
level which corresponds to N of the innermost instances (here, PU is the
innermost object.  It is probably always so. The text below says PUs, but
the code means the innermost object anyway). Let k be the depth
such that an object at depth k contains strictly more than N PUs, while
an object at depth k+1 contains less than, or exactly N PUs. For example,
above, for N=8 we have k=3 (at depth k+1 a single L2Cache is above less
than 8 PUs, but at level k=3, below one L3Cache object we have 16 PUs,
which is more than 8). When specifying an integer N, it must be a
multiple of the number of PUs below an object at level k+1 (here, 8 is
indeed a multiple of 2, which is the number of PUs below an L2Cache), as
well as a divisor of the number of PUs below an object at level k (here,
8 divides 16). The latter requirement is waived if the "no-replicate"
modifier is present.

We bind jobs to hwloc objects at a particular depth by writing simply the
object name (e.g. NUMANode) as a binding level. Names are never
case-sensitive, and match in the same way hwloc tools such as hwloc-calc
recognizes them. In the example above, this is equivalent to the integer
16, but this varies. When specifying an object name, it is possible to
append "*N" or "/N" to loosen or restrict the binding level, provided the
resulting integer satisfies the constraints above.

The "fit" keyword is special. It corresponds to the minimum number of PUs
that must be grouped together so that if the machine is filled with
identical jobs, the memory available on the machine suffices. This relies
on an estimation of the amount of required memory, of course (see also
--job-memory). It is possible to write fit*N, which loosens this binding
(always subject to the same requirement). When appending /XXX, where XXX
denotes an object name, this ensures that the corresponding binding does
not escape the size of one object of that name. See also the "loose" and
"strict" modifiers.

Number of jobs with same binding
--------------------------------

This is an integer, too. We may also use hwloc specifiers as shortcuts.
However the interpretation is reversed compared to the binding level.
Here, "Core" indicates a number which is the number of "Core" objects
within the binding specification.

The "fit" keyword echoes the functionality of the "fit" binding
specifier. Again, it means that within the binding specifier, we put
exactly the number of jobs so that if this is replicated over the whole
machine, the memory fits.

Naturally, when using an object name (or "fit"), it is mandatory that the
binding specifier contains a positive integer number of objects of that
name.

Number of threads per job
-------------------------

Again, an integer. As with the number of jobs, aliases are understood
relative to the job size (which may occupy only a fraction of the binding
level, but always an integer fraction anyway).

The two most useful settings are "PU" and "Core", which are equivalent to
using or not using hyperthreading.

Modifiers
---------

'loose' / 'strict' (default: 'loose'). If a binding specifier of the form
fit*N/XXX or fit/XXX is used, warn (with 'loose') or error out (with
'strict') if the calculated binding is coarser than XXX.

'no-replicate' : by default, we may put several jobs with the same
binding, and the set of jobs that are put in a single binding is
replicated over the whole machine (as many times as this binding fits in
the machine). With the "no-replicate" option, this replication is
disabled (however we still have the same number of sub-jobs with one
binding).

Examples
========

    -t 8    -> resolves to machine,1,8 = run only one job, with no
            particular binding. There is only one job, because there is
            only one "machine" per... machine, so whether or not we use
            the "no-replicate" modifier makes no difference.

    -t numa,fit,pu,strict
            -> binds jobs at the numa level. Set as many jobs as can fit,
            and then run each on exactly the number of pus (hyperthreaded
            cores). Abort if it is not possible to fit on one numa node.

    -t numa,core*2,fit,pu
            -> memory-bind at the numa level, and do cpu-bind on pairs of
            cores. Stick as many different jobs as can fit in each cpu
            binding context, and have each of these jobs use as many
            threads as the number of pus we have (per sub-job in the cpu
            binding context of two cores).

    -t numa,fit,20
            -> same, but force 20-thread jobs (no matter what).

    -t auto -> resolves to )" << default_placement_with_auto << R"( = bind fractions of NUMA
            nodes, and put as many jobs in there as can fit. Then have
            each run the number of threads so that all PUs are busy.

    -t package/2,core,pu
            -> binds to half packages, run 1 job per core, and have each
            job use as many threads as we have PUs per core.
)";
    fputs(os.str().c_str(), stderr);
}/*}}}*/

las_parallel_desc::las_parallel_desc()
    : help (std::make_shared<helper>())
{}

las_parallel_desc::las_parallel_desc(cxx_param_list & pl, double jobram_arg)
    : las_parallel_desc()
{
    jobram = jobram_arg;
    const char * desc_c0 = param_list_lookup_string(pl, "t");

    if (!desc_c0) {
        verbose_output_start_batch();
        verbose_output_print(0, 1, "# No -t option found, running single-threaded.");
#ifdef HAVE_HWLOC
        verbose_output_print(0, 1, " See also -t help");
#endif
        verbose_output_print(0, 1, "\n");
        verbose_output_end_batch();
        /* feed a constant string */
        description_string = "1";
    } else {
        description_string = desc_c0;
    }

    if (description_string == "help") {/*{{{*/
        extended_usage();
        exit(EXIT_SUCCESS);
    }/*}}}*/

    help->replace_aliases(description_string);
    help->set_banner(description_string);

#ifdef HAVE_HWLOC
    param_list_parse_double(pl, "memory-margin", &help->total_ram_margin);
    /* The --job-memory argument will **ALWAYS** win */
    param_list_parse_double(pl, "job-memory", &jobram);
    if (jobram != -1)
        help->min_pu_fit(jobram);
#endif

    try {
        help->parse(description_string);
    } catch (bad_specification & b) {
        throw bad_specification("Cannot interpret specification -t ",
                description_string,
                ": ", b.what(), "."
                " See -t help for extended documentation");
    }
    help->interpret_memory_binding_specifier();
    help->interpret_cpu_binding_specifier();
    help->interpret_njobs_specifier();
    help->interpret_nthreads_specifier();

    nsubjobs_per_cpu_binding_zone = help->nsubjobs_per_cpu_binding_zone;
    nthreads_per_subjob = help->nthreads_per_subjob;

#ifdef HAVE_HWLOC
    if (!help->depth) return;
    memory_binding_size = help->memory_binding_size;
    nmemory_binding_zones = help->number_of(-1, 0) / memory_binding_size;
    if (!help->replicate)
        nmemory_binding_zones = 1;
    cpu_binding_size = help->cpu_binding_size;
    ncpu_binding_zones_per_memory_binding_zone = memory_binding_size / cpu_binding_size;

    /* Need to populate the arrays subjob_binding_cpusets and
     * memory_binding_nodesets. */
    help->compute_binding_bitmaps();

#if 0
    const struct hwloc_topology_support* s = hwloc_topology_get_support(help->topology);
    /* see /usr/share/doc/libhwloc-doc/html/a00234.html */
    /* just examples of flags we may be interested in */
    ASSERT_ALWAYS(s->membind->get_area_membind);
#endif

#endif
}

void las_parallel_desc::display_binding_info() const /*{{{*/
{
    verbose_output_start_batch();
    verbose_fmt_print(0, 1, "{}", help->banner);
    verbose_fmt_print(0, 2, "{}", help->full_diagnostics);
#ifdef HAVE_HWLOC
    verbose_fmt_print(0, 1, "# {} memory binding zones\n", nmemory_binding_zones);
    verbose_fmt_print(0, 1, "# {} cpu binding zones within each memory binding\n", ncpu_binding_zones_per_memory_binding_zone);
    verbose_fmt_print(0, 1, "# {} jobs within each binding context (hence {} in total)\n",
            number_of_subjobs_per_cpu_binding_zone(),
            number_of_subjobs_total());
    verbose_fmt_print(0, 1, "# {} threads per job\n", nthreads_per_subjob);
    if (jobram >= 0) {
        double const all = jobram * number_of_subjobs_total();
        int const physical = help->total_ram() >> 30;
        double const ratio = 100 * all / physical;

        verbose_fmt_print(0, 1, "# Based on an estimate of {:.2f} GB per job, we use {:.2f} GB in total, i.e. {:.1f}% of {} GB\n",
                jobram, all, ratio, physical);
    }
#else
    verbose_fmt_print(0, 1, "# {} jobs in parallel\n", number_of_subjobs_total());
    verbose_fmt_print(0, 1, "# {} threads per job\n", nthreads_per_subjob);
#endif


#ifdef HAVE_HWLOC
    if (help->depth) {
        std::vector<std::string> tm = help->all_textual_descriptions_for_binding(memory_binding_size);
        std::vector<std::string> tc = help->all_textual_descriptions_for_binding(cpu_binding_size);
        std::vector<std::string> const tj = help->all_textual_descriptions_for_binding(cpu_binding_size / number_of_subjobs_per_cpu_binding_zone());
        size_t m = 0, c = 0;
        {
            std::ostringstream pu_app;
            pu_app << "[" << memory_binding_size << " PUs]";
            for(auto & x : tm) {
                x += " ";
                x += pu_app.str();
                m = std::max(m, x.size());
            }
        }
        {
            std::ostringstream pu_app;
            pu_app << "[" << cpu_binding_size << " PUs]";
            for(auto & x : tc) {
                x += " ";
                x += pu_app.str();
                c = std::max(c, x.size());
            }
        }
        size_t const qc = number_of_subjobs_per_cpu_binding_zone();
        size_t const qm = qc * ncpu_binding_zones_per_memory_binding_zone;
        for(size_t i = 0 ; i < tj.size() ; i++) {
            if (!help->replicate && i >= qc) break;
            ASSERT_ALWAYS(i < help->subjob_binding_cpusets.size());
            char * sc;
            hwloc_bitmap_asprintf(&sc, help->subjob_binding_cpusets[i]);
            ASSERT_ALWAYS(i/qm < help->memory_binding_nodesets.size());
            char * sm;
            hwloc_bitmap_asprintf(&sm, help->memory_binding_nodesets[i/qm]);

            verbose_fmt_print(0, 2, "# {:<{}} {:<{}} {} [{} thread(s)] [ m:{} c:{} ]\n",
                    (i % qm) ? "" : tm[i / qm], m,
                    (i % qc) ? "" : tc[i / qc], c,
                    "job", // tj[i].c_str(),
                    nthreads_per_subjob,
                    sm, sc
                    );
            free(sm);
            free(sc);
        }
    }
#endif
    verbose_output_end_batch();
}/*}}}*/


#if 0
int las_parallel_desc::query_memory_binding(void * addr, size_t len);
{
    /* tentative specification. This is most probably a debug function.
     *
     * return an integer k for something that matches precisely the k-th
     * memory binding zone that we have defined.
     *
     * return -1 for something that isn't bound.
     *
     * return -2 if none of the above.
     */
    /* base this on one of the following hwloc api calls ?
    
    int hwloc_get_area_membind_nodeset (hwloc_topology_t topology,
            const void *addr, size_t len,
            hwloc_nodeset_t nodeset,
            hwloc_membind_policy_t *policy,
            int flags);
    int hwloc_get_area_membind (hwloc_topology_t topology,
            const void *addr, size_t len,
            hwloc_bitmap_t set,
            hwloc_membind_policy_t *policy,
            int flags);
    int hwloc_get_area_memlocation (hwloc_topology_t topology,
            const void *addr, size_t len,
            hwloc_bitmap_t set,
            int flags);

     */
    return -2;
}
#endif

int las_parallel_desc::set_loose_binding() const
{
#ifdef HAVE_HWLOC
    if (help->depth == 0)
        return 0;
    hwloc_obj_t root = hwloc_get_root_obj(help->topology);
    hwloc_nodeset_t n = root->nodeset;
    hwloc_cpuset_t c = root->cpuset;
    /* This achieves a **global** binding. It's not entirely clear we
     * ever want this in the normal course of an execution of las, in
     * fact.
     */
    int rc;
#if HWLOC_API_VERSION < 0x010b03
    /* this legacy called remained valid throughout hwloc 1.x */
    rc = hwloc_set_membind_nodeset(help->topology, n, HWLOC_MEMBIND_BIND,
            HWLOC_MEMBIND_THREAD |
            HWLOC_MEMBIND_STRICT);
#else
    /* newer call */
    rc = hwloc_set_membind(help->topology, n, HWLOC_MEMBIND_BIND,
            HWLOC_MEMBIND_THREAD |
            HWLOC_MEMBIND_STRICT |
            HWLOC_MEMBIND_BYNODESET);
#endif
    if (rc < 0) {
        char * s;
        hwloc_bitmap_asprintf(&s, n);
        if (errno == EXDEV) {
            fprintf(stderr, "Error, cannot enforce loose memory binding [ %s ]\n", s);
            /*
        } else {
            fprintf(stderr, "Error while attempting to set loose memory binding [ %s ]\n", s);
            */
        }
        free(s);
        return -1;
    }
    rc = hwloc_set_cpubind(help->topology, c,
            HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
    if (rc < 0) {
        char * s;
        hwloc_bitmap_asprintf(&s, c);
        if (errno == EXDEV) {
            fprintf(stderr, "Error, cannot enforce loose cpu binding [ %s ]\n", s);
        } else {
            fprintf(stderr, "Error while attempting to set loose cpu binding [ %s ]\n", s);
        }
        free(s);
        return -1;
    }
#endif
    return 0;
}

int las_parallel_desc::set_subjob_binding(int k) const
{
#ifdef HAVE_HWLOC
    if (help->depth == 0)
        return -1;
#endif
    set_subjob_cpu_binding(k);
    set_subjob_mem_binding(k);
    return 0;
}

int las_parallel_desc::set_subjob_mem_binding(int k MAYBE_UNUSED) const
{
#ifdef HAVE_HWLOC
    if (help->depth == 0)
        return -1;
    if (help->memory_binding_nodesets.empty())
        return -1;

    ASSERT_ALWAYS(0<= k && k < (int) help->subjob_binding_cpusets.size());
    int const m = k / number_of_subjobs_per_memory_binding_zone();
    ASSERT_ALWAYS(m < (int) help->memory_binding_nodesets.size());
#if HWLOC_API_VERSION < 0x010b03
    /* this legacy called remained valid throughout hwloc 1.x */
    int rc = hwloc_set_membind_nodeset(help->topology,
            help->memory_binding_nodesets[m],
            HWLOC_MEMBIND_BIND,
            HWLOC_MEMBIND_THREAD |
            HWLOC_MEMBIND_STRICT);
#else
    /* newer call */
    int const rc = hwloc_set_membind(help->topology,
            help->memory_binding_nodesets[m],
            HWLOC_MEMBIND_BIND,
            HWLOC_MEMBIND_THREAD |
            HWLOC_MEMBIND_STRICT |
            HWLOC_MEMBIND_BYNODESET);
#endif
    if (rc < 0) {
        char * s;
        hwloc_bitmap_asprintf(&s, help->memory_binding_nodesets[m]);
        if (errno == EXDEV) {
            fprintf(stderr, "Error, cannot enforce memory binding for job %d [ %s ]\n", k, s);
            /*
        } else {
            fprintf(stderr, "Error while attempting to set memory binding for job %d [ %s ]\n", k, s);
            */
        }
        free(s);
        return -1;
    }
#endif
    return 0;
}
int las_parallel_desc::set_subjob_cpu_binding(int k MAYBE_UNUSED) const
{
#ifdef HAVE_HWLOC
    if (help->depth == 0)
        return -1;
    ASSERT_ALWAYS(0<= k && k < (int) help->subjob_binding_cpusets.size());
    int const rc = hwloc_set_cpubind(help->topology, help->subjob_binding_cpusets[k], HWLOC_CPUBIND_THREAD |  HWLOC_CPUBIND_STRICT);
    if (rc < 0) {
        char * s;
        hwloc_bitmap_asprintf(&s, help->subjob_binding_cpusets[k]);
        if (errno == EXDEV) {
            fprintf(stderr, "Error, cannot enforce cpu binding for job %d [ %s ]\n", k, s);
        } else {
            fprintf(stderr, "Error while attempting to set cpu binding for job %d [ %s ]\n", k, s);
        }
        free(s);
        return -1;
    }
#endif
    return 0;
}

#if 0
int las_parallel_desc::number_of_threads_loose() const {
#ifdef HAVE_HWLOC
    if (help->depth)
        return help->replicate ? help->number_of(-1, 0) : help->memory_binding_size;
    else
#endif
        return nthreads_per_subjob;
}
#endif

#ifdef HAVE_HWLOC
cxx_hwloc_nodeset las_parallel_desc::current_memory_binding() const {
    return help->current_memory_binding();
}
#endif

