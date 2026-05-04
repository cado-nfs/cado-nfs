#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mutex>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>

#include "fmt/base.h"
#include "fmt/format.h"

#include "params.hpp"
#include "macros.h"
#include "version_info.h"
#include "verbose.hpp"
#include "utils_cxx.hpp"
#include "cxx_mpz.hpp"
#include "portability.h" // strdup // IWYU pragma: keep

using cado::params::cxx_param_list;

void cxx_param_list::print_usage(FILE *f) const
{
    const auto argv0 = cmdline_argv0[0];

    if (argv0 != nullptr)
        fmt::print(f, "Usage: {} <parameters>\n", argv0);

    if (!usage_header.empty())
        fmt::print(f, "{}\n", usage_header);

    fmt::print(f, "The available parameters are the following:\n");

    /* This is a copy, so that we can mangle the docs a little bit */
    auto full_doc = documentation;

    /* prepend "(switch) " to the documentation of all switch strings */
    for(auto const & [s, n] : switches) {
        std::string & v = full_doc[s];
        if (v.empty()) v = "UNDOCUMENTED";
        v = fmt::format("(switch) {}", v);
    }

    /* prepend "(alias -BLAH) " to the documentation of all alissed
     * options
     */
    for(auto const & a : aliases) {
        std::string & v = full_doc[a.second];
        if (v.empty()) v = "UNDOCUMENTED";
        v = fmt::format("(alias -{}) {}", a.first, v);
    }

    for(auto const & k : documentation_structure) {
        ASSERT_ALWAYS(!k.empty());
        if (k[0] == ' ') {
            fmt::print(f, "\n=== {} ===\n", k.substr(1));
            continue;
        }
        auto const & d = full_doc[k];
        std::string ks = k;
        for(auto const & dl : split(d, "\n")) {
            fmt::print(f, "{:4}{:<20} {}\n", "", ks, dl);
            ks = {};
        }
    }
}

template<typename... Args>
static void add_key(cxx_param_list & pl,
        std::string const & key,
        Args&& ...args)
{

    auto it = pl.p.find(key);
    cxx_param_list::parameter n { std::forward<Args>(args)... };

    // switches always count as parsed, of course.
    if (pl.switches.find(key) != pl.switches.end())
        n.parsed = true;

    if (it != pl.p.end()) {
        if (it->second.from > n.from) {
            /* ignore lower-priority parameers */
            return;
        } else if (it->second.from == n.from) {
            if (pl.p[key].value == n.value) {
                pl.p[key].seen += n.seen;
                return;
            }
        }
    }
    pl.p[key] = n;
}

void cxx_param_list::add_key(
        std::string const & key, std::string const & value, enum origin o)
{
    ::add_key(*this, key, value, o);
}

/* If step_on_empty_line is non-zero, then this function reads the file until a
 * line containing only space caracters (check with isspace) is found. It allows
 * to read a file contaning more than one polynomial.
 * Otherwise the function reads the whole file.
 */
int cxx_param_list::read(std::istream & is, bool stop_on_empty_line)
{
    int all_ok=1;
    for(std::string line ; std::getline(is, line, '\n') ; ) {
        if (line[0] == '#')
            continue;

        // remove possible comment at end of line.
        std::string::size_type p = line.find('#');
        if (p != std::string::npos)
            line.erase(p);

        for( ; !line.empty() && isspace(line[line.size()-1]) ; )
            line.erase(line.size()-1);

        for(p = 0 ; !line.empty() && isspace(line[p]) ; p++);
        if (p) line.erase(0, p);

        // empty ps are ignored (unless stop_on_empty_line is non-zero).
        if (line.empty()) {
            if (stop_on_empty_line)
                break;
            else
                continue;
        }

#if 0
        // look for a left-hand-side. We grok anything that *BEGINS WITH
        // A DIGIT* as something that goes with the "NULL" token in the
        // pl dictionary. That looks like a pretty obscure hack, in fact.
        // Do we ever use it ?
        if (!isalpha((int)(unsigned char)line[0]) && line[0] != '_' && line[0] != '-') {
            param_list_add_key(pl, nullptr, line.c_str(), PARAMETER_FROM_FILE);
            continue;
        }
#endif
        std::string::const_iterator q = line.begin();
        for( ; q != line.end() && (isalnum((int)(unsigned char)*q) || *q == '_' || *q == '-') ; ++q);
        if (q == line.begin()) {
            fprintf(stderr, "Parse error, no usable key for config line:\n%s\n",
                    line.c_str());
            all_ok=0;
            continue;
        }

        const std::string key(line.cbegin(), q);

        /* Now we can match (whitespace+ | whitespace* separator whitespace*) data
         */
        for( ; q != line.end() && isspace((int)(unsigned char)*q) ; ++q);

        /* should we actually allow it, after all ? */
        if (q == line.end()) {
            fprintf(stderr, "Parse error, key with no value in config:\n%s\n",
                    line.c_str());
            continue;
        }

        /* match separator, which is one of : = := */
        if (*q == '=') {
            q++;
        } else if (*q == ':') {
            q++;
            if (*q == '=')
                q++;
        } else if (q == line.begin() + (ptrdiff_t) key.size()) {
            fprintf(stderr, "Parse error, no separator for config line:\n%s\n",
                    line.c_str());
            all_ok=0;
            continue;
        }
        for( ; *q && isspace((int)(unsigned char)*q) ; q++);

        const std::string value(q, line.cend());

        add_key(key, value, origin::FROM_FILE);
    }

    return all_ok;
}



void cxx_param_list::configure_alias(std::string const & key, std::string const & alias)
{
    auto ckey = drop_one_or_two_leading_dashes(key);
    auto calias = drop_one_or_two_leading_dashes(alias);

    ASSERT_ALWAYS(!calias.empty());
    ASSERT_ALWAYS(!ckey.empty());

    if (use_doc && !is_documented(ckey))
        fmt::print(stderr, "# Warning: an alias {} is declared to the key {}, which is undocumented\n", calias, ckey);


    aliases[calias] = ckey;
}

void cxx_param_list::configure_switch(std::string const & switchname)
{
    auto cswitchname = drop_one_or_two_leading_dashes(switchname);
    
    if (use_doc && !is_documented(cswitchname))
        fmt::print(stderr, "# Warning: a switch {} is declared but is undocumented\n", switchname);

    switches.emplace(cswitchname, nullptr);
}

void cxx_param_list::configure_switch_old(std::string const & switchname, int * p)
{
    auto cswitchname = drop_one_or_two_leading_dashes(switchname);
    
    if (use_doc && !is_documented(cswitchname))
        fmt::print(stderr, "# Warning: a switch {} is declared but is undocumented\n", switchname);

    switches.emplace(cswitchname, p);
}

int cxx_param_list::update_cmdline(
        int & argc, char const ** & argv)
{
    ASSERT_ALWAYS(argv);
    if (!cmdline_argv0) {
        cmdline_argv0 = argv-1;
        cmdline_argc0 = argc+1;
    }
    if (argc == 0)
        return 0;

    const char * arg = argv[0];
    if (strcmp(arg, "--") == 0) {
        /* return, but let the caller recognize "--" as the reason for
         * breaking out early */
        return 0;
    }

    bool has_dashes = false;
    for(int i = 0 ; i < 2 && *arg == '-' ; arg++, i++, has_dashes = true) ;
    std::string arg2 = arg;
    std::string rhs;
    for(auto it = arg2.begin() ; it != arg2.end() ; ++it) {
        if (*it == '=') {
            rhs = std::string(it + 1, arg2.end());
            arg2.erase(it, arg2.end());
            break;
        }
    }

    /* Note that we now support --option=X as well */

    bool negate = false;

    if (has_dashes) {
        if (arg2.starts_with("no-") || arg2.starts_with("no_")) {
            arg2 = arg2.substr(3);
            negate = true;
        }

        /* negate + having a rhs does not make sense */
        if (negate && !rhs.empty())
            return 0;
    }

    {
        auto const it = aliases.find(arg2);
        if (it != aliases.end())
            arg2 = it->second;
    }

    auto const it = switches.find(arg2);
    if (it != switches.end()) {
        /* A switch should really follow most of the standard logic for
         * the rest of parsing, so it makes sense to store "true" or "1"
         * in the dictionary, so that the boolean parser is happy.
         * Actually, "1" is more flexible because it still makes it
         * possible to have cumulative switches.
         *
         * And things like --v=false are taken to mean what we imagine.
         */

        int v = 0;
        if (!rhs.empty()) {
            if (!cado::params::parse(rhs, v)) { // E: No matching member function …
                bool w;
                if (!cado::params::parse(rhs, w)) { // E: No matching member funct…
                    fail("Unexpected value {} stored for switch {}", // E: …
                            rhs, arg2);
                }
                v = w;
            }
            add_key(arg2, rhs, origin::FROM_CMDLINE);
        } else if (!negate) {
            /* potentially take the already store value, and increase it
             */
            if (auto const jt = p.find(arg2); jt != p.end()) {
                if (!cado::params::parse(jt->second.value, v)) {
                    bool w;
                    if (!cado::params::parse(jt->second.value, w)) {
                        fail("Unexpected value {} stored for switch {}",
                                jt->second.value, arg2);
                    }
                    v = w;
                }
            }
            v++;
            add_key(arg2, fmt::format("{}", v), origin::FROM_CMDLINE);
        } else {
            add_key(arg2, "0", origin::FROM_CMDLINE);
        }

        /* BEGIN LEGACY BLOCK -- will be removed */
        if (it->second)
            *it->second = v;
        /* END LEGACY BLOCK */
        /* add it to the dictionary, so that
         * param_list_parse_switch can find it later on.
         */

        argv += 1;
        argc -= 1;
        return 1;
    } else {
        /* negate only makes sense for a switch! */
        if (negate)
            return 0;
    }

    if (has_dashes && rhs.empty()) {
        /* If we have no RHS yet, take it from the rest of the command
         * line IF AND ONLY IF our key starts with dashes. No-dash format
         * only groks a rhs introduced with an = sign. */

        if (argc < 2)
            return 0;

        rhs = argv[1];
        argv += 1;
        argc -= 1;
    }

    if (!rhs.empty()) {
        add_key(arg2, rhs, origin::FROM_CMDLINE);
        argv += 1;
        argc -= 1;
        return 1;
    }
    return 0;
}

std::string cado::params::collect_command_line(int argc, char const *argv[])
{
    return join(argv, argv + argc, " ");
}

void cxx_param_list::process_command_line(
        int & argc, char const ** & argv, bool accept_trailing_args)
{
    cmdline_argv0 = argv;
    cmdline_argc0 = argc;

    argc--, argv++;

    for (; argc;) {
        if (update_cmdline(argc, argv)) {
            continue;
        }
        if (strcmp(*argv, "--help") == 0) {
            print_usage(stderr);
            exit(EXIT_SUCCESS);
        } else if (!accept_trailing_args) {
            fail("unexpected argument: {}\n", argv[0]);
        } else {
            /* extra arguments are accepted, we can break here. */
            break;
        }
    }
}

void cxx_param_list::process_command_line_and_extra_parameter_files(
        int & argc, char const ** & argv)
{
    cmdline_argv0 = argv;
    cmdline_argc0 = argc;

    argc--, argv++;

    for (; argc;) {
        if (update_cmdline(argc, argv)) {
            continue;
        }
        if (strcmp(*argv, "--help") == 0) {
            print_usage(stderr);
            exit(EXIT_SUCCESS);
        } else {
            /* Could also be a file */
            std::ifstream f(argv[0]);
            if (!f)
                fail("Cannot read {}", argv[0]);
            read(f);
            argv++,argc--;
            continue;
        }
        fail("unexpected argument: {}\n", argv[0]);
    }
}



/* Look up an entry in a param_list, and update the parsed flag. It does
   mutex locking to make look-ups thread safe; the caller must not access
   any param_list entries by itself. */
std::string const *
cxx_param_list::get_assoc_ptr(std::string const & key0, bool stealth, bool * seen)
{
    auto key = drop_one_or_two_leading_dashes(key0);
    const std::scoped_lock dummy(mutex);
    if (use_doc && !stealth && !is_documented(key))
        fmt::print(stderr, "# Warning: parameter {} is checked by this program but is undocumented.\n", key);
    auto it = p.find(key);
    if (it == p.end())
        return nullptr;

    it->second.parsed = true;
    if (seen) *seen = it->second.seen;
    return &it->second.value;
}

bool cxx_param_list::get_assoc(std::string const & key, std::string & value, bool stealth, bool * seen)
{
    const std::string * t = get_assoc_ptr(key, stealth, seen);
    if (t)
        value = *t;
    return t != nullptr;
}

int cxx_param_list::warn_unused() const
{
    int u = 0;
    for(auto const & [ key, par ] : p) {
        if (!par.parsed && par.from != origin::FROM_FILE) {
            fmt::print(stderr, "Warning: unused command-line parameter {}\n",
                    key);
            u++;
        }
    }
    return u;
}

void cxx_param_list::print_command_line(FILE * stream) const
{
    if (!cmdline_argv0)
        return;

    if (verbose_enabled(CADO_VERBOSE_PRINT_CMDLINE)) {
        /* print command line */
        fmt::print(stream, "# ({}) {}\n", cado_revision_string,
                collect_command_line(cmdline_argc0, cmdline_argv0));
    }
    if (verbose_enabled(CADO_VERBOSE_PRINT_MODIFIED_FILES)) {
        if (strlen(cado_modified_files) > 1)
          fmt::print(stream, "# List of modified files in working directory and "
                   "their SHA1 sum:\n{}", cado_modified_files);
    }
    if (verbose_enabled(CADO_VERBOSE_PRINT_COMPILATION_INFO)) {
#ifdef  __GNUC__
#ifndef __ICC
        fmt::print(stream, "# Compiled with gcc " __VERSION__ "\n");
#else
        /* icc defines __GNUC__ too */
        fmt::print(stream, "# Compiled with icc {}.{}.{} (gcc version {}.{}.{} compatibility)\n",
                 __ICC / 100, __INTEL_COMPILER % 100, __INTEL_COMPILER_UPDATE,
                 __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

        /* Apple Clang (based on llvm) identifies itself as a flavour of GNUC 4.2.1
           but seems to compile CADO-NFS properly 
           $ echo | clang -dD -E -pipe -
# 1 "<stdin>"
# 1 "<stdin>" 1
# 1 "<built-in>" 1
# 1 "<built-in>" 3
#define __llvm__ 1
#define __clang__ 1
#define __clang_major__ 4
#define __clang_minor__ 1
#define __clang_patchlevel__ 0
#define __clang_version__ "4.1 ((tags/Apple/clang-421.11.66))"
#define __GNUC_MINOR__ 2
#define __GNUC_PATCHLEVEL__ 1
#define __GNUC__ 4
...
*/

#if GNUC_VERSION_ATLEAST(4,1,2) && GNUC_VERSION_ATMOST(4,2,2) \
        && ! (__llvm__ || __clang__)
        fmt::print(stream, "# WARNING: this version of GCC is known to miscompile CADO-NFS. See https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/14490\n");
#endif
#endif
        fmt::print(stream, "# Compilation flags (C) " CFLAGS "\n");
        fmt::print(stream, "# Compilation flags (C++) " CXXFLAGS "\n");
    }
}






size_t cxx_param_list::get_list_count(std::string const & key, std::string const & sep)
{
    if (auto const * t = has(drop_one_or_two_leading_dashes(key)))
        return split(*t, sep).size();
    else
        return 0;
}

/***********************************************************************/
/* compatibility calls */

int param_list_parse_mpz(cxx_param_list & pl, const char * key, mpz_ptr x)
{
    cxx_mpz u;
    int const r = pl.parse(key, u);
    if (r)
        mpz_set(x, u);
    return r;
}
int cxx_param_list::read(FILE *f, bool stop_on_empty_line)
{
    int all_ok=1;
    const int linelen = 2048;
    char line[linelen];
    while (!feof(f)) {
        if (fgets(line, linelen, f) == nullptr)
            break;
        if (line[0] == '#')
            continue;
        // remove possible comment at end of line.
        char * hash;
        if ((hash = strchr(line, '#')) != nullptr) {
            *hash = '\0';
        }

        char * p = line;

        // trailing space
        auto l = strlen(p);
        for( ; l && isspace((int)(unsigned char)p[l-1]) ; l--);
        p[l] = '\0';

        // leading space.
        for( ; *p && isspace((int)(unsigned char)*p) ; p++, l--);

        // empty ps are ignored (unless stop_on_empty_line is non-zero).
        if (l == 0)
        {
          if (stop_on_empty_line)
            break;
          else
            continue;
        }

        // look for a left-hand-side. We grok anything that *BEGINS WITH
        // A DIGIT* as something that goes with the "NULL" token in the
        // pl dictionary. That looks like a pretty obscure hack, in fact.
        // Do we ever use it ?
        l = 0;
        if (!isalpha((int)(unsigned char)p[l]) && p[l] != '_' && p[l] != '-') {
            ASSERT_ALWAYS(0);
            add_key("", line, origin::FROM_FILE);
            continue;
        }
        for( ; p[l] && (isalnum((int)(unsigned char)p[l]) || p[l] == '_' || p[l] == '-') ; l++);

        auto const lhs_length = l;

        if (lhs_length == 0) {
            fprintf(stderr, "Parse error, no usable key for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }

        /* Now we can match (whitespace*)(separator)(whitespace*)(data)
         */
        char * q = p + lhs_length;
        for( ; *q && isspace((int)(unsigned char)*q) ; q++);

        /* match separator, which is one of : = := */
        if (*q == '=') {
            q++;
        } else if (*q == ':') {
            q++;
            if (*q == '=')
                q++;
        } else if (q == p + lhs_length) {
            fprintf(stderr, "Parse error, no separator for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }
        for( ; *q && isspace((int)(unsigned char)*q) ; q++);

        add_key(std::string(p, lhs_length), q, origin::FROM_FILE);
    }

    return all_ok;
}

