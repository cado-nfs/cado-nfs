#include "cado.h" // IWYU pragma: keep
#include <ios>
#include <cctype>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>              /* isdigit isspace */
#include <limits.h>             /* INT_MIN INT_MAX */
#include <errno.h>
#include <stdarg.h>
#include <inttypes.h>
#include <mutex>
#include <stdio.h>
#include <gmp.h>
#include <type_traits>
#include <ios>
#include "fmt/core.h"
#include "fmt/format.h"
#include "cxx_mpz.hpp"

#include "params.h"
#include "macros.h"
#include "version_info.h"
#include "verbose.h"
#include "portability.h" // strdup // IWYU pragma: keep

typedef int (*sortfunc_t) (const void *, const void *);

static std::mutex mutex;

void param_list_init(param_list_ptr pl)
{
    pl->pimpl = new param_list_impl();
}

void param_list_set(param_list_ptr pl, param_list_srcptr pl0)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    param_list_impl & pli0 = *static_cast<param_list_impl *>(pl0->pimpl);
    pli = pli0;
}

void param_list_swap(param_list_ptr pl, param_list_ptr pl0)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    param_list_impl & pli0 = *static_cast<param_list_impl *>(pl0->pimpl);
    std::swap(pli, pli0);
}

void param_list_clear(param_list_ptr pl)
{
    delete static_cast<param_list_impl *>(pl->pimpl);
}

int param_list_empty(param_list_srcptr pl)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    return pli.p.empty();
}

void param_list_usage_header(param_list_ptr pl, const char * hdr, ...)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    va_list ap;
    va_start(ap, hdr);
    char * tmp;
    int rc = vasprintf(&tmp, hdr, ap);
    ASSERT_ALWAYS(rc >= 0);
    pli.usage_header = std::string(tmp);
    free(tmp);
    va_end(ap);
}

void param_list_decl_usage(param_list_ptr pl, const char * key, const char * doc)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    pli.documentation[key] = doc;
    /* Note that duplicate calls to param_list_decl_usage for the same
     * key will not trigger two distinct prints of the same
     * documentation string. See the collapsing logic in
     * param_list_print_usage.
     */
    pli.use_doc = true;
}


static int is_documented_key(param_list_srcptr pl, const char *key) {
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    return pli.documentation.find(key) != pli.documentation.end();
}

void param_list_print_usage(param_list_srcptr pl, const char * argv0, FILE *f)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);

    if (argv0 != NULL)
        fprintf(f, "Usage: %s <parameters>\n", argv0);

    if (!pli.usage_header.empty())
        fputs(pli.usage_header.c_str(), f);

    fprintf(f, "The available parameters are the following:\n");

    /* This is a copy, so that we can mangle the docs a little bit */
    auto full_doc = pli.documentation;

    /* prepend "(switch) " to the documentation of all switch strings */
    for(auto const & s : pli.switches) {
        std::string & v = full_doc[s.first];
        if (v.empty()) v = "UNDOCUMENTED";
        v = std::string("(switch) ") + v;
    }

    /* prepend "(alias -BLAH) " to the documentation of all alissed
     * options
     */
    for(auto const & a : pli.aliases) {
        std::string & v = full_doc[a.second];
        if (v.empty()) v = "UNDOCUMENTED";
        v = std::string("(alias -") + a.first + ") " + v;
    }

    for(auto const & d : full_doc) {
        fprintf(f, "    -%-*s %s\n", 20, d.first.c_str(), d.second.c_str());
    }
}

template<typename... Args>
void param_list_add_key(param_list_impl & pli,
        std::string const & key,
        Args&& ...args)
{

    auto it = pli.p.find(key);
    param_list_impl::parameter p { std::forward<Args>(args)... };

    // switches always count as parsed, of course.
    if (pli.switches.find(key) != pli.switches.end())
        p.parsed = true;

    if (it != pli.p.end()) {
        if (it->second.origin > p.origin) {
            /* ignore lower-priority parameers */
            return;
        } else if (it->second.origin == p.origin) {
            if (pli.p[key].value == p.value) {
                pli.p[key].seen += p.seen;
                return;
            }
        }
    }
    pli.p[key] = p;
}

void param_list_add_key(param_list_ptr pl,
        const char * key, const char * value, enum parameter_origin o)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    param_list_add_key(pli, key, value, o);
}


void param_list_remove_key(param_list_ptr pl, const char * key)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    auto it = pli.p.find(key);
    if (it != pli.p.end())
        pli.p.erase(it);
}

/* If step_on_empty_line is non-zero, then this function reads the file until a
 * line containing only space caracters (check with isspace) is found. It allows
 * to read a file contaning more than one polynomial.
 * Otherwise the function reads the whole file.
 */
int param_list_read_stream(param_list_ptr pl, FILE *f, int stop_on_empty_line)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    int all_ok=1;
    const int linelen = 2048;
    char line[linelen];
    while (!feof(f)) {
        if (fgets(line, linelen, f) == NULL)
            break;
        if (line[0] == '#')
            continue;
        // remove possible comment at end of line.
        char * hash;
        if ((hash = strchr(line, '#')) != NULL) {
            *hash = '\0';
        }

        char * p = line;

        // trailing space
        int l = strlen(p);
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
        if (!(isalpha((int)(unsigned char)p[l]) || p[l] == '_' || p[l] == '-')) {
            param_list_add_key(pl, NULL, line, PARAMETER_FROM_FILE);
            continue;
        }
        for( ; p[l] && (isalnum((int)(unsigned char)p[l]) || p[l] == '_' || p[l] == '-') ; l++);

        int lhs_length = l;

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

        param_list_add_key(pli, std::string(p, lhs_length), q, PARAMETER_FROM_FILE);
    }

    return all_ok;
}

int param_list_read_file(param_list_ptr pl, const char * name)
{
    FILE * f;
    f = fopen(name, "r");
    if (f == NULL) {
        fprintf(stderr, "Cannot read %s\n", name);
        exit(1);
    }
    int r = param_list_read_stream(pl, f, 0);
    fclose(f);
    return r;
}

int param_list_configure_alias(param_list_ptr pl, const char * key, const char * alias)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    for(int i = 0 ; i < 2 && *key == '-' ; key++, i++) ;
    for(int i = 0 ; i < 2 && *alias == '-' ; alias++, i++) ;

    if (pli.use_doc) {
        if (!is_documented_key(pl, key))
            fprintf(stderr, "# Warning: an alias %s is declared to the key %s, which is undocumented\n", alias, key);
    }

    ASSERT_ALWAYS(alias != NULL);
    ASSERT_ALWAYS(key != NULL);

    pli.aliases[alias] = key;
    return 0;
}

int param_list_configure_switch(param_list_ptr pl, const char * switchname, int * ptr)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    for(int i = 0 ; i < 2 && *switchname == '-' ; switchname++, i++) ;
    
    if (pli.use_doc)
        if (!is_documented_key(pl, switchname))
            fprintf(stderr, "# Warning: a switch %s is declared but is undocumented\n", switchname);

    pli.switches[switchname] = ptr;
    return 0;
}

int param_list_update_cmdline(param_list_ptr pl,
        int * p_argc, char *** p_argv)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    ASSERT_ALWAYS(*p_argv != NULL);
    if (!pli.cmdline_argv0) {
        pli.cmdline_argv0 = (*p_argv)-1;
        pli.cmdline_argc0 = (*p_argc)+1;
    }
    if (*p_argc == 0)
        return 0;

    const char * arg = (*p_argv)[0];
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

    if (has_dashes && !rhs.empty()) {
        /* Maybe we should support this, after all? Could -v=5 be
         * something useful?
         */
        return 0;
    }

    bool negate = false;

    if (has_dashes) {
        if (arg2.substr(0, 3) == "no-" || arg2.substr(0, 3) == "no_") {
            arg2 = arg2.substr(3);
            negate = true;
        }
    }

    for(auto const & a : pli.aliases) {
        if (a.first == arg2) {
            arg2 = a.second;
        }
    }
    
    if (has_dashes) {
        ASSERT_ALWAYS(rhs.empty());

        for(auto & s : pli.switches) {
            if (s.first == arg2) {
                if (s.second) {
                    /* Note that a switch with a NULL pointer is allowed
                     */
                    if (negate)
                        *s.second = 0;
                    else
                        ++*s.second;
                }
                /* add it to the dictionary, so that
                 * param_list_parse_switch can find it later on.
                 */
                param_list_add_key(pli, arg2, std::string(), PARAMETER_FROM_CMDLINE);
                (*p_argv)+=1;
                (*p_argc)-=1;
                return 1;
            }
        }

        /* We have dashes, we're not a switch. We have no RHS yet, take
         * it from the rest of the command line */

        if (*p_argc < 2)
            return 0;

        rhs = (*p_argv)[1];
        (*p_argv)+=1;
        (*p_argc)-=1;
    }
    if (!rhs.empty()) {
        param_list_add_key(pli, arg2, rhs, PARAMETER_FROM_CMDLINE);
        (*p_argv)+=1;
        (*p_argc)-=1;
        return 1;
    }
    return 0;
}

/* Look up an entry in a param_list, and update the parsed flag. It does
   mutex locking to make look-ups thread safe; the caller must not access
   any param_list entries by itself. */
static const char *
get_assoc_ptr(param_list_ptr pl, const char * key, bool stealth = false, int * const seen = NULL)
{
    param_list_impl & pli = *static_cast<param_list_impl *>(pl->pimpl);
    std::lock_guard<std::mutex> dummy(mutex);
    if (pli.use_doc) {
        for(int i = 0 ; i < 2 && *key == '-' ; key++, i++) ;
        if (!stealth && !is_documented_key(pl, key)) 
            fprintf(stderr, "# Warning: parameter %s is checked by this program but is undocumented.\n", key);
    }
    auto it = pli.p.find(key);
    if (it == pli.p.end())
        return NULL;

    it->second.parsed = 1;
    if (seen) *seen = it->second.seen;
    return it->second.value.data();
}

static int
get_assoc(param_list_ptr pl, const char * key, std::string & value, bool stealth = false, int * const seen = NULL)
{
    const char * t = get_assoc_ptr(pl, key, stealth, seen);
    if (t)
        value = t;
    return t != NULL;
}


struct parameter_exception : public std::runtime_error {
    parameter_exception(std::string const & what) : std::runtime_error(what) {}
};

/* recognize any integer in one of the following forms:
 *      <some integer>
 *      <somd double that casts properly to T **AND BACK**>
 *
 * note that we want the full string to match!
 *
 */

template<typename T> struct parse {
    bool operator()(std::string const & s, T & value) const
    {
        std::istringstream ss(s);
        if (ss >> value && ss.eof()) {
            return true;
        }
        double d;
        ss.str(s);
        if (ss >> d && ss.eof()) {
            value = T(d);
            if (double(value) == d)
                return true;
            else
                throw parameter_exception("incorrect double-to-integer conversion in parameter");
        }
        return false;
    }
};

template<> struct parse<double> { 
    bool operator()(std::string const & s, double & value) const
    {
        std::istringstream ss(s);
        return ss >> value && ss.eof();
    }
};

template<> struct parse<std::string> {
    bool operator()(std::string const & s, std::string & value) const
    {
        value = s;
        return true;
    }
};

template<typename T>
bool
parse_list(std::string const & s, std::vector<T> & value, std::string const & sep)
{
    std::vector<std::string> tokens;
    value.clear();

    for(size_t pos0 = 0, pos ; pos0 < s.size() ; pos0 = pos) {
        pos = s.find(sep, pos0);
        if (pos == std::string::npos) {
            tokens.push_back(s.substr(pos0));
            break;
        }
        tokens.push_back(s.substr(pos0, pos - pos0));
        pos += sep.size();
        if (sep == " ")
            for( ; pos < s.size() && std::isspace(s[pos]) ; pos++) ;
    }
    for(auto const & t : tokens) {
        T v;
        if (!parse<T>()(t, v))
            return false;
        value.push_back(v);
    }
    return true;
}

template<typename T>
struct parse<std::vector<T>> {
    bool
        operator()(std::string const & s, std::vector<T> & value) const
        {
            return parse_list<T>(s, value, ",");
        }
};


template<> struct parse<cxx_mpz_poly> {
    bool operator()(std::string const & s, cxx_mpz_poly & value) const
    {
        if (s.find(',') == std::string::npos) {
            return mpz_poly_set_from_expression(value, s.c_str());
        }
        std::vector<cxx_mpz> blah;
        if (!parse<std::vector<cxx_mpz>>()(s, blah))
            return false;
        mpz_poly_set_zero(value);
        for(int i = 0 ; i < (int) blah.size() ; i++) {
            mpz_poly_setcoeff(value, i, blah[i]);
        }
        return true;
    }
};

template<> struct parse<cxx_mpz> {
    bool operator()(std::string const & s, cxx_mpz & value) const
    {
        unsigned int nread;
        int rc = gmp_sscanf(s.c_str(), "%Zi%n", (mpz_ptr) value, &nread);
        if (rc != 1 || nread != s.size()) {
            /* also recognize integers written as floating-point, like 6.3e8
            */
            if (s[nread] == '.' || s[nread] == 'e') {
                mpf_t zf;
                mpf_init(zf);
                rc = gmp_sscanf(s.c_str(), "%Ff%n", zf, &nread);
                rc = (rc == 1 && nread == s.size() && mpf_integer_p(zf));
                if (rc) mpz_set_f(value, zf);
                mpf_clear(zf);
                if (rc) return true;
            }

            if (mpz_set_from_expression(value, s.c_str()))
                return true;

            throw parameter_exception("incorrect conversion to mpz");
        }
        return true;
    }
};

template<typename T>
int
param_list_parse_inner(param_list_ptr pl, const char * key, T & r, bool stealth = false)
{
    std::string value;
    std::string diagnostic;
    if (!get_assoc(pl, key, value, stealth))
        return 0;
    try {
        if (parse<T>()(value, r))
            return true;
    } catch (parameter_exception const & e) {
        diagnostic = e.what();
    }

    std::ostringstream ss;
    ss << "cannot cast parameter " << key << " to type " << typeid(T).name();
    if (!diagnostic.empty())
        ss << ": " << diagnostic;
    throw parameter_exception(ss.str());
}

template<typename T>
int
param_list_parse(param_list_ptr pl, const char * key, T & r)
{
    return param_list_parse_inner<T>(pl, key, r);
}


template int param_list_parse<int>(param_list_ptr pl, const char * key, int & r);
template int param_list_parse<unsigned int>(param_list_ptr pl, const char * key, unsigned int & r);
#ifndef LONG_IS_EXACTLY_INT
template int param_list_parse<long>(param_list_ptr pl, const char * key, long & r);
template int param_list_parse<unsigned long>(param_list_ptr pl, const char * key, unsigned long & r);
#endif

#ifndef INT64_T_IS_EXACTLY_LONG
template int param_list_parse<int64_t>(param_list_ptr pl, const char * key, int64_t & r);
template int param_list_parse<uint64_t>(param_list_ptr pl, const char * key, uint64_t & r);
#endif

template int param_list_parse<double>(param_list_ptr pl, const char * key, double & r);
template int param_list_parse<std::vector<int>>(param_list_ptr pl, const char * key, std::vector<int> & r);
template int param_list_parse<std::vector<std::string>>(param_list_ptr pl, const char * key, std::vector<std::string> & r);

template int param_list_parse<std::string>(param_list_ptr pl, const char * key, std::string & r);

int param_list_parse_long(param_list_ptr pl, const char * key, long * r)
{
    return param_list_parse<long>(pl, key, *r);
}

int param_list_parse_ulong(param_list_ptr pl, const char * key, unsigned long * r)
{
    return param_list_parse<unsigned long>(pl, key, *r);
}

int param_list_parse_size_t(param_list_ptr pl, const char * key, size_t * r)
{
    return param_list_parse<size_t>(pl, key, *r);
}

int param_list_parse_int64(param_list_ptr pl, const char * key, int64_t * r)
{
    return param_list_parse<int64_t>(pl, key, *r);
}

int param_list_parse_uint64(param_list_ptr pl, const char * key, uint64_t * r)
{
    return param_list_parse<uint64_t>(pl, key, *r);
}

int param_list_parse_int(param_list_ptr pl, const char * key, int * r)
{
    return param_list_parse<int>(pl, key, *r);
}

int param_list_parse_uint(param_list_ptr pl, const char * key, unsigned int * r)
{
    return param_list_parse<unsigned int>(pl, key, *r);
}

int param_list_parse_double(param_list_ptr pl, const char * key, double * r)
{
    return param_list_parse<double>(pl, key, *r);
}

int param_list_parse_string(param_list_ptr pl, const char * key, char * r, size_t n)
{
    std::string value;
    if (!get_assoc(pl, key, value))
        return 0;
    if (r && strlcpy(r, value.c_str(), n) > n) {
        throw parameter_exception { fmt::format(FMT_STRING(
                " parameter for key {} does not fit within string buffer"
                " of length {}\n"), key, n) };
    }
    return 1;
}

int param_list_parse_mpz(param_list_ptr pl, const char * key,
    mpz_ptr f)
{
    cxx_mpz P;
    int b = param_list_parse(pl, key, P);
    if (b)
        mpz_set(f, P);
    return b;
}

int param_list_parse_mpz_poly(param_list_ptr pl, const char * key,
    mpz_poly_ptr f)
{
    cxx_mpz_poly P;
    int b = param_list_parse(pl, key, P);
    if (b)
        mpz_poly_set(f, P);
    return b;
}


template<typename T>
int param_list_parse_pair_with_and(param_list_ptr pl, const char * key, T * res, const char * sep)
{
    std::vector<T> r;
    std::string value;
    std::string diagnostic;
    if (!get_assoc(pl, key, value))
        return 0;
    try {
        if (parse_list<T>(value, r, sep) && r.size() == 2) {
            if (res) res[0] = r[0];
            if (res) res[1] = r[1];
            return true;
        }
    } catch (parameter_exception const & e) {
        diagnostic = e.what();
    }
    std::ostringstream ss;
    ss << "cannot cast parameter " << key << " to pairs of elements of type " << typeid(T).name();
    if (!diagnostic.empty())
        ss << ": " << diagnostic;
    throw parameter_exception(ss.str());
}

int param_list_parse_uint_and_uint(param_list_ptr pl, const char * key, unsigned int * r, const char * sep)
{
    return param_list_parse_pair_with_and<unsigned int>(pl, key, r, sep);
}

int param_list_parse_int_and_int(param_list_ptr pl, const char * key, int * r, const char * sep)
{
    return param_list_parse_pair_with_and<int>(pl, key, r, sep);
}

int param_list_parse_ulong_and_ulong(param_list_ptr pl, const char * key, unsigned long * r, const char * sep)
{
    return param_list_parse_pair_with_and<unsigned long>(pl, key, r, sep);
}

int param_list_parse_uint64_and_uint64(param_list_ptr pl, const char * key, uint64_t * r, const char * sep)
{
    return param_list_parse_pair_with_and<uint64_t>(pl, key, r, sep);
}

int param_list_parse_intxint(param_list_ptr pl, const char * key, int * r)
{
    return param_list_parse_int_and_int(pl, key, r, "x");
}

int param_list_parse_double_and_double(param_list_ptr pl, const char * key, double * r, const char * sep)
{
    return param_list_parse_pair_with_and<double>(pl, key, r, sep);
}


template<typename T>
int param_list_parse_list(param_list_ptr pl, const char * key, std::vector<T> & res, const char * sep, bool stealth = false)
{
    std::string value;
    std::string diagnostic;
    if (!get_assoc(pl, key, value, stealth))
        return 0;
    try {
        if (parse_list<T>(value, res, sep))
            return true;
    } catch (parameter_exception const & e) {
        diagnostic = e.what();
    }
    std::ostringstream ss;
    ss << "cannot cast parameter " << key << " to list of elements of type " << typeid(T).name();
    if (!diagnostic.empty())
        ss << ": " << diagnostic;
    throw parameter_exception(ss.str());
}

#if 0
/* shouldn't exist, really */
int param_list_parse_string_list_alloc(param_list_ptr pl, const char * key, char *** r, int * n, const char * sep)
{
    char * value;
    *r = NULL;
    *n = 0;
    int parsed = 0;
    if (!get_assoc(pl, key, &value, NULL))
        return 0;
    for( ; ; ) {
        char * v = strstr(value, sep);
        int itemsize;
        if (v == NULL) {
            itemsize = strlen(value);
        } else {
            itemsize = v - value;
        }
        *r = realloc(*r, (parsed + 1) * sizeof(char *));
        (*r)[parsed++] = strndup(value, itemsize);
        if (!v)
            break;
        value = v + strlen(sep);
    }
    *n = parsed;
    return parsed;
}
#endif

size_t param_list_get_list_count(param_list_ptr pl, const char * key)
{
    std::vector<std::string> v;
    if (!param_list_parse(pl, key, v))
        return 0;
    return v.size();
}

/* this returns the size of the list, it's a bit awkward */
template<typename T>
int param_list_parse_raw_fixed_list(param_list_ptr pl, const char * key, T * r, size_t n, const char * sep)
{
    std::vector<T> v;
    if (!param_list_parse_list(pl, key, v, sep))
        return 0;

    if (v.size() > n) {
        throw parameter_exception { fmt::format(FMT_STRING(
                    " parameter for key {} does not fit in fixed list"
                    " of length {}\n"), key, n) };
    }
    std::copy(v.begin(), v.end(), r);
    return v.size();
}

template<typename T>
int param_list_parse_raw_variable_list(param_list_ptr pl, const char * key, T ** r, size_t * n, const char * sep)
{
    std::vector<T> v;
    if (!param_list_parse_list(pl, key, v, sep))
        return 0;

    *r = new T[v.size()];
    *n = v.size();
    std::copy(v.begin(), v.end(), *r);
    return v.size();
}

int param_list_parse_int_list(param_list_ptr pl, const char * key,
    int * r, size_t n, const char * sep)
{
    return param_list_parse_raw_fixed_list(pl, key, r, n, sep);
}

int param_list_parse_uint_list(param_list_ptr pl, const char * key,
    unsigned int * r, size_t n, const char * sep)
{
    return param_list_parse_raw_fixed_list(pl, key, r, n, sep);
}

int param_list_parse_double_list(param_list_ptr pl, const char * key,
    double * r, size_t n, const char * sep)
{
    return param_list_parse_raw_fixed_list(pl, key, r, n, sep);
}

int param_list_parse_uint64_list(param_list_ptr pl, const char * key,
    uint64_t * r, size_t n, const char * sep)
{
    return param_list_parse_raw_fixed_list(pl, key, r, n, sep);
}

int param_list_parse_uint_list_size(param_list_ptr pl, const char * key , unsigned int ** r ,
        size_t *n)
{
    return param_list_parse_raw_variable_list (pl, key, r, n, ",");
}

const char * param_list_lookup_string(param_list_ptr pl, const char * key)
{
    return get_assoc_ptr(pl, key);
}

/**********************************************/
int param_list_parse_switch(param_list_ptr pl, const char * key)
{
    std::string value;
    int seen;
    for(int i = 0 ; i < 2 && *key == '-' ; key++, i++) ;
    if (!get_assoc(pl, key, value, false, &seen))
        return 0;
    if (!value.empty())
        throw parameter_exception(fmt::format(FMT_STRING("Parse error: option {} accepts no argument\n"), key));
    return seen;
}

int param_list_all_consumed(param_list_srcptr pl, const char ** extraneous)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    for(auto const & p : pli.p) {
        if (!p.second.parsed) {
            if (extraneous)
                *extraneous = p.first.c_str();
            return 0;
        }
    }
    return 1;
}

int param_list_warn_unused(param_list_srcptr pl)
{
    int u = 0;
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    for(auto const & p : pli.p) {
        if (!p.second.parsed && p.second.origin != PARAMETER_FROM_FILE) {
            fprintf(stderr, "Warning: unused command-line parameter %s\n",
                    p.first.c_str());
            u++;
        }
    }
    return u;
}

void param_list_display(param_list_srcptr pl, FILE *f)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    for(auto const & p : pli.p)
        fprintf(f,"%s=%s\n", p.first.c_str(), p.second.value.c_str());
}

void param_list_save(param_list_srcptr pl, const char * filename)
{
    FILE * f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }

    param_list_display(pl, f);
    fclose(f);
}

int param_list_save_parameter(param_list_ptr pl, enum parameter_origin o, 
        const char * key, const char * format, ...)
{
    va_list ap;
    va_start(ap, format);

    char * tmp;
    int rc;
    rc = vasprintf(&tmp, format, ap);
    param_list_add_key(pl, key, tmp, o);
    free(tmp);
    va_end(ap);

    return rc;
}


void param_list_print_command_line(FILE * stream, param_list_srcptr pl)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    if (!pli.cmdline_argv0)
        return;

    char **argv = pli.cmdline_argv0;
    int argc = pli.cmdline_argc0;

    if (verbose_enabled(CADO_VERBOSE_PRINT_CMDLINE)) {
        /* print command line */
        fprintf (stream, "# (%s) %s", cado_revision_string, argv[0]);
        for (int i = 1; i < argc; i++)
            fprintf (stream, " %s", argv[i]);
        fprintf (stream, "\n");
    }
    if (verbose_enabled(CADO_VERBOSE_PRINT_MODIFIED_FILES)) {
        if (strlen(cado_modified_files) > 1)
          fprintf (stream, "# List of modified files in working directory and "
                   "their SHA1 sum:\n%s", cado_modified_files);
    }
    if (verbose_enabled(CADO_VERBOSE_PRINT_COMPILATION_INFO)) {
#ifdef  __GNUC__
#ifndef __ICC
        fprintf (stream, "# Compiled with gcc " __VERSION__ "\n");
#else
        /* icc defines __GNUC__ too */
        fprintf (stream, "# Compiled with icc %d.%d.%d (gcc version %d.%d.%d compatibility)\n",
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
        fprintf (stream, "# WARNING: this version of GCC is known to miscompile CADO-NFS. See https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/14490\n");
#endif
#endif
        fprintf(stream, "# Compilation flags (C) " CFLAGS "\n");
        fprintf(stream, "# Compilation flags (C++) " CXXFLAGS "\n");
    }
}

void param_list_generic_failure(param_list_srcptr pl, const char *missing)
{
    param_list_impl const & pli = *static_cast<param_list_impl const *>(pl->pimpl);
    if (missing)
    {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, pli.cmdline_argv0[0], stderr);
    exit(EXIT_FAILURE);
}

/* This allows a mixed syntax for the specification of a list of values.
 * We accept either:
 *  - a comma-separated list of values, e.g., lpb=27,29
 *    (for historical reasons, lpbs=27,29 also works).
 *  - OR values being specified one by one, e.g., lpb0=27 lpb1=29
 *
 * This function definitely wants all resulting values to be set, so that
 * a policy has to be defined for the missing values (if the items are
 * specified in the "direct" way as with lpb1).
 *
 *  - ARGS_PER_SIDE_DEFAULT_AS_IS: we leave the missing values exactly as
 *    they were on input
 *  - ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS: if lpb<n> missing and lpb<n-1>
 *    is set, use the latter as a default value.
 *
 * RETURN VALUE: the number of sides for which data was left as is in the
 * input vector.
 */
template<typename T>
int param_list_parse_per_side(param_list_ptr pl, const char * key, T * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    for(int i = 0 ; i < 2 && *key == '-' ; key++, i++) ;
    ASSERT_ALWAYS(*key);
    int has_lpb01 = 0;
    for(int side = 0 ; side < n ; side++) {
        char * keyi;
        int rc = asprintf(&keyi, "%s%u", key, side);
        ASSERT_ALWAYS(rc >= 0);
        int gotit = param_list_parse_inner<T>(pl, keyi, lpb_arg[side], true);
        free(keyi);
        if (!gotit && side == 0 && policy == ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS)
            break;
        if (gotit)
            has_lpb01 = side + 1;
    }
    if (has_lpb01 && policy == ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS) {
        /* at this point we know that key0 has a value */
        for(int side = 0 ; side < n ; side++) {
            char * keyi;
            int rc = asprintf(&keyi, "%s%u", key, side);
            ASSERT_ALWAYS(rc >= 0);
            if (side)
                lpb_arg[side] = lpb_arg[side - 1];
            int gotit = param_list_parse_inner<T>(pl, keyi, lpb_arg[side], true);
            ASSERT_ALWAYS(side > 0 || gotit);
            free(keyi);
        }
    }

    int has_nlpbs;
    std::vector<T> temp;
    {

        if (!(has_nlpbs = param_list_parse_list(pl, key, temp, ","))) {
            char c = std::string(key).back();
            if (c != 'x' && c != 's') {
                /* try with s. */
                char * keys;
                int rc = asprintf(&keys, "%ss", key);
                ASSERT_ALWAYS(rc >= 0);
                has_nlpbs = param_list_parse_list(pl, keys, temp, ",", true);
                free(keys);
            }
        }
        has_nlpbs = temp.size();
        for(int i = 0 ; i < n && i < has_nlpbs ; i++)
            lpb_arg[i] = temp[i];
        if (has_nlpbs > n)
            throw parameter_exception(fmt::format(
                        FMT_STRING(
                            "Number of values for parameter {}"
                            " exceeds the maximum {}.")
                        , key, n));

    }

    if (has_nlpbs && has_lpb01) {
        fprintf(stderr, "Error, %1$s[01] and %1$s are incompatible\n", key);
        exit(EXIT_FAILURE);
    }

    if (!has_nlpbs && !has_lpb01) {
        return 0;
    }

    has_nlpbs = has_lpb01 = has_nlpbs + has_lpb01;

    if (policy == ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS) {
        /* Default to the last explicitly given lpb */
        for( ; has_nlpbs < n ; has_nlpbs++)
            lpb_arg[has_nlpbs] = lpb_arg[has_lpb01-1];

        if (has_nlpbs != n) {
            fprintf(stderr,
                    "Error, the number of values given for %ss does not "
                    "correspond to the number of polynomials\n", key);
        }
    }

    /* if the policy is to leave the unparsed parameters as is, then so
     * be it. Pretend we're happy, in that case. (it's quite likely that
     * we don't check the return value in this case anyway)
     */

    return n;
}

int param_list_parse_int_args_per_side(param_list_ptr pl, const char * key, int * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    return param_list_parse_per_side(pl, key, lpb_arg, n, policy);
}

int param_list_parse_uint_args_per_side(param_list_ptr pl, const char * key, unsigned int * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    return param_list_parse_per_side(pl, key, lpb_arg, n, policy);
}

/* We're only using these from C++ code at the moment, so that it is
 * sufficient to instantiate them here.
 */
template int param_list_parse_per_side<double>(param_list_ptr pl, const char * key, double * lpb_arg, int n, enum args_per_side_policy_t policy);
template int param_list_parse_per_side<int>(param_list_ptr pl, const char * key, int * lpb_arg, int n, enum args_per_side_policy_t policy);
template int param_list_parse_per_side<unsigned int>(param_list_ptr pl, const char * key, unsigned int * lpb_arg, int n, enum args_per_side_policy_t policy);
#ifndef UNSIGNED_LONG_IS_EXACTLY_UNSIGNED
template int param_list_parse_per_side<unsigned long>(param_list_ptr pl, const char * key, unsigned long * lpb_arg, int n, enum args_per_side_policy_t policy);
#endif
template int param_list_parse_per_side<std::string>(param_list_ptr pl, const char * key, std::string * lpb_arg, int n, enum args_per_side_policy_t policy);
