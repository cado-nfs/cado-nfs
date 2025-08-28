#ifndef CADO_PARAMS_H
#define CADO_PARAMS_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
#include <istream>
#include <map>
#include <string>
#include <utility>
#include <stdexcept>
#endif

#include <gmp.h>

#include "macros.h"
#ifdef __cplusplus
#include "cxx_mpz.hpp"
#include "prime_power_factorization.hpp"
#endif
#include "mpz_poly.h" // TODO: modify this.

/* This is by increasing order of priority */
enum parameter_origin { PARAMETER_FROM_FILE, PARAMETER_FROM_CMDLINE };

struct param_list_s {
    /* it's not a real pimpl idiom, for the simple reason that we're also
     * called from C.
     */
    void * pimpl;
};
typedef struct param_list_s param_list[1];
typedef struct param_list_s * param_list_ptr;
typedef struct param_list_s const * param_list_srcptr;

enum args_per_side_policy_t {
    ARGS_PER_SIDE_DEFAULT_AS_IS,
    ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS,
};

#ifdef __cplusplus
extern "C" {
#endif

// in any case, calls to param_list functions overwrite the previously
// set parameters in the parameter list.

extern void param_list_init(param_list_ptr pl);
extern void param_list_clear(param_list_ptr pl);
extern void param_list_set(param_list_ptr pl, param_list_srcptr pl0);
extern void param_list_swap(param_list_ptr pl, param_list_ptr pl0);

// document the usage of a parameter.
extern void param_list_decl_usage(param_list_ptr pl, const char * key,
        const char * doc);
extern void param_list_print_usage(param_list_srcptr pl, const char * argv0, FILE *f);
extern void param_list_usage_header(param_list_ptr pl, const char * hdr);

// takes a file, in the Cado-NFS params format, and stores the dictionary
// of parameters to pl.
extern int param_list_read_stream(param_list_ptr pl, FILE *f, int stop_on_empty_line);
extern int param_list_read_file(param_list_ptr pl, const char * name);

// sees whether the arguments pointed to by argv[0] and (possibly)
// argv[1] correspond to either -<key> <value>, --<key> <value> or
// <key>=<value> ; configured switches and aliases for the param list are
// also checked.
extern int param_list_update_cmdline(param_list_ptr pl,
        int * p_argc, char const *** p_argv) ATTRIBUTE_NONNULL((2,3));

extern void param_list_generic_failure(param_list_srcptr pl, const char *missing);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
template<typename T>
int param_list_parse(param_list_ptr pl, std::string const & key, T & r);

/* this returns a default constructed T if the key is absent */
template<typename T>
T
param_list_parse(param_list_ptr pl, std::string const & key)
{
    T r;
    param_list_parse<T>(pl, key, r);
    return r;
}

template<typename T>
T
param_list_parse_mandatory(param_list_ptr pl, std::string const & key)
{
    T r;
    if (!param_list_parse<T>(pl, key, r))
        param_list_generic_failure(pl, key.c_str());

    return r;
}


template<typename T>
int param_list_parse_per_side(param_list_ptr pl, std::string const & key, T * lpb_arg, int n, enum args_per_side_policy_t policy);

extern template int param_list_parse_per_side<double>(param_list_ptr pl, std::string const & key, double * lpb_arg, int n, enum args_per_side_policy_t policy);
extern template int param_list_parse_per_side<int>(param_list_ptr pl, std::string const & key, int * lpb_arg, int n, enum args_per_side_policy_t policy);
extern template int param_list_parse_per_side<unsigned int>(param_list_ptr pl, std::string const & key, unsigned int * lpb_arg, int n, enum args_per_side_policy_t policy);
#ifndef UNSIGNED_LONG_IS_EXACTLY_UNSIGNED
extern template int param_list_parse_per_side<unsigned long>(param_list_ptr pl, std::string const & key, unsigned long * lpb_arg, int n, enum args_per_side_policy_t policy);
#endif
extern template int param_list_parse_per_side<std::string>(param_list_ptr pl, std::string const & key, std::string * lpb_arg, int n, enum args_per_side_policy_t policy);

/* We have all of these defined in params.cpp, and they can be used from
 * c++ code only.
 */
extern template int param_list_parse<bool>(param_list_ptr pl, std::string const & key, bool & r);
extern template int param_list_parse<int>(param_list_ptr pl, std::string const & key, int & r);
extern template int param_list_parse<unsigned int>(param_list_ptr pl, std::string const & key, unsigned int & r);
#ifndef LONG_IS_EXACTLY_INT
extern template int param_list_parse<long>(param_list_ptr pl, std::string const & key, long & r);
#endif
#ifndef UNSIGNED_LONG_IS_EXACTLY_UNSIGNED
extern template int param_list_parse<unsigned long>(param_list_ptr pl, std::string const & key, unsigned long & r);
#endif
#ifndef INT64_T_IS_EXACTLY_LONG
extern template int param_list_parse<int64_t>(param_list_ptr pl, std::string const & key, int64_t & r);
#endif
#ifndef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
extern template int param_list_parse<uint64_t>(param_list_ptr pl, std::string const & key, uint64_t & r);
#endif
extern template int param_list_parse<double>(param_list_ptr pl, std::string const & key, double & r);
extern template int param_list_parse<std::vector<int>>(param_list_ptr pl, std::string const & key, std::vector<int> & r);
extern template int param_list_parse<std::vector<unsigned int>>(param_list_ptr pl, std::string const & key, std::vector<unsigned int> & r);
extern template int param_list_parse<std::vector<std::string>>(param_list_ptr pl, std::string const & key, std::vector<std::string> & r);
extern template int param_list_parse<cxx_mpz>(param_list_ptr pl, std::string const & key, cxx_mpz & r);
extern template int param_list_parse<cxx_mpz_poly>(param_list_ptr pl, std::string const & key, cxx_mpz_poly & r);
extern template int param_list_parse<cado::prime_power_factorization>(param_list_ptr pl, std::string const & key, cado::prime_power_factorization & r);
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern int param_list_empty(param_list_srcptr);
extern int param_list_parse_int(param_list_ptr, const char *, int *);
extern int param_list_parse_long(param_list_ptr, const char *, long *);
extern int param_list_parse_uint(param_list_ptr, const char *, unsigned int *);
extern int param_list_parse_ulong(param_list_ptr, const char *, unsigned long *);
extern int param_list_parse_int64(param_list_ptr, const char *, int64_t *);
extern int param_list_parse_uint64(param_list_ptr, const char *, uint64_t *);
extern int param_list_parse_double(param_list_ptr, const char *, double *);
extern int param_list_parse_double_and_double(param_list_ptr, const char *,
    double *, const char *);
extern int param_list_parse_mpz(param_list_ptr, const char *, mpz_ptr);
extern int param_list_parse_intxint(param_list_ptr, const char * key, int * r);
extern int param_list_parse_int_and_int(param_list_ptr, const char * key, int * r, const char * sep);
int param_list_parse_uint_and_uint(param_list_ptr, const char * key, unsigned int * r, const char * sep);
extern int param_list_parse_long_and_long(param_list_ptr, const char * key, long * r, const char * sep);
extern int param_list_parse_ulong_and_ulong(param_list_ptr pl, const char * key, unsigned long * r, const char * sep);
extern size_t param_list_get_list_count(param_list_ptr pl, const char * key);
extern int param_list_parse_int_list(param_list_ptr, const char * key, int * r, size_t n, const char * sep);
int param_list_parse_uint64_and_uint64(param_list_ptr, const char * key,
    uint64_t * r, const char * sep);
extern int param_list_parse_uint_list(param_list_ptr, const char * key,
    unsigned int * r, size_t n, const char * sep);
int param_list_parse_uint64_list(param_list_ptr, const char * key,
    uint64_t * r, size_t n, const char * sep);
extern int param_list_parse_uchar_list(param_list_ptr, const char * key,
    unsigned char * r, size_t n, const char * sep);
int param_list_parse_double_list(param_list_ptr, const char * key,
    double * r, size_t n, const char * sep);

/*
 * Return an mpz_poly f. Can be given either as a comma-separated list of
 * coefficients, or as a polynomial expression.
 *
 * pl: parameter list.
 * key: key in the parameter list.
 * f: the polynomial.
 */
extern int param_list_parse_mpz_poly(param_list_ptr, const char * key,
    mpz_poly_ptr f);

extern int param_list_parse_size_t(param_list_ptr, const char * key, size_t * r);
extern int param_list_parse_switch(param_list_ptr, const char * key);

extern const char * param_list_lookup_string(param_list_ptr, const char * key);

// This one allows shorthands. Notice that the alias string has to
// contain the exact form of the wanted alias, which may be either "-x",
// "--x", or "x=" (a terminating = tells the program that the option is
// wanted all in one go, like in ./a.out m=42, in contrast to ./a.out -m
// 42).
extern int param_list_configure_alias(param_list_ptr, const char * key, const char * alias);

// A switch is a command-line argument which sets a value by its mere
// presence. Could be for instance --verbose, or --use-smart-algorithm
extern int param_list_configure_switch(param_list_ptr, const char * key, int * ptr);

// tells whether everything has been consumed. Otherwise, return the key
// of the first unconsumed argument.
extern int param_list_all_consumed(param_list_srcptr, const char ** extraneous);

// warns against unused command-line parameters. This normally indicates
// a user error. parameters ignored from config files are considered
// normal (although note that in some cases, it could be bad as well).
extern int param_list_warn_unused(param_list_srcptr pl);

// this one is the ``joker'' call. Return type is for internal use.
extern void param_list_add_key(param_list_ptr,
        const char *, const char *, enum parameter_origin);
// removing a key can be handy before savin g a config file. Some options
// are relevant only for one particular invokation, and not for saving.
extern void param_list_remove_key(param_list_ptr, const char * key);

// for debugging.
extern void param_list_display(param_list_srcptr, FILE *f);

#if 0
/* these two are never used */
extern void param_list_save(param_list_srcptr, const char * filename);

// quick way to reinject parameters in the param_list (presumably before
// saving)
extern int param_list_save_parameter(param_list_ptr, enum parameter_origin o, 
        const char * key, const char * format, ...) ATTR_PRINTF(4,5);
#endif

// This function is a shorthand which does employ some hackery put into
// param lists, which remember their oldest argv, argc pair.
extern void param_list_print_command_line(FILE * stream, param_list_srcptr);

extern int param_list_parse_uint_args_per_side(param_list_ptr pl, const char * key, unsigned int * lpb_arg, int n, enum args_per_side_policy_t policy);
extern int param_list_parse_int_args_per_side(param_list_ptr pl, const char * key, int * lpb_arg, int n, enum args_per_side_policy_t policy);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern int param_list_read(param_list_ptr pl, std::istream & is, bool stop_on_empty_line = false);
#endif

#ifdef __cplusplus

/* The param_list implementation is only exposed to c++ code
 */
struct param_list_impl {
    // documented parameters
    std::string usage_header;
    struct collate {
        /* compare two strings, intentionally collating - and _ (except at the
         * beginning of the string) */
        bool operator()(std::string const & a, std::string const & b) const {
            size_t k;
            for(k = 0 ; k < a.size() && k < b.size() ; k++) {
                int r = (a[k] > b[k]) - (b[k] > a[k]);
                if (k && (a[k] == '-' || a[k] == '_') && (b[k] == '-' || b[k] == '_')) r = 0;
                if (r) return r < 0;
            }
            return (a[k] > b[k]) < (b[k] > a[k]);
        }
    };
    std::map<std::string, std::string, collate> documentation;   /* for each key */
    struct parameter {
        std::string value;
        enum parameter_origin origin;
        bool parsed;
        int seen;
        explicit parameter(std::string value = std::string(), 
                enum parameter_origin origin = PARAMETER_FROM_FILE,
                bool parsed = false,
                int seen = 1)
            : value(std::move(value))
            , origin(origin)
            , parsed(parsed)
            , seen(seen) {}
        parameter(parameter const &) = default;
        parameter& operator=(parameter const &) = default;
        parameter(parameter &&) = default;
        parameter& operator=(parameter &&) = default;
        ~parameter() = default;
    };
    std::map<std::string, parameter, collate> p;
    // aliases
    std::map<std::string, std::string, collate> aliases;
    // switches
    std::map<std::string, int *, collate> switches;
    /* We use this to remember the first command line pointer which have
     * been given to us */
    int cmdline_argc0;
    char const ** cmdline_argv0;
    // did the user use the doc functionality ?
    bool use_doc;
};
#endif

#ifdef __cplusplus
/* Same idea as for cxx_mpz and friends */
struct cxx_param_list {
    param_list x;
    param_list_impl& impl() { return *static_cast<param_list_impl *>(x->pimpl); }
    param_list_impl const & impl() const { return *static_cast<param_list_impl const *>(x->pimpl); }

    cxx_param_list() { param_list_init(x); }
    cxx_param_list(param_list_srcptr f) { param_list_init(x); param_list_set(x, f); }
    ~cxx_param_list() { param_list_clear(x); }
    cxx_param_list(cxx_param_list const & o) {
        param_list_init(x);
        param_list_set(x, o.x);
    }
    cxx_param_list & operator=(cxx_param_list const & o) {
        param_list_set(x, o.x);
        return *this;
    }
    cxx_param_list(cxx_param_list && o) noexcept {
        param_list_init(x);
        param_list_swap(x, o.x);
    }
    cxx_param_list& operator=(cxx_param_list && o) noexcept {
        param_list_swap(x, o.x);
        return *this;
    }
    operator param_list_ptr() { return x; }
    operator param_list_srcptr() const { return x; }
    param_list_ptr operator->() { return x; }
    param_list_srcptr operator->() const { return x; }
};

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void param_list_init(cxx_param_list & pl) __attribute__((error("param_list_init must not be called on a param_list reference -- it is the caller's business (via a ctor)")));
extern void param_list_clear(cxx_param_list & pl) __attribute__((error("param_list_clear must not be called on a param_list reference -- it is the caller's business (via a dtor)")));
#endif

struct parameter_error : public std::runtime_error {
    explicit parameter_error(std::string const & arg)
        : std::runtime_error(arg)
    {}
};

#endif


#endif	/* CADO_PARAMS_H */
