#ifndef PARAMS_H_
#define PARAMS_H_

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

#include "macros.h"
#include "mpz_poly.h" // TODO: modify this.

struct param_list_doc_s {
    char * key;
    char * doc;
};
typedef struct param_list_doc_s param_list_doc[1];
typedef struct param_list_doc_s * param_list_doc_ptr;
typedef struct param_list_doc_s const * param_list_doc_srcptr;

/* This is by increasing order of priority */
enum parameter_origin { PARAMETER_FROM_FILE, PARAMETER_FROM_CMDLINE };

struct parameter_s {
    char * key;
    char * value;
    enum parameter_origin origin;
    int parsed;
    int seen;
};
typedef struct parameter_s parameter[1];
typedef struct parameter_s * parameter_ptr;
typedef struct parameter_s const * parameter_srcptr;

struct param_list_alias_s {
    char * alias;
    const char * key;
};

typedef struct param_list_alias_s param_list_alias[1];

struct param_list_switch_s {
    char * switchname;
    int * ptr;
};

typedef struct param_list_switch_s param_list_switch[1];

struct param_list_s {
    // documented parameters
    char * usage_hdr;
    int ndocs;
    int ndocs_alloc;
    param_list_doc * docs;
    // parameters given by user
    unsigned int alloc;
    unsigned int size;
    parameter * p;
    int consolidated;
    // aliases
    param_list_alias * aliases;
    int naliases;
    int naliases_alloc;
    // switches
    param_list_switch * switches;
    int nswitches;
    int nswitches_alloc;
    /* We use this to remember the first command line pointer which have
     * been given to us */
    int cmdline_argc0;
    char ** cmdline_argv0;
    // did the user use the doc functionality ?
    int use_doc;
};

typedef struct param_list_s param_list[1];
typedef struct param_list_s * param_list_ptr;
typedef struct param_list_s const * param_list_srcptr;

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
extern void param_list_print_usage(param_list_ptr pl, const char * argv0, FILE *f);
extern void param_list_usage_header(param_list_ptr pl, const char * hdr, ...)
    ATTR_PRINTF(2, 3);

// takes a file, in the Cado-NFS params format, and stores the dictionary
// of parameters to pl.
extern int param_list_read_stream(param_list_ptr pl, FILE *f, int stop_on_empty_line);
extern int param_list_read_file(param_list_ptr pl, const char * name);

// sees whether the arguments pointed to by argv[0] and (possibly)
// argv[1] correspond to either -<key> <value>, --<key> <value> or
// <key>=<value> ; configured switches and aliases for the param list are
// also checked.
extern int param_list_update_cmdline(param_list_ptr pl,
        int * p_argc, char *** p_argv);

#if 0
// This one allows shorthands. Notice that the alias string has to
// contain the exact form of the wanted alias, which may be either "-x",
// "--x", or "x=" (a terminating = tells the program that the option is
// wanted all in one go, like in ./a.out m=42, in contrast to ./a.out -m
// 42).
extern int param_list_update_cmdline_alias(param_list_ptr pl, const char * key,
        const char * alias, int * p_argc, char *** p_argv);
#endif

extern int param_list_parse_int(param_list_ptr, const char *, int *);
extern int param_list_parse_long(param_list_ptr, const char *, long *);
extern int param_list_parse_uint(param_list_ptr, const char *, unsigned int *);
extern int param_list_parse_ulong(param_list_ptr, const char *, unsigned long *);
extern int param_list_parse_int64(param_list_ptr, const char *, int64_t *);
extern int param_list_parse_uint64(param_list_ptr, const char *, uint64_t *);
extern int param_list_parse_uchar(param_list_ptr, const char * , unsigned char *);
extern int param_list_parse_double(param_list_ptr, const char *, double *);
extern int param_list_parse_double_and_double(param_list_ptr, const char *,
    double *, const char *);
extern int param_list_parse_string(param_list_ptr, const char *, char *, size_t);
extern int param_list_parse_mpz(param_list_ptr, const char *, mpz_ptr);
extern int param_list_parse_intxint(param_list_ptr, const char * key, int * r);
extern int param_list_parse_int_and_int(param_list_ptr, const char * key, int * r, const char * sep);
int param_list_parse_uint_and_uint(param_list_ptr, const char * key, unsigned int * r, const char * sep);
extern int param_list_parse_long_and_long(param_list_ptr, const char * key, long * r, const char * sep);
extern int param_list_parse_string_list_alloc(param_list_ptr, const char * key, char *** r, int * n, const char * sep);
extern int param_list_get_list_count(param_list_ptr pl, const char * key);
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
  Return an array r with its size t. The array is initialised with the string
  separate by sep.
  Usage: "5,5,5" gives r = [5,5,5] and t = 3.

  pl: parameter list.
  key: key in the parameter list.
  r: array with the integer value contained in the key.
  t: number of element in array.
*/
extern int param_list_parse_int_list_size(param_list_ptr, const char * key,
    int ** r, unsigned int * t);
/*
  Return a mpz_poly f. The polynomial is initialised with the string separate by
  sep.
  Usage: "5,4,3" gives f = 5 + 4 * x + 3 * x^2.

  pl: parameter list.
  key: key in the parameter list.
  f: the polynomial.
*/
extern int param_list_parse_mpz_poly(param_list_ptr, const char * key,
    mpz_poly_ptr f);

extern int param_list_parse_size_t(param_list_ptr, const char * key, size_t * r);
extern int param_list_parse_switch(param_list_ptr, const char * key);

extern const char * param_list_lookup_string(param_list_ptr, const char * key);

extern void param_list_save(param_list_ptr, const char * filename);

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
extern int param_list_all_consumed(param_list_ptr, char ** extraneous);

// warns against unused command-line parameters. This normally indicates
// a user error. parameters ignored from config files are considered
// normal (although note that in some cases, it could be bad as well).
extern int param_list_warn_unused(param_list pl);

// this one is the ``joker'' call. Return type is for internal use.
extern int param_list_add_key(param_list_ptr,
        const char *, const char *, enum parameter_origin);
// removing a key can be handy before savin g a config file. Some options
// are relevant only for one particular invokation, and not for saving.
extern void param_list_remove_key(param_list_ptr, const char * key);

// for debugging.
extern void param_list_display(param_list_ptr, FILE *f);

// quick way to reinject parameters in the param_list (presumably before
// saving)
extern int param_list_save_parameter(param_list_ptr, enum parameter_origin o, 
        const char * key, const char * format, ...) ATTR_PRINTF(4,5);

// This function is a shorthand which does employ some hackery put into
// param lists, which remember their oldest argv, argc pair.
extern void param_list_print_command_line(FILE * stream, param_list_ptr);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* Same idea as for cxx_mpz and friends */
struct cxx_param_list {
    param_list x;
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
#if __cplusplus >= 201103L
    cxx_param_list(cxx_param_list && o) {
        param_list_init(x);
        param_list_swap(x, o.x);
    }
    cxx_param_list& operator=(cxx_param_list && o) {
        param_list_swap(x, o.x);
        return *this;
    }
#endif
    operator param_list_ptr() { return x; }
    operator param_list_srcptr() const { return x; }
    param_list_ptr operator->() { return x; }
    param_list_srcptr operator->() const { return x; }
};

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void param_list_init(cxx_param_list & pl) __attribute__((error("param_list_init must not be called on a param_list reference -- it is the caller's business (via a ctor)")));
extern void param_list_clear(cxx_param_list & pl) __attribute__((error("param_list_clear must not be called on a param_list reference -- it is the caller's business (via a dtor)")));
#endif

#endif


#endif	/* PARAMS_H_ */
