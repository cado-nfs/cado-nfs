/* dup2: 2nd pass

   Usage: dup2 -nrels <nrels> -renumber <renumberfile> [ -outdir <dir> ]
               [ -outfmt <fmt> ] [ -dl ]
               [ -filelist <fl> [ -basepath <dir> ] | file1 ... filen ]

   Input files can be given on command line, or via a filelist file.

   By default dup2 overwrites input files. To modify this behaviour an output
   directory can be given with '-outdir'.
   In case a filelist is given, the -basepath option enables to tell in which
   directory those files are.

   By default, the output will be bzipped/gzipped according to the status
   of the input file. The output format can be chosen with '-outfmt'.

   Allocates a hash-table of size next_prime(100 + 1.2*nrels).

   The relations are renumbered according to the file given via the
   '-renumberfile' argument. Input relations are of the following format:

       a,b[@s1,s2]:p1,p2,...,pj:q1,q2,...,qk                         (*)

   where p1,p2,...,pj are ideals on side s1 (possibly duplicate), and
   q1,q2,...,qk are ideals on side s2 (possibly duplicate). The part '@s1,s2'
   is optional and will default to s1=0 and s2=1 if not present. Output is:

       a,b:r1,r2,...,rm                                              (**)

   where each index r1,r2,...,rm refering to ideals via the renumbering table.
   By default valuations are reduced modulo 2. This can be changed by using the
   '-dl' argument.

   The format of each file is recognized by counting the number of ':' in the
   first line: if two we have the raw format (*), if only one we have the
   renumbered format (**). It is assumed that all files in renumbered format
   come first.

   Algorithm: for each (a,b) pair, we compute h(a,b) = (CA*a+CB*b) % 2^64.

   h has 64 bits. We know bits 2..6 are identical for a given slice, thus
   we remove them, and we obtain a 59-bit value h'.

   Let h' = j * 2^k + i.

   We store j % 2^32 at the next empty cell after index i in the hash-table.
*/

#include "cado.h"     // IWYU pragma: keep
#include <cerrno>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <memory>

#ifdef HAVE_MINGW
#include <fcntl.h>
#endif

#include "fmt/base.h"

#include "cado_poly.h"
#include "filter_config.h"
#include "filter_io.h"
#include "gmp_aux.h"
#include "gzip.h"
#include "macros.h"
#include "misc.h"
#include "params.h"
#include "portability.h"
#include "relation-tools.h"
#include "renumber.hpp"
#include "typedefs.h"
#include "verbose.h"
#include "portability.h"   // strdup       // IWYU pragma: keep

#define DEBUG 0

/* TODO: get rid of this gazillion globals! */

char const * argv0; /* = argv[0] */

/* Renumbering table to convert from (p,r) to an index */
static renumber_t renumber_tab;

static uint32_t * H;   /* H contains the hash table */
static uint64_t K = 0; /* Size of the hash table */
static unsigned long nrels_expected = 0;
static double cost = 0.0; /* Cost to insert all rels in the hash table */
/* Number of duplicates and rels on the current file */
static uint64_t ndup, nrels;
/* Number of duplicates and rels on all read files */
static uint64_t ndup_tot = 0, nrels_tot = 0;

/* sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
   and check for hash collisions */
std::vector<std::pair<cxx_mpz, cxx_mpz>> sanity_ab;
unsigned long sanity_checked = 0;
unsigned long sanity_collisions = 0;
/* end sanity check */

static double factor = 1.0;

static int is_for_dl; /* Do we reduce mod 2 or not */

/* The number of threads used to compute roots mod p is hard-coded to 3,
   since it was observed that using a large value gives some contention.
   The value of 3 is optimal for a c130 on a 64-core node. */
static int nthreads_for_roots = 3;

/* For debugging */
// #define TRACE_HASH_TABLE
#ifdef TRACE_HASH_TABLE
#define TRACE_I 42
#define TRACE_J 17
#endif

template <typename Ta, typename Tb>
static inline void sanity_check(uint64_t i, Ta a, Tb b)
{
    sanity_checked++;
    if (sanity_ab[i].first == 0) {
        sanity_ab[i].first = a;
        sanity_ab[i].second = b;
    } else if (sanity_ab[i].first != a) {
        sanity_collisions++;
        fmt::print(stderr,
                "Collision between ({}, {}) and ({}, {})\n",
                sanity_ab[i].first, sanity_ab[i].second,
                cxx_mpz(a), cxx_mpz(b));
    }
}

static inline void print_warning_size()
{
    const uint64_t nodup = nrels_tot - ndup_tot;
    const double full_table = 100.0 * double_ratio(nodup, K);
    fprintf(stderr, "Warning, hash table is %1.0f%% full (avg cost %1.2f)\n",
            full_table, double_ratio(cost, nrels_tot));
    if (full_table >= 99.0) {
        fprintf(stderr, "Error, hash table is full\n");
        exit(1);
    }
    factor += 1.0;
}

/* Print the relation 'rel' in a line of the form:
    a,b:h_1,h_2,...,h_k
   with a (signed) and b (unsigned) written in hexa and
   and i_1 ... i_k (hexadecimal) are the indices of the ideals

    The function adds a column of 1 if necessary, which is always
    column 0.
*/
template<filter_io_config cfg>
static inline void print_relation(FILE * file, typename cfg::rel_t rel)
{
    char buf[1 << 12], *p, *op;
    uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */

    p = d64toa16(buf, rel->a);
    *p++ = ',';
    p = u64toa16(p, rel->b);
    *p++ = ':';

    /* write everything with trailing commas, and strip the last one in
     * the end.  */
    for (unsigned int i = 0; i < rel->nb; i++) {
        if (rel->primes[i].e > 0) {
            op = p;
            p = u64toa16(p, (uint64_t)rel->primes[i].h);
            *p++ = ',';
            const ptrdiff_t t = p - op;
            for (int j = 1 ; j < rel->primes[i].e ; j++, p += t)
                memcpy(p, op, t);
            nonvoidside |= ((uint64_t)1) << rel->primes[i].side;
        }
    }

    /* If needed, print the additional columns (they are always at the beginning
     * of the renumbering table).
     * if naddcols == 0:
     *    do nothing.
     * else:
     *    we add the columns i if and only if the polynomial on side i is non
     *    monic and the relation contains at least one prime on side i.
     *
     * nb_polys==2 is special (see also  renumber_t::index_from_p_r) ; in
     * that case, we only use a single combined additional column. (except
     * for free relations, of course).
     */

    if (renumber_tab.number_of_additional_columns()) {
        size_t n = renumber_tab.get_nb_polys();
        if (n == 2) {
            /* Possible cases when we have two sides:
             *  - both are monic, there is no J ideal
             *  - only one is monic, there is only one J ideal
             *  - neither is monic, we have two J ideals, but since n==2
             *  it's okay to have only one additional column to count both
             *  of them together.
             * Either way, we only have one additional column (if we have
             * any, of course), and it's number zero
             */
            p = u64toa16(p, (uint64_t)0);
            *p++ = ',';
        } else {
            auto sides = renumber_tab.get_sides_of_additional_columns();
            for (size_t idx = 0; idx < sides.size(); idx++) {
                int side = sides[idx];
                if ((nonvoidside & (((uint64_t)1) << side))) {
                    p = u64toa16(p, (uint64_t)idx);
                    *p++ = ',';
                }
            }
        }
    }

    *(--p) = '\n';
    p[1] = '\0';
    if (fputs(buf, file) == EOF) {
        perror("Error writing relation");
        abort();
    }
}

static inline uint64_t
compute_hash(int64_t a, uint64_t b)
{
    return CA_DUP2 * (uint64_t)a + CB_DUP2 * b;
}

static inline uint64_t
compute_hash(mpz_srcptr a, mpz_srcptr b)
{
    cxx_mpz r;
    mpz_mul_uint64(r, a, CA_DUP2);
    mpz_addmul_uint64(r, b, CB_DUP2);
    return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
}

/* if duplicate is_dup = 1, else is_dup = 0
 * return i for sanity check, which is the place in the hash table where
 * the relation is stored: more precisely we store in H[i] the value
 * of floor(h/2^32), where h(a,b) is a 64-bit value.
 */
template<filter_io_config cfg>
static inline uint64_t
insert_relation_in_dup_hashtable(
        typename cfg::rel_srcptr rel,
        unsigned int * is_dup)
{
    uint64_t h, i;
    uint32_t j;
    double local_cost = 0;

    h = compute_hash(rel->a, rel->b);

    /* We put j = floor(h/2^32) in cell H[i] where i = h % K (or the first
       available cell after i if H[i] is already occupied). We have extraneous
       duplicates when:
       (a) either we have a collision on h: this should happen with probability
           2^(-64)
       (b) j = j' and |i-i'| is small, in particular i=i'. If K is at least
       twice the number of relations, then |i-i'| is bounded by say 2 on
       average, thus this happens when h' = h + n*K or h' = h + n*K + 1. This
       happens with probability 2/K. Moreover we should have j = j', which
       happens with probability 2^(-32). Since K is an odd prime, the global
           probability is about 2^(-31)/K. */

    i = h % K;
    j = (uint32_t)(h >> 32);
#ifdef TRACE_HASH_TABLE
    uint64_t old_i = i;
#endif
    while (H[i] != 0 && H[i] != j) {
        i++;
        if (UNLIKELY(i == K))
            i = 0;
        local_cost++;
    }

    /* The largest block for a linear probing table of size m with l entries
       behaves like (log(m)-1.5*log(log(m)))/(beta-1-log(beta)) where beta =
       l/m. See B. Pittel, Linear probing: the probable largest search time
       grows logarithmically with the number of records, Journal of Algorithms
       8, number 2, 236-249, 1987. Here we have m = K and l <= nrels_expected.
       Since we take K >= 6/5 * nrels_expected, we have beta <= 5/6,
       the factor 1/(beta-1-log(beta)) is bounded by 64.
       For m=2e8, the largest block can have up to length 938,
       thus it is not surprising to have local_cost > 100. */

    cost += local_cost;

    /* Note: since we use 0 for uninitialized entries, entries with j=0
       will get always marked as 'duplicate' and be lost. */

    if (H[i] == j)
        *is_dup = 1;
    else {
        H[i] = j;
        *is_dup = 0;
    }

#ifdef TRACE_HASH_TABLE
    if (i == TRACE_I && j == TRACE_J) {
        fmt::print(stderr,
                "TRACE: a = {}\nTRACE: b = {}\nTRACE: i = {}\nTRACE: j = {}\n"
                "TRACE: initial value of i was {}\nTRACE: h = {}\n"
                "TRACE: is_dup = {}\n",
                cxx_mpz(rel->a), cxx_mpz(rel->b), i, j, old_i, h, *is_dup);
    }
#endif
    return i;
}

/* modify in place the relation rel to take into account:
 *  - the renumbering
 *  - the bad ideals
 */
template<filter_io_config cfg>
static inline void compute_index_rel(typename cfg::rel_t rel)
{
    unsigned int i;
    p_r_values_t r;
    prime_t * pr = rel->primes;
    const weight_t len = rel->nb; // rel->nb can be modified by bad ideals

    for (i = 0; i < len; i++) {
        if (pr[i].e > 0) {
            if (pr[i].side != renumber_tab.get_rational_side()) {
#if DEBUG >= 1
                // Check for this bug : [#15897] [las] output "ideals" that are
                // not prime
                if (!modul_isprime(&(pr[i].p))) {
                    fprintf(stderr,
                            "Error, relation with a=%" PRId64 " b=%" PRIu64 " "
                            "contains %" PRpr " which is not prime.\nRemove "
                            "this relation from the file and re-run dup2.\n",
                            rel->a, rel->b, pr[i].p);
                    abort();
                }
#endif
                r = (p_r_values_t)relation_compute_r(rel->a, rel->b, pr[i].p);
            } else
                r = 0; // on the rational side we need not compute r, which is m
                       // mod p.

            const p_r_values_t p = pr[i].p;
            const int side = pr[i].side;
            renumber_t const & R = renumber_tab;
            renumber_t::p_r_side const prside { p, r, side };
            if (R.is_bad(prside)) {
                auto [ first_index, exps ] = R.indices_from_p_a_b(
                        prside, pr[i].e, rel->a, rel->b);

                /* allocate room for (nb) more valuations */
                for (; rel->nb + exps.size() - 1 > rel->nb_alloc;) {
                    realloc_buffer_primes<cfg>(rel);
                    pr = rel->primes;
                }
                /* the first is put in place, while the other are put at the end
                 * of the relation. As a side-effect, the relations produced are
                 * unsorted. Anyway, given that we're mixing sides when
                 * renumbering, we're bound to do sorting downhill. */
                pr[i].h = first_index;
                pr[i].e = !is_for_dl ? exps[0] & 1 : exps[0];

                for (size_t n = 1; n < exps.size(); n++) {
                    pr[rel->nb].h = first_index + n;
                    pr[rel->nb].e = exps[n];
                    if (!is_for_dl) { /* Do we reduce mod 2 */
                        pr[rel->nb].e &= 1;
                    }
                    rel->nb++;
                }
            } else {
                pr[i].h = R.index_from_p_r(prside);
                int const f = R.inertia_from_p_r(prside);
                if (f > 1) {
                    /* XXX there's a catch here. A non-bad ideal can still have
                     * non-trivial inertia (say f=2), in which case we must
                     * divide the valuation (which comes from the norm) by the
                     * inertia in order to obtain the valuation at the prime
                     * ideal
                     */
                    ASSERT_ALWAYS(pr[i].e % f == 0);
                    pr[i].e /= f;
                }
            }
        }
        if (!is_for_dl) { /* Do we reduce mod 2 */
            /* XXX should we compress as well ? */
            pr[i].e &= 1;
        }
    }
}

/* return in *oname and *oname_tmp two file names for writing the output
 * of processing the given input file infilename. Both files are placed
 * in the directory outdir if not NULL, otherwise in the current
 * directory.  The parameter outfmt specifies the output file extension
 * and format (semantics are as for fopen_maybe_compressed).
 *
 * proper use requires that data be first written to the file whose name
 * is *oname_tmp, and later on upon successful completion, that file must
 * be renamed to *oname. Otherwise disaster may occur, as there is a slim
 * possibility that *oname == infilename on return.
 *
 * TODO: rewrite with std::string
 */
static void get_outfilename_from_infilename(char const * infilename,
                                            char const * outfmt,
                                            char const * outdir, char ** oname,
                                            char ** oname_tmp)
{
    char const * suffix_in;
    char const * suffix_out;
    get_suffix_from_filename(infilename, &suffix_in);
    suffix_out = outfmt ? outfmt : suffix_in;

    char * newname = strdup(infilename);
    ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
    newname[strlen(newname) - strlen(suffix_in)] = '\0';

#define chkrcp(x)                                                              \
    do {                                                                       \
        const int rc = x;                                                      \
        ASSERT_ALWAYS(rc >= 0);                                                \
    } while (0)
    if (outdir) {
        char const * basename = path_basename(newname);
        chkrcp(
            asprintf(oname_tmp, "%s/%s.tmp%s", outdir, basename, suffix_out));
        chkrcp(asprintf(oname, "%s/%s%s", outdir, basename, suffix_out));
    } else {
        chkrcp(asprintf(oname_tmp, "%s.tmp%s", newname, suffix_out));
        chkrcp(asprintf(oname, "%s%s", newname, suffix_out));
    }
#undef chkrcp

#if DEBUG >= 1
    fprintf(stderr,
            "DEBUG: Input file name: %s,\nDEBUG: temporary output file "
            "name: %s,\nDEBUG: final output file name: %s\n",
            infilename, *oname_tmp, *oname);
#endif
    free(newname);
}

static void dup_print_stat(char const * s, uint64_t nrels, uint64_t ndup)
{
    const uint64_t nrem = nrels - ndup;
    const double pdup = 100.0 * double_ratio(ndup, nrels);
    fprintf(stderr,
            "%s: nrels=%" PRIu64 " dup=%" PRIu64 " (%.2f%%) rem=%" PRIu64 "\n",
            s, nrels, ndup, pdup, nrem);
}

template<filter_io_config cfg>
static void * hash_renumbered_rels(void *, typename cfg::rel_ptr rel)
{
    unsigned int is_dup;

    nrels++;
    nrels_tot++;
    uint64_t i = insert_relation_in_dup_hashtable<cfg>(rel, &is_dup);

    static unsigned long count = 0;

    // They should be no duplicate in already renumbered file
    if (is_dup && count++ < 10) {
        fmt::print(stderr,
                "Warning, duplicate relation in already renumbered files:"
                "\na = {}\nb = {}\ni = {}\nj = {}\n"
                "This warning may be due to a collision on the hash function "
                "or to an actual duplicate\nrelation. If it appears often you "
                "should check the input set of relations.\n\n",
                cxx_mpz(rel->a), cxx_mpz(rel->b), i, H[i]);
    }

    if (i < sanity_ab.size()) {
        if constexpr (std::is_same_v<cfg, filter_io_large_ab_cfg>) {
            /* Needed to avoid that the compiler infer mpz_struct * for the
             * type of rel->a and rel->b, which is not compatible with
             * operator== of the cxx_mpz class used in sanity_check.
             */
            sanity_check<mpz_srcptr, mpz_srcptr>(i, rel->a, rel->b);
        } else {
            sanity_check(i, rel->a, rel->b);
        }
    }

    if (cost >= factor * (double)(nrels_tot - ndup_tot))
        print_warning_size();

    return NULL;
}

template<filter_io_config cfg>
static void * thread_dup2(void * context_data, typename cfg::rel_ptr rel)
{
    unsigned int is_dup;
    uint64_t i;
    FILE * output = (FILE *)context_data;
    nrels++;
    nrels_tot++;
    i = insert_relation_in_dup_hashtable<cfg>(rel, &is_dup);
    if (!is_dup) {
        if (i < sanity_ab.size())
        {
            if constexpr (std::is_same_v<cfg, filter_io_large_ab_cfg>) {
                sanity_check<mpz_srcptr, mpz_srcptr>(i, rel->a, rel->b);
            } else {
                sanity_check(i, rel->a, rel->b);
            }
        }
        if (cost >= factor * (double)(nrels_tot - ndup_tot))
            print_warning_size();

        print_relation<cfg>(output, rel);
    } else {
        ndup++;
        ndup_tot++;
    }

    return NULL;
}

template<filter_io_config cfg>
static void * thread_root(void *, typename cfg::rel_ptr rel)
{
    /* We used to reduce exponents here. However, it's not a good idea if
     * we want to get the valuations at bad ideals.
     * (maybe we could reduce exponents for primes which we know are
     * above the bad ideal bound...).
     *
     * Anyway. Reduction is now done at the _end_ of compute_index_rel.
     */
    compute_index_rel<cfg>(rel);

    return NULL;
}

int check_whether_file_is_renumbered(char const * filename)
{
    unsigned int count = 0;
    char s[1024];
    FILE * f_tmp = fopen_maybe_compressed(filename, "rb");

    if (!f_tmp) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        abort();
    }

    if (feof(f_tmp)) /* file is empty */
    {
        fclose_maybe_compressed(f_tmp, filename);
        return 1; /* an empty file might be considered as renumbered */
    }

    /* Look for first non-comment line */
    while (1) {
        const char * ret = fgets(s, 1024, f_tmp);
        if (ret == NULL) {
            /* fgets returns NULL when the end of file occurs or when there is
               an error */
            if (feof(f_tmp)) {
                fclose_maybe_compressed(f_tmp, filename);
                return 1;
            }
            fprintf(stderr, "Error while reading %s\n", filename);
            exit(1);
        }
        if (strlen(s) >= 1023) {
            fprintf(stderr, "Too long line while reading %s\n", filename);
            exit(1);
        }
        size_t i = 0;
        while (s[i] == ' ')
            i++;
        if (s[i] != '#')
            break;
    }
    for (unsigned int i = 0; i < strlen(s); i++)
        count += s[i] == ':';
    fclose_maybe_compressed(f_tmp, filename);

    if (count == 1)
        return 1;
    else if (count == 2)
        return 0;
    else {
      fprintf (stderr, "Error: invalid line in %s (line has %u colons but 1 "
                       "or 2 were expected):\n %s", filename, count, s);
      abort();
    }
}

template<filter_io_config cfg>
static void
filter_new_rels(
        std::vector<std::string> const & files_new,
        char const * outdir,
        char const * outfmt)
{
    typename cfg::description_t desc[3] = {
        { .f = thread_root<cfg>, .arg = nullptr, .n = nthreads_for_roots, },
        { .f = thread_dup2<cfg>, .arg = nullptr, .n = 1, },
        { .f = nullptr,          .arg = nullptr, .n = 1, },
    };
    fmt::print(stderr, "Reading new files (using {} auxiliary threads for "
                       "roots mod p):\n", desc[0].n);

    /* TODO: use std::string for filename mangling */
    for (std::string const & p: files_new) {
        FILE * output = NULL;
        char *oname, *oname_tmp;

        get_outfilename_from_infilename(p.c_str(), outfmt, outdir, &oname,
                                        &oname_tmp);
        output = fopen_maybe_compressed(oname_tmp, "w");
        if (output == NULL) {
            fprintf(stderr,
                    "Error, could not open file to write the "
                    "relations. Check that the directory %s exists\n",
                    outdir);
            abort();
        }
        desc[1].arg = (void *)output;

        nrels = ndup = 0;

        uint64_t loc_nrels = filter_rels2<cfg>(
                std::vector<std::string> { p }, desc,
                EARLYPARSE_NEED_AB_DECIMAL | EARLYPARSE_NEED_PRIMES, NULL,
                NULL);

        ASSERT_ALWAYS(loc_nrels == nrels);

        fclose_maybe_compressed(output, oname_tmp);

        int rc;

#ifdef HAVE_MINGW /* For MinGW, rename cannot overwrite an existing file */
        rc = remove(oname);
#endif
        rc = rename(oname_tmp, oname);
        if (rc < 0) {
            fprintf(stderr, "rename(%s -> %s): %s\n", oname_tmp, oname,
                    strerror(errno));
            abort();
        }

        // stat for the current file
        dup_print_stat(path_basename(p.c_str()), nrels, ndup);
        // stat for all the files already read
        dup_print_stat("Total so far", nrels_tot, ndup_tot);

        free(oname);
        free(oname_tmp);
    }
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "input polynomial file");
    param_list_decl_usage(pl, "filelist",
                          "file containing a list of input files");
    param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
    param_list_decl_usage(pl, "renumber", "input file for renumbering table");
    param_list_decl_usage(pl, "nrels",
                          "number of relations to be found in the slice");
    param_list_decl_usage(pl, "outdir",
                          "by default, input files are overwritten");
    param_list_decl_usage(pl, "outfmt",
                          "format of output file (default same as input)");
    param_list_decl_usage(pl, "dl", "do not reduce exponents modulo 2");
    param_list_decl_usage(pl, "large-ab", "enable support for a and b larger "
                                          "than 64 bits");
    param_list_decl_usage(pl, "force-posix-threads",
                          "force the use of posix threads, do not rely on "
                          "platform memory semantics");
    param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
    verbose_decl_usage(pl);
}

static void usage(param_list pl, char const * argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int main(int argc, char const * argv[])
{
    argv0 = argv[0];

    cxx_param_list pl;
    cxx_cado_poly cpoly;

    declare_usage(pl);
    argv++, argc--;

    param_list_configure_switch(pl, "force-posix-threads",
                                &filter_rels_force_posix_threads);

    is_for_dl = 0; /* By default we do dup2 for factorization */
    param_list_configure_switch(pl, "dl", &is_for_dl);

    int largeab = 0; /* By default, do not use mpz for a,b */
    param_list_configure_switch(pl, "large-ab", &largeab);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY; /* Binary open for all files */
#endif

    if (argc == 0)
        usage(pl, argv0);

    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);

    char const * polyfilename = param_list_lookup_string(pl, "poly");
    char const * outfmt = param_list_lookup_string(pl, "outfmt");
    char const * filelist = param_list_lookup_string(pl, "filelist");
    char const * basepath = param_list_lookup_string(pl, "basepath");
    char const * outdir = param_list_lookup_string(pl, "outdir");
    char const * renumberfilename = param_list_lookup_string(pl, "renumber");
    char const * path_antebuffer =
        param_list_lookup_string(pl, "path_antebuffer");
    param_list_parse_ulong(pl, "nrels", &nrels_expected);

    if (param_list_warn_unused(pl)) {
        fprintf(stderr, "Error, unused parameters are given\n");
        usage(pl, argv0);
    }

    if (polyfilename == NULL || !cado_poly_read(cpoly, polyfilename)) {
        fprintf(stderr, "Error, missing -poly command line argument\n");
        usage(pl, argv0);
    }
    if (renumberfilename == NULL) {
        fprintf(stderr, "Error, missing -renumber command line argument\n");
        usage(pl, argv0);
    }
    if (nrels_expected == 0) {
        fprintf(stderr, "Error, missing -nrels command line argument "
                        "(or nrels = 0)\n");
        usage(pl, argv0);
    }
    K = 100 + nrels_expected + (nrels_expected / 5);
    K = uint64_nextprime(K);

    if (basepath && !filelist) {
        fprintf(stderr, "Error, -basepath only valid with -filelist\n");
        usage(pl, argv0);
    }
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "Error, output compression format unsupported\n");
        usage(pl, argv0);
    }
    if ((filelist != NULL) + (argc != 0) != 1) {
        fprintf(stderr,
                "Error, provide either -filelist or freeform file names\n");
        usage(pl, argv0);
    }

    set_antebuffer_path(argv0, path_antebuffer);

    renumber_tab = renumber_t(cpoly);
    renumber_tab.read_from_file(renumberfilename, is_for_dl);

    /* sanity check: since we allocate two 64-bit words for each, instead of
       one 32-bit word for the hash table, taking K/100 will use 2.5% extra
       memory */
    sanity_ab.resize(1u + (K / 100u));
    fmt::print(stderr, "[checking true duplicates on sample of {} cells]\n",
            sanity_ab.size());

    auto Hmem = std::make_unique<uint32_t[]>(K);
    H = Hmem.get();
    memset(H, 0, K * sizeof(uint32_t));
    fprintf(stderr,
            "Allocated hash table of %" PRIu64 " entries (%" PRIu64 "MiB)\n", K,
            (K * sizeof(uint32_t)) >> 20);

    /* Construct the two filelists : new files and already renumbered files */
    std::vector<std::string> files_already_renumbered;
    std::vector<std::string> files_new;
    {
        unsigned int nb_files = 0;
        fprintf(stderr, "Constructing the two filelists...\n");
        char const ** files =
            filelist ? filelist_from_file(basepath, filelist, 0) : argv;
        for (const auto * p = files; *p; p++)
            nb_files++;

        /* separate already processed files
         * check if f_tmp is in raw format a,b:...:... or
         *            in renumbered format a,b:...
         */
        for (const auto * p = files; *p; p++) {
            /* always strdup these, so that we can safely call
             * filelist_clear in the end */
            if (check_whether_file_is_renumbered(*p)) {
                files_already_renumbered.emplace_back(*p);
            } else {
                files_new.emplace_back(*p);
            }
        }
        fmt::print(stderr, "{} files ({} new and {} already renumbered)\n",
                nb_files, files_new.size(), files_already_renumbered.size());
        ASSERT_ALWAYS(nb_files == files_new.size() + files_already_renumbered.size());
        /* if filelist was not given, then files == argv, which of course
         * must not be cleared */
        if (filelist)
            filelist_clear(files);
    }

    fprintf(stderr, "Reading files already renumbered:\n");
    if (!largeab) {
        filter_rels<filter_io_default_cfg>(
                files_already_renumbered,
                &hash_renumbered_rels<filter_io_default_cfg>, NULL,
                EARLYPARSE_NEED_AB_HEXA, NULL, NULL);
    } else {
        filter_rels<filter_io_large_ab_cfg>(
                files_already_renumbered,
                &hash_renumbered_rels<filter_io_large_ab_cfg>, NULL,
                EARLYPARSE_NEED_AB_HEXA, NULL, NULL);
    }

    if (!largeab) {
        filter_new_rels<filter_io_default_cfg>(files_new, outdir, outfmt);
    } else {
        filter_new_rels<filter_io_large_ab_cfg>(files_new, outdir, outfmt);
    }

    fprintf(stderr, "At the end: %" PRIu64 " remaining relations\n",
            nrels_tot - ndup_tot);

    fprintf(stderr,
            "At the end: hash table is %1.2f%% full\n"
            "            hash table cost: %1.2f per relation\n",
            100.0 * double_ratio(nrels_tot - ndup_tot, K),
            1.0 + double_ratio(cost, nrels_tot));
    fprintf(stderr,
            "  [found %lu true duplicates on sample of %lu relations]\n",
            sanity_collisions, sanity_checked);

    if (files_already_renumbered.empty()) {
        if (nrels_tot != nrels_expected) {
            fprintf(stderr,
                    "Warning: number of relations read (%" PRIu64
                    ") does not match the number of relations expected (%lu)\n",
                    nrels_tot, nrels_expected);
        }
    } else {
        /* when we have renumbered files, we know that we won't have the
         * total number of relations... */
        if (nrels_tot > nrels_expected) {
            fprintf(stderr,
                    "Warning: number of relations read (%" PRIu64
                    ") exceeds the number of relations expected (%lu)\n",
                    nrels_tot, nrels_expected);
        }
    }

    return 0;
}
