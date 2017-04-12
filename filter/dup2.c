/* dup2: 2nd pass

   Usage: dup2 -nrels <nrels> -renumber <renumberfile> [ -outdir <dir> ]
               [ -outfmt <fmt> ] [ -dl ] [ -badidealinfo <file> ]
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

       a,b:p1,p2,...,pj:q1,q2,...,qk                                 (*)

   where p1,p2,...,pj are rational (side 0) ideals (possibly duplicate), and
   q1,q2,...,qk are algebraic (side 1) ideals (possibly duplicate). Output is:

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

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "filter_config.h"
#include "utils_with_io.h"

#include "filter_badideals.h"
#include "relation-tools.h"

#define DEBUG 0

static const uint64_t constants_ab_dup2[EARLYPARSED_RELATION_MAX_AB] = {
    /* the first two used to be called CA_DUP2 and CB_DUP2 -- well,
     * except that now that we've negated the semantics of b, the hash
     * values differ anyway... */
    UINT64_C(271828182845904523),
    UINT64_C(577215664901532889),
    /* These are just obtained with RandomPrime(64) */
    UINT64_C(11362315840839667127),
    UINT64_C(5093043268839177131),
};



char *argv0; /* = argv[0] */

/* Renumbering table to convert from (p,r) to an index */
renumber_t renumber_tab;

static uint32_t *H; /* H contains the hash table */
static uint64_t K = 0; /* Size of the hash table */
static unsigned long nrels_expected = 0;
static double cost = 0.0; /* Cost to insert all rels in the hash table */
/* Number of duplicates and rels on the current file */
static uint64_t ndup, nrels;
/* Number of duplicates and rels on all read files */
static uint64_t ndup_tot = 0, nrels_tot = 0;

/* sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
   and check for hash collisions */
uint64_t sanity_size;
int64_t *sanity_a;
int64_t *sanity_b;
unsigned long sanity_checked = 0;
unsigned long sanity_collisions = 0;
/* end sanity check */

static double factor = 1.0;

static int is_for_dl; /* Do we reduce mod 2 or not */

/* For debugging */
//#define TRACE_HASH_TABLE
#ifdef TRACE_HASH_TABLE
#define TRACE_I 42
#define TRACE_J 17
#endif

static inline void
sanity_check (uint64_t i, int64_t a, int64_t b)
{
  sanity_checked++;
  if (sanity_a[i] == 0)
  {
    sanity_a[i] = a;
    sanity_b[i] = b;
  }
  else if (sanity_a[i] != a)
  {
    sanity_collisions++;
    fprintf(stderr, "Collision between (%" PRId64 ",%" PRIu64 ") and "
                    "(%" PRId64 ",%" PRIu64 ")\n", sanity_a[i], sanity_b[i],
                    a, b);
  }
}

static inline void
print_warning_size ()
{
  uint64_t nodup = nrels_tot - ndup_tot;
  double full_table = 100.0 * (double) nodup / (double) K;
  fprintf (stderr, "Warning, hash table is %1.0f%% full (avg cost %1.2f)\n",
           full_table, cost / (double) nrels_tot);
  if (full_table >= 99.0)
  {
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
static inline void
print_relation (FILE * file, earlyparsed_relation_srcptr rel)
{
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;
  uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */

#if 0
  p = d64toa16(buf, rel->a);
  *p++ = ',';
  p = u64toa16(p, rel->b);
  *p++ = ':';
#else
  int nab = EARLYPARSED_RELATION_MAX_AB;
  for( ; nab > 2 && !rel->ab[nab-1] ; nab--);
  p = buf;
  *p++ = 'X';
  for(int i = 0 ; i < nab ; i++) {
      *p++ = i ? ',' : ' ';
      p = d64toa16(p, rel->ab[i]);
  }
  *p++ = ':';
#endif

  for (i = 0; i < rel->nb; i++)
  {
    if (rel->primes[i].e > 0)
    {
      op = p;
      p = u64toa16(p, (uint64_t) rel->primes[i].h);
      *p++ = ',';
      t = p - op;
      for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
        memcpy(p, op, t);
      nonvoidside |= ((uint64_t) 1) << rel->primes[i].side;
    }
  }

  /* If needed, print the additional columns (they are always at the beginning
   * of the renumbering table).
   * if naddcols == 0:
   *    do nothing.
   * else: 
   *    if nb_polys == 2:
   *      we add the columns 0 (in this case there is always 1 additional column
   *      and it is always necessary)
   *    if nb_polys != 2:
   *      we add the columns i if and only if the ith bit of
   *      renumber_tab->nonmonic is 1 (i.e., the polynomial on side i is non
   *      monic) and the relation contains at least one prime on side i.
   */
  if (renumber_tab->naddcols)
  {
    if (renumber_tab->nb_polys == 2)
    {
      p = u64toa16(p, (uint64_t) 0);
      *p++ = ',';
    }
    else
    {
      index_t index_add_col = 0;
      for (uint64_t b = renumber_tab->nonmonic, side = 0; b != 0;
                                                b>>=1, nonvoidside>>=1, side++)
      {
        if (b & ((uint64_t) 1))
        {
          if (nonvoidside & ((uint64_t) 1))
          {
            p = u64toa16(p, (uint64_t) index_add_col);
            *p++ = ',';
          }
          index_add_col++;
        }
      }
    }
  }


  *(--p) = '\n';
  p[1] = 0;
  if (fputs(buf, file) == EOF) {
    perror("Error writing relation");
    abort();
  }
}

/* if duplicate is_dup = 1, else is_dup = 0
 * return i for sanity check, which is the place in the hash table where
 * the relation is stored: more precisely we store in H[i] the value
 * of floor(h/2^32), where h(a,b) is a 64-bit value.
 */
static inline uint64_t
insert_relation_in_dup_hashtable (earlyparsed_relation_srcptr rel, unsigned int *is_dup)
{
  uint64_t h=0, i;
  uint32_t j;
  double local_cost = 0;

  for(int i = 0 ; i < EARLYPARSED_RELATION_MAX_AB ; i++) {
      h += constants_ab_dup2[i] * rel->ab[i];
  }

  /* We put j = floor(h/2^32) in cell H[i] where i = h % K (or the first
     available cell after i if H[i] is already occupied). We have extraneous
     duplicates when:
     (a) either we have a collision on h: this should happen with probability
         2^(-64)
     (b) j = j' and |i-i'| is small, in particular i=i'. If K is at least twice
         the number of relations, then |i-i'| is bounded by say 2 on average,
         thus this happens when h' = h + n*K or h' = h + n*K + 1. This happens
         with probability 2/K. Moreover we should have j = j', which happens
         with probability 2^(-32). Since K is an odd prime, the global
         probability is about 2^(-31)/K. */

  i = h % K;
  j = (uint32_t) (h >> 32);
#ifdef TRACE_HASH_TABLE
  uint64_t old_i = i;
#endif
  while (H[i] != 0 && H[i] != j)
  {
    i++;
    if (UNLIKELY(i == K))
      i = 0;
    local_cost++;
  }

  if (local_cost > 100)
  {
      static int count = 0;
      if (count++ < 10) {
          mpz_poly ab;
          mpz_poly_init(ab, EARLYPARSED_RELATION_MAX_AB - 1);
          mpz_poly_setcoeffs_int64(ab, rel->ab, EARLYPARSED_RELATION_MAX_AB - 1);
          fprintf (stderr, "Warning, insertion cost %1.0f for ab=", local_cost);
          mpz_poly_fprintf_coeffs(stderr, ab, ',');
          fprintf (stderr, " h=%" PRIu64 " i=%" PRIu64 " j=%u\n", h, i, j);
          mpz_poly_clear(ab);
      }
    }

  cost += local_cost;

  /* Note: since we use 0 for uninitialized entries, entries with j=0
     will get always marked as 'duplicate' and be lost. */

  if (H[i] == j)
    *is_dup = 1;
  else
  {
    H[i] = j;
    *is_dup = 0;
  }

#ifdef TRACE_HASH_TABLE
  if (i == TRACE_I && j == TRACE_J)
  {
    fprintf (stderr, "TRACE: a = %s%" PRIx64 "\nTRACE: b = %" PRIx64 "\n"
                     "TRACE: i = %" PRIu64 "\nTRACE: j = %" PRIu32 "\n"
                     "TRACE: initial value of i was %" PRIu64 "\n"
                     "TRACE: h = %" PRIx64 "\nTRACE: is_dup = %u\n",
                     (rel->a < 0) ? "-" : "",
                     (uint64_t) ((rel->a < 0) ? -rel->a : rel->a),
                     rel->b, i, j, old_i, h, *is_dup);
  }
#endif
  return i;
}

/* modify in place the relation rel to take into account:
 *  - the renumbering
 *  - the bad ideals
 */
static inline void
compute_index_rel (earlyparsed_relation_ptr rel, allbad_info_t info)
{
  unsigned int i;
  p_r_values_t r;
  prime_t *pr = rel->primes;
  weight_t len = rel->nb; // rel->nb can be modified by bad ideals

  int nab = EARLYPARSED_RELATION_MAX_AB;
  for( ; nab > 2 && !rel->ab[nab-1] ; nab--);
  ASSERT_ALWAYS(nab == 2);

  for (i = 0; i < len; i++)
  {
    if (pr[i].e > 0)
    {
      if (pr[i].side != renumber_tab->rat)
      {
#if DEBUG >= 1
  // Check for this bug : [#15897] [las] output "ideals" that are not prime
        if (!modul_isprime(&(pr[i].p)))
        {
          fprintf (stderr, "Error, relation with a=%" PRId64" b=%" PRIu64 " "
                           "contains %" PRpr " which is not prime.\nRemove "
                           "this relation from the file and re-run dup2.\n",
                           rel->a, rel->b, pr[i].p);
          abort();
        }
#endif
        r = (p_r_values_t) relation_compute_r (rel->ab[0], rel->ab[1], pr[i].p);
      }
      else
        r = 0; // on the rational side we need not compute r, which is m mod p.
      
      int nb; //number of ideals above the bad ideal
      index_t first_index; // first index of the ideals above a bad ideal
      if (pr[i].p <= renumber_tab->bad_ideals.max_p
          && renumber_is_bad (&nb, &first_index, renumber_tab, pr[i].p, r,
                                                                    pr[i].side))
      {
        int exp_above[RENUMBER_MAX_ABOVE_BADIDEALS] = {0,};
        handle_bad_ideals (exp_above, rel->ab[0], rel->ab[1], pr[i].p, pr[i].e,
                           pr[i].side, info);
        
        /* allocate room for (nb) more valuations */
        if (rel->nb + nb - 1 > rel->nb_alloc)
        {
           realloc_buffer_primes(rel);
           pr = rel->primes;
        }

        /* the first is put in place, while the other are put at the end
         * of the relation. As a side-effect, the relations produced are
         * unsorted. Anyway, given that we're mixing sides when
         * renumbering, we're bound to do sorting downhill. */
        pr[i].h = first_index;
        pr[i].e = exp_above[0];
        for (int n = 1; n < nb; n++)
        {
          pr[rel->nb].h = first_index + n;
          pr[rel->nb].e = exp_above[n];
          rel->nb++;
        }
      }
      else
        pr[i].h = renumber_get_index_from_p_r(renumber_tab, pr[i].p, r,
                                              pr[i].side);
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
 */
static void
get_outfilename_from_infilename (char *infilename, const char *outfmt,
                                 const char *outdir, char **oname,
                                 char **oname_tmp)
{
    const char * suffix_in;
    const char * suffix_out;
    get_suffix_from_filename (infilename, &suffix_in);
    suffix_out = outfmt ? outfmt : suffix_in;

    char * newname = strdup(infilename);
    ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
    newname[strlen(newname)-strlen(suffix_in)]='\0';

#define chkrcp(x) do { int rc = x; ASSERT_ALWAYS(rc>=0); } while (0)
    if(outdir) {
      const char * basename = path_basename(newname);
      chkrcp(asprintf(oname_tmp, "%s/%s.tmp%s", outdir, basename, suffix_out));
      chkrcp(asprintf(oname, "%s/%s%s", outdir, basename, suffix_out));
    } else {
      chkrcp(asprintf(oname_tmp, "%s.tmp%s", newname, suffix_out));
      chkrcp(asprintf(oname, "%s%s", newname, suffix_out));
    }
#undef  chkrcp

#if DEBUG >= 1
  fprintf (stderr, "DEBUG: Input file name: %s,\nDEBUG: temporary output file "
                   "name: %s,\nDEBUG: final output file name: %s\n", infilename,
                   *oname_tmp, *oname);
#endif
  free(newname);
}

static void
dup_print_stat (const char *s, uint64_t nrels, uint64_t ndup)
{
  uint64_t nrem = nrels - ndup;
  double pdup = 100.0 * ((double) ndup) / ((double) nrels);
  fprintf (stderr, "%s: nrels=%" PRIu64 " dup=%" PRIu64 " (%.2f%%) rem=%" PRIu64 "\n",
                   s, nrels, ndup, pdup, nrem);
}

static void *
hash_renumbered_rels (void * context_data MAYBE_UNUSED, earlyparsed_relation_ptr rel)
{
    unsigned int is_dup;

    nrels++;
    nrels_tot++;
    uint64_t i = insert_relation_in_dup_hashtable (rel, &is_dup);

    static unsigned long count = 0;

    // They should be no duplicate in already renumbered file
    if (is_dup && count++ < 10)
    {
        mpz_poly ab;
        mpz_poly_init(ab, EARLYPARSED_RELATION_MAX_AB - 1);
        mpz_poly_setcoeffs_int64(ab, rel->ab, EARLYPARSED_RELATION_MAX_AB - 1);
        fprintf (stderr, "Warning, duplicate relation in already renumbered files:"
                "\nab = ");
        mpz_poly_fprintf_coeffs(stderr, ab, ',');
        fprintf(stderr, "\ni = %" PRIu64 "\nj = %" PRIu32 "\n", i, H[i]);
        fprintf (stderr, "This warning may be due to a collision on the hash "
                "function or to an actual duplicate\nrelation. If it "
                "appears often you should check the input set of "
                "relations.\n\n");
        mpz_poly_clear(ab);
    }

    if (i < sanity_size)
        sanity_check (i, rel->ab[0], rel->ab[1]);

    if (cost >= factor * (double) (nrels_tot - ndup_tot))
        print_warning_size ();

    return NULL;
}

static void *
thread_dup2 (void * context_data, earlyparsed_relation_ptr rel)
{
    unsigned int is_dup;
    uint64_t i;
    FILE * output = (FILE*) context_data;
    nrels++;
    nrels_tot++;
    i = insert_relation_in_dup_hashtable (rel, &is_dup);
    if (!is_dup) {
        if (i < sanity_size)
            sanity_check(i, rel->ab[0], rel->ab[1]);
        if (cost >= factor * (double) (nrels_tot - ndup_tot))
            print_warning_size ();

        print_relation (output, rel);
    } else {
        ndup++;
        ndup_tot++;
    }

    return NULL;
}


void *
thread_root(void * context_data, earlyparsed_relation_ptr rel)
{
    if (!is_for_dl) { /* Do we reduce mod 2 */
        /* XXX should we compress as well ? */
        for (unsigned int i = 0; i < rel->nb; i++)
            rel->primes[i].e &= 1;
    }

    compute_index_rel (rel, context_data);

    return NULL;
}

int check_whether_file_is_renumbered(const char * filename, unsigned int npoly)
{
    unsigned int count = 0;
    char s[1024];
    FILE *f_tmp = fopen_maybe_compressed (filename, "rb");

    if (!f_tmp) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        abort();
    }

    if (feof (f_tmp)) /* file is empty */
      {
        fclose_maybe_compressed (f_tmp, filename);
        return 1; /* an empty file might be considered as renumbered */
      }

    /* Look for first non-comment line */
    while (1) {
      char *ret = fgets (s, 1024, f_tmp);
      if (ret == NULL)
      {
        /* fgets returns NULL when the end of file occurs or when there is an
           error */
        if (feof (f_tmp))
          {
            fclose_maybe_compressed (f_tmp, filename);
            return 1;
          }
        fprintf (stderr, "Error while reading %s\n", filename);
        exit (1);
      }
      if (strlen (s) >= 1023)
      {
        fprintf (stderr, "Too long line while reading %s\n", filename);
        exit (1);
      }
      size_t i = 0;
      while (s[i] == ' ')
        i++;
      if (s[i] != '#')
        break;
    }
    for (unsigned int i = 0; i < strlen (s); i++)
      count += s[i] == ':';
    fclose_maybe_compressed (f_tmp, filename);
    
    if (count == 1)
        return 1;
    else if (count == npoly)
        return 0;
    else
    {
      fprintf (stderr, "Error: invalid line in %s (line has %u colons but %u "
                       "were expected):\n %s", filename, count, npoly, s);
      abort();
    }
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "nrels",
                              "number of relations to be found in the slice");
  param_list_decl_usage(pl, "outdir", "by default, input files are overwritten");
  param_list_decl_usage(pl, "outfmt",
                               "format of output file (default same as input)");
  param_list_decl_usage(pl, "dl", "(switch) do not reduce exponents modulo 2");
  param_list_decl_usage(pl, "badidealinfo", "file containing info about bad ideals");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
    argv0 = argv[0];

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

    is_for_dl = 0; /* By default we do dup2 for factorization */
    param_list_configure_switch(pl, "dl", &is_for_dl);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char * outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * outdir = param_list_lookup_string(pl, "outdir");
    const char * renumberfilename = param_list_lookup_string(pl, "renumber");
    const char * badidealinfofile = param_list_lookup_string(pl, "badidealinfo");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    param_list_parse_ulong(pl, "nrels", &nrels_expected);

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (renumberfilename == NULL)
    {
      fprintf (stderr, "Error, missing -renumber command line argument\n");
      usage(pl, argv0);
    }
    if (badidealinfofile == NULL && is_for_dl) {
      fprintf (stderr, "Error, missing -badidealinfo command line argument\n");
      usage(pl, argv0);
    }
    if (nrels_expected == 0)
    {
      fprintf(stderr, "Error, missing -nrels command line argument "
                      "(or nrels = 0)\n");
      usage(pl, argv0);
    }
    K = 100 + nrels_expected + (nrels_expected / 5);
    K = uint64_nextprime (K);

    if (basepath && !filelist)
    {
      fprintf(stderr, "Error, -basepath only valid with -filelist\n");
      usage(pl, argv0);
    }
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "Error, output compression format unsupported\n");
        usage(pl, argv0);
    }
    if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
      usage(pl, argv0);
    }

    allbad_info_t badidealinfo;
    if (is_for_dl)
        read_bad_ideals_info(badidealinfofile, badidealinfo);

    set_antebuffer_path (argv0, path_antebuffer);

    renumber_init_for_reading (renumber_tab);
    renumber_read_table (renumber_tab, renumberfilename);

  /* sanity check: since we allocate two 64-bit words for each, instead of
     one 32-bit word for the hash table, taking K/100 will use 2.5% extra
     memory */
  sanity_size = 1 + (K / 100);
  fprintf (stderr, "[checking true duplicates on sample of %" PRIu64
           " cells]\n", sanity_size);
  sanity_a = (int64_t*)  malloc (sanity_size * sizeof (int64_t));
  if (sanity_a == NULL)
    {
      fprintf (stderr, "Error, cannot allocate sanity_a\n");
      exit (1);
    }
  memset (sanity_a, 0, sanity_size * sizeof (int64_t));

  sanity_b = (int64_t*) malloc (sanity_size * sizeof (uint64_t));
  if (sanity_b == NULL)
    {
      fprintf (stderr, "Error, cannot allocate sanity_b\n");
      exit (1);
    }

  H = (uint32_t*) malloc (K * sizeof (uint32_t));
  if (H == NULL)
    {
      fprintf (stderr, "Error, cannot allocate hash table\n");
      exit (1);
    }
  memset (H, 0, K * sizeof (uint32_t));
  fprintf (stderr, "Allocated hash table of %" PRIu64 " entries (%" PRIu64 "Mb)\n",
           K, (K * sizeof (uint32_t)) >> 20);

  /* Construct the two filelists : new files and already renumbered files */
  char ** files_already_renumbered, ** files_new;
  {
      unsigned int nb_files = 0;
      fprintf(stderr, "Constructing the two filelists...\n");
      char ** files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;
      for (char ** p = files; *p; p++)
          nb_files++;

      files_already_renumbered = malloc((nb_files + 1) * sizeof(char*));
      files_new = malloc((nb_files + 1) * sizeof(char*));

      /* separate already processed files
       * check if f_tmp is in raw format a,b:...:... or 
       *            in renumbered format a,b:... 
       */
      unsigned int nb_f_new = 0;
      unsigned int nb_f_renumbered = 0;
      for (char ** p = files; *p; p++) {
          /* always strdup these, so that we can safely call
           * filelist_clear in the end */
          if (check_whether_file_is_renumbered(*p, renumber_tab->nb_polys)) {
              files_already_renumbered[nb_f_renumbered++] = strdup(*p);
          } else {
              files_new[nb_f_new++] = strdup(*p);
          }
      }
      files_new[nb_f_new] = NULL;
      files_already_renumbered[nb_f_renumbered] = NULL;
      fprintf (stderr, "%u files (%u new and %u already renumbered)\n", nb_files, 
              nb_f_new, nb_f_renumbered);
      ASSERT_ALWAYS (nb_f_new + nb_f_renumbered == nb_files);
      /* if filelist was not given, then files == argv, which of course
       * must not be cleared */
      if (filelist) filelist_clear(files);
  }


  fprintf (stderr, "Reading files already renumbered:\n");
  filter_rels(files_already_renumbered,
          (filter_rels_callback_t) &hash_renumbered_rels,
          NULL,
          EARLYPARSE_NEED_AB, NULL, NULL);

  {
      struct filter_rels_description desc[3] = {
          { .f = thread_root, .arg=0, .n=4, },
          { .f = thread_dup2, .arg=0, .n=1, },
          { .f = NULL, },
      };
      if (is_for_dl)
          desc[0].arg = (void *) &badidealinfo[0];
      fprintf (stderr, "Reading new files"
              " (using %d auxiliary threads for roots mod p):\n",
              desc[0].n);

      for (char **p = files_new; *p ; p++) {
          FILE * output = NULL;
          char * oname, * oname_tmp;
          char * local_filelist[] = { *p, NULL};

          get_outfilename_from_infilename (*p, outfmt, outdir, &oname, &oname_tmp);
          output = fopen_maybe_compressed(oname_tmp, "w");
          desc[1].arg = (void*) output;

          nrels = ndup = 0;

          uint64_t loc_nrels = filter_rels2(local_filelist, desc,
                  EARLYPARSE_NEED_AB | EARLYPARSE_NEED_PRIMES,
                  NULL, NULL);

          ASSERT_ALWAYS(loc_nrels == nrels);

          fclose_maybe_compressed(output, oname_tmp);

#ifdef HAVE_MINGW /* For MinGW, rename cannot overwrite an existing file */
          remove (oname);
#endif
          if (rename(oname_tmp, oname))
          {
              fprintf(stderr, "Error while renaming %s into %s\n", oname_tmp, oname);
              abort();
          }

          // stat for the current file
          dup_print_stat (path_basename(*p), nrels, ndup);
          // stat for all the files already read
          dup_print_stat ("Total so far", nrels_tot, ndup_tot);

          free(oname);
          free(oname_tmp);
      }
  }

  fprintf (stderr, "At the end: %" PRIu64 " remaining relations\n",
                   nrels_tot - ndup_tot);

  fprintf (stderr, "At the end: hash table is %1.2f%% full\n"
                   "            hash table cost: %1.2f per relation\n",
                   100.0 * (double) (nrels_tot - ndup_tot) / (double) K,
                   1.0 + cost / (double) nrels_tot);
  fprintf (stderr, "  [found %lu true duplicates on sample of %lu relations]\n",
           sanity_collisions, sanity_checked);

  if (!*files_already_renumbered) {
      if (nrels_tot != nrels_expected) {
          fprintf(stderr, "Warning: number of relations read (%" PRIu64") does not match with the number of relations expected (%lu)\n", nrels_tot, nrels_expected);
      }
  } else {
      /* when we have renumbered files, we know that we won't have the
       * total number of relations... */
      if (nrels_tot > nrels_expected) {
          fprintf(stderr, "Warning: number of relations read (%" PRIu64") exceeds the number of relations expected (%lu)\n", nrels_tot, nrels_expected);
      }
  }

  if (is_for_dl)
      free(badidealinfo->badid_info);
  free (H);
  free (sanity_a);
  free (sanity_b);
  filelist_clear(files_already_renumbered);
  filelist_clear(files_new);

  param_list_clear(pl);
  renumber_clear (renumber_tab);
  return 0;
}
