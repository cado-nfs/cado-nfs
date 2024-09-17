#include "cado.h" // IWYU pragma: keep
#include <cinttypes>        // for PRId64, PRIu64
#include <cstdint>          // for int64_t, uint64_t, uint32_t
#include <cstring>          // for strcmp, strlen, memcpy, memset, strdup
#include <cstdio>            // for fprintf, stderr, NULL, asprintf, FILE
#include <cstdlib>           // for free, abort, exit, malloc, EXIT_FAILURE
#include <unordered_map>    // for unordered_map
#include <tuple>            // for tuple, get
#include <vector>           // for vector
#include <string>           // for string
#ifdef HAVE_MINGW
#include <fcntl.h>   /* for _O_BINARY */
#endif
#include "cado_poly.h"       // for cado_poly_clear, cado_poly_init, cado_po...
#include "filter_config.h"   // for CA_DUP2, CB_DUP2
#include "filter_io.h"       // for earlyparsed_relation_s, filter_rels_desc...
#include "galois_utils.hpp"  // for galois_action
#include "gzip.h"            // for fclose_maybe_compressed, fopen_maybe_com...
#include "macros.h"          // for ASSERT_ALWAYS, UNLIKELY
#include "misc.h"            // for filelist_clear, filelist_from_file
#include "mod_ul.h"          // for modul_clear, modul_clearmod, modul_get_ul
#include "params.h"          // for param_list_decl_usage, param_list_lookup...
#include "portability.h" // strdup  // IWYU pragma: keep
#include "relation-tools.h"  // for u64toa16, d64toa16
#include "renumber.hpp"      // for renumber_t, renumber_t::p_r_side
#include "timing.h"          // timingstats_dict_t
#include "typedefs.h"        // for p_r_values_t, index_t, prime_t, PRpr
#include "verbose.h"         // for verbose_decl_usage, verbose_interpret_pa...


char *argv0; /* = argv[0] */

static uint32_t *H; /* H contains the hash table */
static unsigned long K = 0; /* Size of the hash table */
static unsigned long noutrels = 0; // Number of output relations

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

static inline uint32_t
insert_relation_in_dup_hashtable (earlyparsed_relation_srcptr rel,
				  unsigned int *is_dup, const galois_action &G)
{
  uint64_t h;
  uint32_t i, j;

  h = G.hash_ab(rel->a, rel->b, CA_DUP2, CB_DUP2);
  i = h % K;
  j = (uint32_t) (h >> 32);
  while (H[i] != 0 && H[i] != j) {
    i++;
    if (UNLIKELY(i == K))
      i = 0;
  }

  if (H[i] == j) {
    *is_dup = 1;
  } else {
    H[i] = j;
    *is_dup = 0;
  }
  return i;
}

static void *
thread_galois (void * context_data, earlyparsed_relation_ptr rel)
{
  unsigned int is_dup;
  auto data = *((std::tuple<std::vector<index_t> const &, galois_action const &, FILE *> *) context_data);
  const std::vector<index_t> &sigma = std::get<0>(data);
  const galois_action &G = std::get<1>(data);
  FILE *output = std::get<2>(data);
  insert_relation_in_dup_hashtable (rel, &is_dup, G);
  if (is_dup)
    return NULL;

  noutrels++;
  
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;

  p = d64toa16(buf, rel->a);
  *p++ = ',';
  p = u64toa16(p, rel->b);
  *p++ = ':';

  for (i = 0; i < rel->nb; i++)
  {
    ASSERT_ALWAYS(rel->primes[i].e != 0);
    index_t h = rel->primes[i].h;
    index_t hrep = sigma[h];
    // The new sign of the exponent is the XOR of the original sign and of 
    // the fact that we change the prime ideal for its conjugate.
    int neg = (rel->primes[i].e < 0) ^ (hrep != h);
    op = p;
    if (neg) { *p++ = '-'; }
    p = u64toa16(p, (uint64_t) hrep);
    *p++ = ',';
    t = p - op;
    if (rel->primes[i].e < 0) {
      j = (unsigned int) ((-rel->primes[i].e) - 1);
    } else {
      j = (unsigned int) (rel->primes[i].e - 1);
    }
    while (j != 0) {
      memcpy(p, op, t);
      p += t;
      j--;
    }
  }

  *(--p) = '\n';
  p[1] = 0;
  if (fputs(buf, output) == EOF) {
    perror("Error writing relation");
    abort();
  }
  return NULL;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "dl", "for DL (untested)");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "outdir", "by default, input files are overwritten");
  param_list_decl_usage(pl, "outfmt",
      "format of output file (default same as input)");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "nrels", "(approximate) number of input relations");
  param_list_decl_usage(pl, "galois", "Galois action among 1/y or _y");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

// coverity[root_function]
int
main (int argc, char *argv[])
{
  argv0 = argv[0];
  cado_poly cpoly;
  unsigned long nrels_expected = 0;
  int for_dl = 0;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "force-posix-threads",
      &filter_rels_force_posix_threads);
  param_list_configure_switch(pl, "dl", &for_dl);

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
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * outfmt = param_list_lookup_string(pl, "outfmt");
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * outdir = param_list_lookup_string(pl, "outdir");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
  const char * action = param_list_lookup_string(pl, "galois");
  param_list_parse_ulong(pl, "nrels", &nrels_expected);

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }
  if (polyfilename == NULL)
  {
    fprintf(stderr, "Error, missing -poly command line argument\n");
    usage(pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage(pl, argv0);
  }
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
  if (nrels_expected == 0)
  {
    fprintf(stderr, "Error, missing -nrels command line argument "
                    "(or nrels = 0)\n");
    usage(pl, argv0);
  }
  K = 100 + 1.2 * nrels_expected;

  if(action == NULL)
  {
    fprintf(stderr, "Error, missing -galois command line argument\n");
    usage(pl, argv0);
  }

  H = (uint32_t*) malloc (K * sizeof (uint32_t));
  ASSERT_ALWAYS(H);
  memset (H, 0, K * sizeof (uint32_t));

  cado_poly_init (cpoly);
  if (!cado_poly_read (cpoly, polyfilename)) {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  set_antebuffer_path (argv0, path_antebuffer);

  /* */
  galois_action gal_action(action);
  std::cout << "# Using " << gal_action << std::endl;

  if (gal_action.get_order() > 2) {
      std::cerr << "Error, Galois action of order > 2 are not supported "
                   "yet. The missing piece of code is the one that takes "
                   "care of computing the new valuation of the the "
                   "rewritten ideals.\n";
      exit(EXIT_FAILURE);
  }

  /* Renumbering table to convert from (p,r) to an index */
  renumber_t renumber_tab(cpoly);
  renumber_tab.read_from_file(renumberfilename, for_dl);

  if (renumber_tab.number_of_bad_ideals() > 0 && gal_action.get_order() > 1) {
      std::cout << "\n/!\\/!\\/!\\/!\\\nWARNING, bad ideals will be left "
                   "unchanged, the output may not be usable depending on "
                   "your use case\nSee comments in utils/galois_utils.cpp "
                   "for more info\n/!\\/!\\/!\\/!\\\n\n";
  }

  std::cout << "Computing Galois action on ideals\n";
  std::vector<index_t> ga_id_cache;
  size_t norb = gal_action.compute_action_on_index(ga_id_cache, renumber_tab);
  std::cout << "Found " << norb << " orbits of length "
            << gal_action.get_order() << ", "
            << ga_id_cache.size() - norb*gal_action.get_order() << " columns "
            << "were left unchanged (among which "
            <<  renumber_tab.number_of_additional_columns() << " additional "
            << "column(s) and " << renumber_tab.number_of_bad_ideals()
            << " column(s) corresponding to badideals)\n";

  std::cout << "Rewriting relations files\n";
  char ** files;
  files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;

  std::cout << "Reading files (using 1 auxiliary thread):\n";
  timingstats_dict_t stats;
  timingstats_dict_init(stats);

  for (char **p = files; *p ; p++) {
    FILE * output = NULL;
    char * oname, * oname_tmp;
    char * local_filelist[] = { *p, NULL};

    get_outfilename_from_infilename (*p, outfmt, outdir, &oname, &oname_tmp);
    output = fopen_maybe_compressed(oname_tmp, "w");
    if (output == NULL){
      fprintf (stderr, "Error, could not open file to write the relations. "
                       "Check that the directory %s exists\n", outdir);
      abort();
    }
    auto arg = std::make_tuple(std::cref(ga_id_cache), std::cref(gal_action),
                                                       output);

    filter_rels(local_filelist, (filter_rels_callback_t) &thread_galois,
                (void *) &arg, EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX,
                NULL, stats);

    fclose_maybe_compressed(output, oname_tmp);

#ifdef HAVE_MINGW /* For MinGW, rename cannot overwrite an existing file */
    remove (oname);
#endif
    if (rename(oname_tmp, oname)) {
      fprintf(stderr, "Error while renaming %s into %s\n", oname_tmp, oname);
      abort();
    }

    free(oname);
    free(oname_tmp);
  }

  fprintf(stderr, "Number of output relations: %lu\n", noutrels);

  if (filelist)
    filelist_clear(files);

  param_list_clear(pl);
  cado_poly_clear (cpoly);
  free(H);

  timingstats_dict_add_mythread(stats, "main");
  timingstats_dict_disp(stats);
  timingstats_dict_clear(stats);
  return 0;
}
