#include "cado.h" // IWYU pragma: keep

#include <cstdint>          // for int64_t, uint64_t, uint32_t
#include <cstdio>            // for fprintf, stderr, asprintf, FILE
#include <cstdlib>           // for free, abort, exit, malloc, EXIT_FAILURE
#include <cstring>          // for strcmp, strlen, memcpy, memset, strdup

#include <algorithm>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>           // for string
#include <vector>           // for vector

#ifdef HAVE_MINGW
#include <fcntl.h>   /* for _O_BINARY */
#endif

#include "cado_poly.h"       // for cado_poly_clear, cado_poly_init, cado_po...
#include "filter_config.h"   // for CA_DUP2, CB_DUP2
#include "filter_io.h"       // for earlyparsed_relation_s, filter_rels_desc...
#include "fmt/base.h"        // for fmt::print
#include "galois_action.hpp"  // for galois_action
#include "gzip.h"            // for fclose_maybe_compressed, fopen_maybe_com...
#include "macros.h"          // for ASSERT_ALWAYS, UNLIKELY
#include "misc.h"            // for filelist_clear, filelist_from_file
#include "params.h"          // for param_list_decl_usage, param_list_lookup...
#include "portability.h" // strdup  // IWYU pragma: keep
#include "relation-tools.h"  // for u64toa16, d64toa16
#include "renumber.hpp"      // for renumber_t, renumber_t::p_r_side
#include "timing.h"          // timingstats_dict_t
#include "typedefs.h"        // for p_r_values_t, index_t, prime_t, PRpr
#include "verbose.h"         // for verbose_decl_usage, verbose_interpret_pa...


// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
static char const * argv0; /* = argv[0] */

static std::unique_ptr<uint32_t[]> H; /* H contains the hash table */
static unsigned long K = 0; /* Size of the hash table */
static unsigned long noutrels = 0; // Number of output relations
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

/* return in oname and oname_tmp two file names for writing the output
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
static std::string
get_outfilename_from_infilename (std::string const & infilename,
        std::string const & outfmt,
    std::string const & outdir)
{
  const char * suffix_in;
  get_suffix_from_filename (infilename.c_str(), &suffix_in);
  std::string suffix_out = outfmt;
  if (suffix_out.empty()) suffix_out = suffix_in;

  std::string newname = infilename.substr(0, infilename.size() - strlen(suffix_in));

  std::string prefix = outdir;

  if (!prefix.empty()) {
      if (prefix.back() != '/')
          prefix += '/';
      newname = path_basename(newname.c_str());
  }

  return prefix + newname + suffix_out;
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

struct thread_galois_arg {
  std::vector<index_t> const & ga_id_cache;
  galois_action const & gal_action;
  std::ostream & os;
};

static void *
thread_galois (void * context_data, earlyparsed_relation_ptr rel)
{
  unsigned int is_dup;
  auto const & data = *(thread_galois_arg const *) context_data;

  const std::vector<index_t> &sigma = data.ga_id_cache;
  const galois_action &G = data.gal_action;
  std::ostream & output(data.os);
  insert_relation_in_dup_hashtable (rel, &is_dup, G);
  if (is_dup)
    return nullptr;

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
    index_t const h = rel->primes[i].h;
    index_t const hrep = sigma[h];
    // The new sign of the exponent is the XOR of the original sign and of 
    // the fact that we change the prime ideal for its conjugate.
    int const neg = (rel->primes[i].e < 0) ^ (hrep != h);
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
  if (!(output << buf)) {
    perror("Error writing relation");
    abort();
  }
  return nullptr;
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
usage (param_list pl, char const * argv0)
{
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

// coverity[root_function]
int
main (int argc, char const * argv[])
{
  argv0 = argv[0];
  cxx_cado_poly cpoly;
  unsigned long nrels_expected = 0;
  int for_dl = 0;

  cxx_param_list pl;
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
  std::string outfmt;
  param_list_parse(pl, "outfmt", outfmt);
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  std::string outdir;
  param_list_parse(pl, "outdir", outdir);
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
  const char * action = param_list_lookup_string(pl, "galois");
  param_list_parse_ulong(pl, "nrels", &nrels_expected);

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }
  if (polyfilename == nullptr)
  {
    fprintf(stderr, "Error, missing -poly command line argument\n");
    usage(pl, argv0);
  }
  if (renumberfilename == nullptr)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage(pl, argv0);
  }
  if (basepath && !filelist)
  {
    fprintf(stderr, "Error, -basepath only valid with -filelist\n");
    usage(pl, argv0);
  }
  if (!outfmt.empty() && !is_supported_compression_format(outfmt.c_str())) {
    fprintf(stderr, "Error, output compression format unsupported\n");
    usage(pl, argv0);
  }
  if ((filelist != nullptr) + (argc != 0) != 1) {
    fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
    usage(pl, argv0);
  }
  if (nrels_expected == 0)
  {
    fprintf(stderr, "Error, missing -nrels command line argument "
                    "(or nrels = 0)\n");
    usage(pl, argv0);
  }
  K = 100 + 1.2 * double(nrels_expected);

  if(action == nullptr)
  {
    fprintf(stderr, "Error, missing -galois command line argument\n");
    usage(pl, argv0);
  }

  H = std::unique_ptr<uint32_t[]>(new uint32_t[K]);
  ASSERT_ALWAYS(H);
  std::fill(H.get(), H.get() + K, 0);

  if (!cado_poly_read (cpoly, polyfilename)) {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  set_antebuffer_path (argv0, path_antebuffer);

  /* */
  galois_action gal_action(action);
  fmt::print("# Using {}\n", gal_action);

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
                   "your use case\nSee comments in utils/galois_action.cpp "
                   "for more info\n/!\\/!\\/!\\/!\\\n\n";
  }

  std::cout << "Computing Galois action on ideals\n";
  std::vector<index_t> ga_id_cache;
  size_t norb = gal_action.compute_action_on_index(ga_id_cache, renumber_tab);
  fmt::print("Found {} orbits of length {}, {} columns were left unchanged "
             "(among which {} additional column(s) and {} column(s) "
             "corresponding to badideals)\n", norb, gal_action.get_order(),
             ga_id_cache.size() - norb*gal_action.get_order(),
             renumber_tab.number_of_additional_columns(),
             renumber_tab.number_of_bad_ideals());

  std::cout << "Rewriting relations files\n";

  char const ** files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;

  std::cout << "Reading files (using 1 auxiliary thread):\n";
  timingstats_dict_t stats;
  timingstats_dict_init(stats);

  for (char const **p = files; *p ; p++) {
    // FILE * output = nullptr;
    std::string oname;
    const char * local_filelist[] = { *p, nullptr};

    oname = get_outfilename_from_infilename (*p, outfmt, outdir);

    ofstream_maybe_compressed output(oname);
    // output = fopen_maybe_compressed(oname_tmp, "w");
    if (!output) {
        fmt::print (stderr, "Error, could not open file to write the relations. "
                       "Check that the directory {} exists\n", outdir);
      abort();
    }

    thread_galois_arg foo { ga_id_cache, gal_action, output };

    filter_rels(local_filelist, (filter_rels_callback_t) &thread_galois,
                (void *) &foo, EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX,
                nullptr, stats);

    // fclose_maybe_compressed(output, oname_tmp);

  }

  fprintf(stderr, "Number of output relations: %lu\n", noutrels);

  if (filelist)
    filelist_clear(files);


  timingstats_dict_add_mythread(stats, "main");
  timingstats_dict_disp(stats);
  timingstats_dict_clear(stats);
  return 0;
}
