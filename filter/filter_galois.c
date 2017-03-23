#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "filter_config.h"
#include "utils_with_io.h"
#include "mod_ul.h"

/* from dup2.c */
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

// Global variable for the table of Galois action
// For an ideal of index idx, Gal[idx] gives the index of a
// representative of the class under Galois action. It might be idx
// itself, or another index.
index_t *Gal;

/* Compute a table with exactly one representative per orbit under the
 * Galois action.
 */
#define c_swap(type, x, y) do { type _tmp = (x); (x) = (y); (y) = _tmp; } while (0)
static void
compute_galois_action (renumber_t tab, cado_poly_srcptr cpoly, galois_automorphism_srcptr sigma)
{
  Gal = (index_t *) malloc(tab->size * sizeof(index_t));
  ASSERT_ALWAYS(Gal != NULL);

  renumber_iterator it;
  renumber_iterator_init(it, tab, cpoly);
  for( ; !renumber_iterator_done(it) ; ) {
      index_t i0 = it->i;
      /* check special values. Those form their own orbit under Galois,
       * at least for now (well, theoretically we could have Galois
       * action on bad ideals after all...).
       */
      if (it->p == 0 || it->side < 0) {
          Gal[it->i] = it->i;
          renumber_iterator_next(it);
          continue;
      }
      /* We want to read a complete range of roots for a given (p,side) */
      p_r_values_t p = it->p;
      int side = it->side;
      p_r_values_t r[32];
      index_t ind[32];
      int nr = 0;
      for( ; p == it->p && side == it->side ; renumber_iterator_next(it)) {
          /* we begin with the simple ordering */
          ASSERT_ALWAYS(i0 + nr == it->i);
          Gal[it->i] = it->i;
          ind[nr] = it->i;
          r[nr] = it->r;
          nr++;
      }
      /* At this point we've advanced our iterator to the next (p,side) */
      if (nr % sigma->order) {
          fprintf(stderr,
                  "Warning: number of roots above"
                  " p=0x%" PRpr " not divisible by %d, skipping\n",
                  p, sigma->order);
          continue;
      }
      /* sort into orbits. There's rather similar code in
       * skip_galois_roots in las.cpp ; here, we're keeping orbits
       * contiguous, although I doubt it really matters.
       */
      /* So far the Gal[] table for this (p,side) is still trivial, and
       * ind[] points to the linearly-ordered range. */
      mpz_t pz, rz;
      mpz_init_set_ui(pz, p);
      mpz_init(rz);
      for(int i = 0 ; i < nr; ) {
          /* Let's say this root is the privileged one. We'll find all
           * roots in its orbit, and stow them close by */
          mpz_set_ui(rz, r[i]);
          int osize = 1;
          for( ; osize < sigma->order ; osize++) {
              galois_automorphism_apply_root(sigma, rz, rz, pz);
              /* where is it ? */
              int found;
              for(found = 0 ; found < nr ; found++) {
                  if (mpz_cmp_ui(rz, r[found]) == 0)
                      break;
              }
              ASSERT_ALWAYS(found < nr);
              /* If already stored, then this orbits wraps around
               * soonish, meaning that we're done.  */
              if (found <= i)
                  break;
              /* we want to put this root in position [i+osize]. This
               * entails changing our ind[] table */
              c_swap(p_r_values_t, r[found],   r[i+osize]);
              c_swap(index_t,      ind[found], ind[i+osize]);
              Gal[ind[i+osize]]=ind[i];
          }
          i += osize;
      }
      mpz_clear(rz);
      mpz_clear(pz);
  }
  renumber_iterator_clear(it);
}

static inline uint32_t
insert_relation_in_dup_hashtable (earlyparsed_relation_srcptr rel,
				  unsigned int *is_dup,
                                  galois_automorphism_srcptr sigma)
{
  uint32_t i, j;
  int64_t ab[EARLYPARSED_RELATION_MAX_AB];
  memcpy(ab, rel->ab, EARLYPARSED_RELATION_MAX_AB * sizeof(int64_t));

  uint64_t h = 0;
  for(int i = 0 ; i < sigma->order ; i++) {
      uint64_t dh=0;
      for(int i = 0 ; i < EARLYPARSED_RELATION_MAX_AB ; i++) {
          dh += constants_ab_dup2[i] * ab[i];
      }
      h ^= dh;
      galois_automorphism_apply_ab(sigma, (uint64_t*)&(ab[0]), &ab[1]);
  }
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

struct thread_galois_context_data {
    FILE * output;
    galois_automorphism_srcptr sigma;
};

/* I think that this is doomed. Just getting by with only a change of
 * sign works only for the narrowish situation where the automorphism has
 * order 2.
 */
static void *
thread_galois (void * context_data, earlyparsed_relation_ptr rel)
{
  unsigned int is_dup;
  struct thread_galois_context_data* context = context_data;
  FILE * output = context->output;
  galois_automorphism_srcptr sigma = context->sigma;

  /* see comment above. */
  ASSERT_ALWAYS(sigma->order == 2);

  insert_relation_in_dup_hashtable (rel, &is_dup, sigma);
  if (is_dup)
    return NULL;

  noutrels++;
  
  /* We're going to do ugly things.  */

  /* taken from dup2.c -- I wish there were something more generic. The
   * function in dup2 does extra stuff with J columns which is completely
   * irrelevant here.
   */
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;

  int nab = EARLYPARSED_RELATION_MAX_AB;
  for( ; nab > 2 && !rel->ab[nab-1] ; nab--);
  p = buf;
  *p++ = 'X';
  for(int i = 0 ; i < nab ; i++) {
      *p++ = i ? ',' : ' ';
      p = d64toa16(p, rel->ab[i]);
  }
  *p++ = ':';

  for (i = 0; i < rel->nb; i++)
  {
    ASSERT_ALWAYS(rel->primes[i].e != 0);
    index_t h = rel->primes[i].h;
    index_t hrep = Gal[h];
    // The new sign of the exponent is the XOR of the original sign and of 
    // the fact that we change the prime ideal for its conjugate.
    int neg = (rel->primes[i].e < 0) ^ (hrep != h);
    op = p;
    /* We prepend the prime with a "-" to indicate that the exponent is
     * to be understood differently...
     */
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
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "outdir", "by default, input files are overwritten");
  param_list_decl_usage(pl, "outfmt",
      "format of output file (default same as input)");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
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

int
main (int argc, char *argv[])
{
  argv0 = argv[0];
  cado_poly cpoly;
  unsigned long nrels_expected = 0;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "force-posix-threads",
      &filter_rels_force_posix_threads);

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

  if( action == NULL 
      || (strcmp(action, "1/y") 
	  && strcmp(action, "_y") 
	  && strcmp(action, "autom3.1g")
	  && strcmp(action, "autom3.2g")
	 )
    )
  {
    fprintf(stderr, "Error, missing -action command line argument\n");
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

  renumber_init_for_reading (renumber_tab);
  renumber_read_table (renumber_tab, renumberfilename);

  fprintf(stderr, "Computing Galois action %s on ideals\n", action);

  galois_automorphism_srcptr sigma = galois_automorphism_get(action);

  compute_galois_action(renumber_tab, cpoly, sigma);

  fprintf(stderr, "Rewriting relations files\n");
  char ** files;
  unsigned int nb_files = 0;
  files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;
  for (char ** p = files; *p; p++)
    nb_files++;

  struct filter_rels_description desc[2] = {
    { .f = thread_galois, .arg = NULL, .n=1, },
    { .f = NULL, },
  };


  fprintf (stderr, "Reading files (using %d auxiliary threads):\n", desc[0].n);
  for (char **p = files; *p ; p++) {
    FILE * output = NULL;
    char * oname, * oname_tmp;
    char * local_filelist[] = { *p, NULL};

    get_outfilename_from_infilename (*p, outfmt, outdir, &oname, &oname_tmp);
    output = fopen_maybe_compressed(oname_tmp, "w");
    
    struct thread_galois_context_data ctx = {
        .output = output,
        .sigma = sigma,
    };

    desc[0].arg = &ctx;

    filter_rels2(local_filelist, desc,
        EARLYPARSE_NEED_AB | EARLYPARSE_NEED_INDEX,
        NULL, NULL);

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
  renumber_clear (renumber_tab);
  cado_poly_clear (cpoly);
  free(H);
  return 0;
}
