/* replay_dblemat -- Kleinjung's double-matrix strategy.

Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2015
          Francois Morain, Emmanuel Thome, Paul Zimmermann,
          Cyril Bouvier, Pierrick Gaudry

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/


/* Notations:
   
   We have in input the matrix MM corresponding to the .purged file.
   Merge has computed a sequence of row operations on MM, stored in
   the .his file, such that the resulting matrix is the sparse matrix
   M_sp given to Wiedemann. In terms of matrix operations, we have
        M_sp = Y . (Prod mu_i)  MM . X,
   where
        - each mu_i is an elementary matrix encoding one line of the .his file,
        - X is a matrix that extracts columns,
        - Y is a matrix that extracts rows.
   Kleinjung's double-matrix strategy consist in writing M_sp as a product 
   of two low-weight matrices:
        M_sp = M1 . M2,
   where
        - M2 is basically MM to which we have applied the first merges.
        - M1 is the product of all the elementary matrices corresponding
        to other merges.

   Sizes:
        - MM is an nr x nc matrix.
        - M_sp ....
        - M1 ...
        - M2 ...
*/


#include "cado.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>   /* for _O_BINARY */
#include "filter_config.h"
#include "filter_io.h"  // earlyparsed_relation_ptr
#include "fix-endianness.h" // fwrite32_little
#include "gzip.h"       // fopen_maybe_compressed
#include "merge_replay_matrix.h"
#include "misc.h"       // derived_filename
#include "params.h" // param_list
#include "sparse.h"
#include "stats.h"                     // stats_data_t
#include "portability.h"

/* An alias for the main matrix type */
typedef typerow_t ** matrix_t;



/******************************************/
/* Functions for reading the .purged file */
/******************************************/

typedef struct {
  matrix_t mat;
  uint64_t ncols;
} replay_read_data_t;

void * fill_in_rows (void *context_data, earlyparsed_relation_ptr rel) {
  replay_read_data_t *data = (replay_read_data_t *) context_data;

  data->mat[rel->num] = (typerow_t*) malloc((rel->nb + 1)*sizeof (typerow_t));
  ASSERT_ALWAYS(data->mat[rel->num] != NULL);

  rowLength(data->mat, rel->num) = rel->nb;
  for(unsigned int j = 0; j < rel->nb; j++) {
    index_t h = rel->primes[j].h;
    exponent_t e = rel->primes[j].e;
    ASSERT_ALWAYS (e == 1);
    ASSERT_ALWAYS (h < data->ncols);
    setCell(data->mat[rel->num], j+1, h, e);
  }
  qsort(&(data->mat[rel->num][1]), rel->nb, sizeof(typerow_t), cmp_typerow_t);
  return NULL;
}

/* Reads purged file; fill-in the matrix mat.
 * This assumes that mat is pre-allocated, and that everything
 * is consistent. */
static void
read_purgedfile(matrix_t mat, const char* filename, uint64_t nrows,
                 uint64_t ncols)
{
  uint64_t nread;
  printf("Reading sparse matrix from %s\n", filename);
  fflush(stdout);
  char *fic[2] = {(char *) filename, NULL};
  replay_read_data_t tmp = (replay_read_data_t) {.mat= mat, .ncols = ncols};
  nread = filter_rels(fic, (filter_rels_callback_t) &fill_in_rows, &tmp,
      EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS (nread == nrows);
}

/*********************/
/* Replaying history */
/*********************/

// A line is "i i1 ... ik [#j]".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
// If given, j is the index of the column used for pivoting (used in DL).

static void
doAllAdds(matrix_t M, const char *str, index_data_t index_data)
{
  index_t j;
  int32_t ind[MERGE_LEVEL_MAX], i0;
  int ni, destroy;

  ni = parse_hisfile_line(ind, str, &j);

  if (ind[0] < 0) {
    destroy = 0;
    i0 = -ind[0]-1;
  } else {
    destroy = 1;
    i0 = ind[0];
  }

  for (int k = 1; k < ni; k++)
    addRowsUpdateIndex(M, index_data, ind[k], i0, j);

  if (destroy) {
    //destroy initial row
    free(M[i0]);
    M[i0] = NULL;

    if (index_data != NULL) // ie we want an index
    {
      index_data[i0].n = 0;
      free(index_data[i0].rels);
      index_data[i0].rels = NULL;
    }
  }
}

// M1 := M2
static void
copy_mat(matrix_t M1, const matrix_t M2, uint64_t nr)
{
  for (uint64_t i = 0; i < nr; ++i) {
    if (M2[i] == NULL) {
      M1[i] = NULL;
    } else {
      M1[i] = (typerow_t *) realloc(M1[i],
          ((rowLength(M2, i))+1)*sizeof(typerow_t));
      ASSERT_ALWAYS(M1[i] != NULL);
      memcpy(M1[i], M2[i], ((rowLength(M2, i))+1)*sizeof(typerow_t));
    }
  }
}


static index_data_t init_index_data(uint64_t nr) {
  // At the beginning, the index_data consists of relsets that
  // are just single relations.
  index_data_t index_data = (relset_t *) malloc (nr*sizeof(relset_t));
  ASSERT_ALWAYS(index_data);
  for (uint64_t i = 0; i < nr; ++i) {
    index_data[i].n = 1;
    index_data[i].rels = (multirel_t *)malloc(sizeof(multirel_t));
    index_data[i].rels[0].ind_row = i;
    index_data[i].rels[0].e = 1;
  }
  return index_data;
}

static void 
apply_hisfile(matrix_t MM, matrix_t M1, matrix_t M2, uint64_t nr, 
    const char *hisname, index_data_t index_data, uint64_t split_point)
{
  FILE *hisfile;
  hisfile = fopen_maybe_compressed (hisname, "r");
  ASSERT_ALWAYS(hisfile != NULL);

#define STRLENMAX 2048
  uint64_t addread = 0;
  char str[STRLENMAX];
  int startedM1 = 0;

  printf("Reading row additions\n");
  stats_data_t stats; /* struct for printing progress */
  /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
   * 2^23 primes after that */
  stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "line");
  while(fgets(str, STRLENMAX, hisfile)) {
    if (str[0] == '#') continue;
    addread++;
    if (stats_test_progress(stats))
      stats_print_progress (stats, addread, 0, 0, 0);
    /* If incomplete line, there is a bug somewhere. Let's crash. */
    ASSERT_ALWAYS(str[strlen(str)-1] == '\n'); 

    if (!startedM1) {
      if (split_point == 0) {
        // Try to detect first 3-merge
        /* The first 3-merge is detected by the fact that we don't erase a
         * row. TODO: This is not robust! 
         * Anyway, using this default is really suboptimal...*/
        if (str[0] == '-') {
          fprintf(stderr, "Reached the first 3-merge\n");
          startedM1 = 1;
          // Save the current matrix in M2
          copy_mat(M2, MM, nr);
        }
      } else {
        if (addread == split_point) {
          fprintf(stderr, "Reached stop point\n");
          startedM1 = 1;
          // Save the current matrix in M2
          copy_mat(M2, MM, nr);
        }
      }
    }

    /* Do the operation on the main matrix. */
    doAllAdds(MM, str, index_data);

    /* If we are after the stop point, record operations in M1 */
    if (startedM1)
      doAllAdds(M1, str, NULL);
  }

  /* If split_point was beyond the end of history file, the user
   * does not seem to want dblemat. We still have to copy MM to M2.
   * In that case, M1 is just the identity matrix. */
  if (!startedM1) {
    copy_mat(M2, MM, nr);
  }

  stats_print_progress (stats, addread, 0, 0, 1);
#undef STRLENMAX
}

/***********************************/
/* Renumbering / counting weights  */
/***********************************/

static uint64_t
total_weight(const matrix_t M, const uint64_t nr, const uint32_t skip)
{
  uint64_t wt = 0;
  for (uint64_t i = 0; i < nr; i++) {
    if (M[i] != NULL) {
      for (uint32_t k = 1; k <= rowLength(M, i); k++) {
        if (rowCell(M[i], k) >= skip) {
          wt++;
        }
      }
    }
  }
  return wt;
}

// Compute the column weights of M, store them in col_weight, which must
// have been pre-allocated with size nc (the number of columns of M).
static void
compute_col_weights(const matrix_t M, uint32_t *col_weight,
    const uint64_t nr, const uint64_t nc)
{
  memset(col_weight, 0, nc * sizeof(uint32_t));
  for (uint64_t i = 0; i < nr; i++) {
    if (M[i] != NULL)
      for (uint32_t k = 1; k <= rowLength(M, i); k++) {
        ASSERT(rowCell(M[i], k) < nc);
        col_weight[rowCell(M[i], k)] += 1;
      }
  }
}

// Count non-empty rows in matrix M
static uint64_t
count_non_empty_rows(const matrix_t M, const uint64_t nr)
{
  uint64_t c = 0;
  for (uint64_t i = 0; i < nr; ++i) {
    if (M[i] != NULL)
      c++;
  }
  return c;
}

// Count non-empty cols in matrix M
static uint64_t
count_non_empty_cols(MAYBE_UNUSED const matrix_t M,
    const uint32_t *col_weight, const uint64_t nc)
{
  uint64_t c = 0;
  for (uint64_t i = 0; i < nc; ++i) {
    if (col_weight[i] != 0)
      c++;
  }
  return c;
}


// Compare pairs of uint64_t, lexicographical order.
// This is inverted: return 1 if *p is less than *q.
static inline int
cmp_uint64_2 (const void *p, const void *q)
{
  uint64_t *x = (uint64_t*) p;
  uint64_t *y = (uint64_t*) q;

  if (x[0] > y[0])
    return -1;
  else if (x[0] < y[0])
    return 1;
  else
    return (x[1] > y[1]) ? 1 : -1;
}

// Compute the permutation corresponding to sorting 
// the columns according to their weight.
// Rem: col_weight contains weights (hence 32 bits), while the resulting 
// array contains indices (hence 64 bits).
static uint64_t *
renumber_decr_weights(const uint32_t *col_weight, const uint64_t nc)
{
  uint64_t *tmp = (uint64_t *)malloc(2*nc*sizeof(uint64_t));
  ASSERT_ALWAYS(tmp != NULL);
  memset(tmp, 0, 2*nc*sizeof(uint64_t));

  for (uint64_t i = 0; i < nc; ++i) {
    tmp[2*i] = col_weight[i];
    tmp[2*i+1] = i;
  }
  qsort(tmp, nc, 2*sizeof(uint64_t), cmp_uint64_2);
  for (uint64_t i = 0; i < nc; ++i)
    tmp[i] = tmp[2*i+1];
  tmp = (uint64_t *)realloc(tmp, nc*sizeof(uint64_t));
  ASSERT_ALWAYS(tmp != NULL);
  uint64_t *col_renum = (uint64_t *)malloc(nc*sizeof(uint64_t));
  ASSERT_ALWAYS(col_renum != NULL);
  for (uint64_t i = 0; i < nc; ++i)
    col_renum[i] = UINT64_MAX;

  for (uint64_t i = 0; i < nc; ++i)
    col_renum[tmp[i]] = i;

  // Now, col_renum[x] contains the new index for the column number x.
  free(tmp);
  return col_renum;
}

/*****************************/
/* Writing matrices to files */
/*****************************/

static void
write_matrix(const char *basename, const matrix_t M, int nr, int nc,
    unsigned int skip, int bin, int zip) {

  const struct {
    const char * ext;
    const char * smat;
    const char * srw;
    const char * scw;
    const char * dmat;
    const char * drw;
    const char * dcw;
  } suffixes[2] = {
    {
      .ext = ".txt",
      .smat = "sparse.txt",
      .srw = "sparse.rw.txt",
      .scw = "sparse.cw.txt",
      .dmat = "dense.txt",
      .drw = "dense.rw.txt",
      .dcw = "dense.cw.txt",
    },
    {
      .ext = ".bin",
      .smat = "sparse.bin",
      .srw = "sparse.rw.bin",
      .scw = "sparse.cw.bin",
      .dmat = "dense.bin",
      .drw = "dense.rw.bin",
      .dcw = "dense.cw.bin",
    },
  }, * suf = &(suffixes[bin]);

  char wmode[3] = "w";
#ifdef HAVE_MINGW
  if (bin)
    strcpy (wmode, "wb");
#endif

  const char * gz = (zip) ? ".gz" : NULL;

  char * smatname = NULL;
  FILE * smatfile = NULL;
  char * srwname  = NULL;
  FILE * srwfile  = NULL;
  char * scwname  = NULL;
  FILE * scwfile  = NULL;

  smatname = derived_filename(basename, suf->smat, gz);
  srwname  = derived_filename(basename, suf->srw, gz);
  scwname = derived_filename(basename, suf->scw, gz);
  smatfile = fopen_maybe_compressed(smatname, wmode);
  srwfile  = fopen_maybe_compressed(srwname, wmode);
  if (!bin) fprintf(smatfile, "%d %d\n", nr, nc - skip);
  ASSERT_ALWAYS(smatfile);
  ASSERT_ALWAYS(srwfile);

  char * dmatname = NULL;
  FILE * dmatfile = NULL;
  char * drwname  = NULL;
  FILE * drwfile  = NULL;
  char * dcwname  = NULL;
  FILE * dcwfile  = NULL;

  if (skip) {
    dmatname = derived_filename(basename, suf->dmat, gz);
    drwname  = derived_filename(basename, suf->drw, gz);
    dcwname = derived_filename(basename, suf->dcw, gz);
    dmatfile = fopen_maybe_compressed(dmatname, wmode);
    drwfile  = fopen_maybe_compressed(drwname, wmode);
    if (!bin) fprintf(dmatfile, "%d %d\n", nr, skip);
    ASSERT_ALWAYS(dmatfile);
    ASSERT_ALWAYS(drwfile);
  }

  /* Variables to store information while we are writing to file.
   * In particular, we reconstruct the weights on the fly */
  unsigned long W = 0;
  unsigned long DW = 0;
  uint32_t * weights = malloc(nc * sizeof(uint32_t));
  ASSERT_ALWAYS(weights != NULL);
  memset(weights, 0, nc * sizeof(uint32_t));

  /* Writing the matrix, with possibly the first columns in the
   * so-called "dense" part, if the parameter skip is given */
  for(int i = 0; i < nr; i++){
    if(M[i] == NULL) {
      /* this should not happen, unless the corresponding combination of
         relations yields a square, thus we have found a dependency in the
         merge process */
      if (bin) {
        const uint32_t x = 0;
        fwrite32_little(&x, 1, smatfile);
        if (skip) fwrite32_little(&x, 1, dmatfile);
      } else {
        fprintf(smatfile, "0");
        if (skip) fprintf(dmatfile, "0");
      }
    } else {
      uint32_t dw = 0;
      uint32_t sw = 0;
      for(unsigned int j = 1; j <= rowLength(M, i); j++){
        if (rowCell(M, i, j) < skip) {
          dw++;
          DW++;
        } else {
          sw++;
          W++;
        }
      }
      if (bin) {
        fwrite32_little(&sw, 1, smatfile);
        if (srwfile) fwrite32_little(&sw, 1, srwfile);
        if (skip) fwrite32_little(&dw, 1, dmatfile);
        if (skip) fwrite32_little(&dw, 1, drwfile);
      } else {
        fprintf(smatfile, "%" PRIu32 "", sw);
        if (srwfile) fprintf(srwfile, "%" PRIu32 "\n", sw);
        if (skip) fprintf(dmatfile, "%" PRIu32 "", dw);
        if (skip) fprintf(drwfile, "%" PRIu32 "\n", dw);
      }
      for(unsigned int j = 1; j <= rowLength(M, i); j++){
        uint32_t x = rowCell(M, i, j);
        if (srwfile) weights[x]++;
        if (x < skip) {
          ASSERT_ALWAYS(skip);
          if (bin) {
            fwrite32_little(&x, 1, dmatfile);
          } else {
            fprintf(dmatfile, " %" PRIu32 "", x);
          }
        } else {
          x-=skip;
          if (bin) {
            fwrite32_little(&x, 1, smatfile);
#ifdef FOR_DL
            uint32_t e = (uint32_t) M[i][j].e;
            fwrite32_little(&e, 1, smatfile);
#endif
          } else {
            fprintf(smatfile, " %" PRIu32 "", x);
#ifdef FOR_DL
            fprintf(smatfile, ":%d", M[i][j].e);
#endif
          }
        }
      }
    }
    if (!bin) {
      fprintf(smatfile, "\n");
      if (skip) fprintf(dmatfile, "\n");
    }
  }
  
  /* finish, and write the column weights */
  fclose_maybe_compressed(smatfile, smatname);
  if (srwfile) fclose_maybe_compressed(srwfile, srwname);

  if (skip) {
    printf("%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
        DW, DW+W, dmatname, 100.0 * (double) DW / (DW+W));
    fflush(stdout);
    fclose_maybe_compressed(dmatfile, dmatname);
    fclose_maybe_compressed(drwfile, drwname);
    
    dcwfile = fopen_maybe_compressed(dcwname, wmode);
    for(uint32_t j = 0; j < skip; j++){
      uint32_t x = weights[j];
      if (bin) {
        fwrite32_little(&x, 1, dcwfile);
      } else {
        fprintf(dcwfile, "%" PRIu32 "\n", x);
      }
    }
    fclose_maybe_compressed(dcwfile, dcwname);
  }

  scwfile = fopen_maybe_compressed(scwname, wmode);
  for(int j = skip; j < nc; j++){
    uint32_t x = weights[j];
    if (bin) {
      fwrite32_little(&x, 1, scwfile);
    } else {
      fprintf(scwfile, "%" PRIu32 "\n", x);
    }
  }
  fclose_maybe_compressed(scwfile, scwname);

  free(smatname);
  free(srwname);
  free(scwname);
  if (skip) {
    free(dmatname);
    free(drwname);
    free(dcwname);
  }
  free(weights);
}


static void
write_index(const char *indexname, index_data_t index_data,
        uint64_t nr, uint64_t nc)
{
  FILE *indexfile;
  indexfile = fopen_maybe_compressed(indexname, "w");
  fprintf(indexfile, "%" PRIu64 " %"PRIu64"\n", nr, nc);

  for (uint64_t i = 0; i < nr; ++i) {
    ASSERT (index_data[i].n > 0);
    fprintf(indexfile, "%d", index_data[i].n);
    for (unsigned int j = 0; j < index_data[i].n; ++j) {
#ifdef FOR_DL
      fprintf(indexfile, " %" PRIx32 ":%d",
          index_data[i].rels[j].ind_row,
          index_data[i].rels[j].e);
#else
      ASSERT (index_data[i].rels[j].e == 1);
      fprintf(indexfile, " %" PRIx32 "",
          index_data[i].rels[j].ind_row);
#endif
    }
    fprintf(indexfile, "\n");
  }
  fclose_maybe_compressed(indexfile, indexname);
}

/********************/
/* Main and friends */
/********************/

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "purged", "input purged file");
  param_list_decl_usage(pl, "his", "input history file");
  param_list_decl_usage(pl, "split_point", "point in history where to split matrix.\n        By default this is 0, meaning that it splits when the first 3-merge\n" "        occurs. Passing -1 means no splitting at all.");
  param_list_decl_usage(pl, "out", "basename for output matrices");
  param_list_decl_usage(pl, "skip", "number of heaviest columns that go to the "
      "dense matrix (default " STR(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "index", "file containing description of rows "
      "(relations-sets) of the matrix");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage(param_list pl, char *argv0)
{
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  double cpu0 = seconds ();
  double wct0 = wct_seconds ();

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  /* Handle parameters on command line */
  char * argv0 = argv[0];
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "force-posix-threads",
      &filter_rels_force_posix_threads);
  if (argc == 0)
    usage (pl, argv0);
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unknown option: %s\n", argv[0]);
    usage(pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char * purgedname = param_list_lookup_string(pl, "purged");
  const char * hisname = param_list_lookup_string(pl, "his");
  const char * sparsename = param_list_lookup_string(pl, "out");
  const char * indexname = param_list_lookup_string(pl, "index");
  int skip = DEFAULT_MERGE_SKIP;
  param_list_parse_int(pl, "skip", &skip);
  const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
  set_antebuffer_path (argv0, path_antebuffer);
  int64_t split_point = 0;
  param_list_parse_int64(pl, "split_point", &split_point);

  /* Some checks on command line arguments */
  if (param_list_warn_unused(pl)) {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }
  if (purgedname == NULL) {
    fprintf(stderr, "Error, missing -purged command line argument\n");
    usage(pl, argv0);
  }
  if (sparsename == NULL) {
    fprintf(stderr, "Error, missing -out command line argument\n");
    usage(pl, argv0);
  }
  if (hisname == NULL) {
    fprintf(stderr, "Error, missing -his command line argument\n");
    usage(pl, argv0);
  }
  int zip = 0;
  if (has_suffix(sparsename, ".gz"))
    zip=1;
  int bin = 1;
  if (has_suffix(sparsename, ".txt") || has_suffix(sparsename, ".txt.gz")) {
    bin = 0;
  }
  if (bin) {
    printf ("# Output matrices will be written in binary format\n");
  } else {
    printf ("# Output matrices will be written in text format\n");
  }
  char * basename = strdup(sparsename);
  ASSERT_ALWAYS(basename != NULL);
  if (has_suffix(basename, ".gz"))
    basename[strlen(basename)-3]='\0';
  if (has_suffix(basename, ".bin") || has_suffix(basename, ".txt"))
    basename[strlen(basename)-4]='\0';
  if (has_suffix(basename, ".sparse"))
    basename[strlen(basename)-7]='\0';
  /* End of handling of parameters */


  /* Real work starts here */


  /* Dimensions of the purged matrix (the one that comes from purge).
   * They are read from the first line of the purged file.
   * We won't modify them thereafter. */
  uint64_t nr, nc;
  purgedfile_read_firstline(purgedname, &nr, &nc);
  printf("Purged matrix has %" PRIu64 " rows and %" PRIu64 " cols\n", nr, nc);

  /* The main matrix we'll work with: starts as the purged matrix, and at
   * the end, contains the sparse matrix.
   * We allocate it first and read it initialize it with the content of
   * the purged file. */
  matrix_t MM;
  MM = (matrix_t) malloc(nr*sizeof(typerow_t *));
  ASSERT_ALWAYS(MM != NULL);
  read_purgedfile(MM, purgedname, nr, nc);


  /* Pre-allocate M1 and M2 with appropriate sizes.
   *   M1 is an nr x nr matrix.
   *   M2 is an nr x nc matrix.
   * Also initialize M1 with the identity matrix and M2 to zero. */
  matrix_t M1, M2;
  M1 = (matrix_t) malloc(nr*sizeof(typerow_t *));
  M2 = (matrix_t) malloc(nr*sizeof(typerow_t *));
  ASSERT_ALWAYS(M1 != NULL && M2 != NULL);
  for(uint64_t i = 0; i < nr; ++i) {
    M1[i] = (typerow_t*) malloc(2*sizeof(typerow_t));
    ASSERT_ALWAYS(M1[i] != NULL);
    rowLength(M1, i) = 1;
    setCell(M1[i], 1, i, 1);
  }
  for(uint64_t i = 0; i < nr; ++i) {
    M2[i] = NULL;
  }

  /* If required, pre-allocate index_data */
  index_data_t index_data = NULL;
  if (indexname != NULL)
    index_data = init_index_data(nr);

  /* Apply .his file to MM, and record M1 and M2 on the way. */
  apply_hisfile(MM, M1, M2, nr, hisname, index_data, (uint64_t) split_point);

  /* Compute column weights in the resulting MM, to deduce extractor
   * matrices (row weights are already known). */
  uint32_t * MM_col_wt = (uint32_t*) malloc(nc * sizeof(uint32_t));
  ASSERT_ALWAYS(MM_col_wt != NULL);
  compute_col_weights(MM, MM_col_wt, nr, nc);

  /* Delete columns of M2 that are empty in MM.
   * This corresponds to applying the X matrix of the description in the 
   * introduction. */
  for (uint64_t i = 0; i < nr; i++) {
    if (M2[i] == NULL)
      continue;
    uint32_t len = rowLength(M2, i);
    uint32_t l = 1;
    for (uint32_t k = 1; k <= len; k++) {
      if (MM_col_wt[rowCell(M2[i], k)] != 0) {
        rowCell(M2[i], l) = rowCell(M2[i], k);
        l++;
      }
    }
    rowLength(M2, i) = l-1;
    if (l == 1) {
      free(M2[i]);
      M2[i] = NULL;
    }
  }

  /* Delete rows of M1 that are empty in MM.
   * This corresponds to applying the Y matrix. */
  for (uint64_t i = 0; i < nr; ++i) {
    if (MM[i] == NULL) {
      free(M1[i]);
      M1[i] = NULL;
    }
  }

  /* Delete relsets of index_data that are empty in MM. */
  if (index_data) {
    uint64_t j = 0;
    for (uint64_t i = 0; i < nr; ++i) {
      if (MM[i] != NULL) {
        ASSERT(index_data[i].n != 0);
        index_data[j] = index_data[i];
        if (j < i) {
          index_data[i].n = 0;
          index_data[i].rels = NULL;
        }
        j++;
      } else {
        ASSERT(index_data[i].n == 0);
      }
    }
  }

  /* Find sizes */
  uint64_t final_nr = 0;
  uint64_t final_nc = 0;
  final_nr = count_non_empty_rows(M1, nr);
  ASSERT(final_nr == count_non_empty_rows(MM, nr));
  final_nc = count_non_empty_cols(MM, MM_col_wt, nc);
  printf("Final sizes of matrix are %" PRIu64 " rows and %" PRIu64 " cols\n",
      final_nr, final_nc);

  /* At this point:
        M1 is a (final_nr x nr) matrix
        M2 is a (nr x final_nc) matrix
        MM = M1 * M2 is a  (final_nr x final_nc)  matrix
     (at least if we don't count empty rows/cols coming from MM.)

     However, M2 still contains empty rows, and the corresponding columns
     in M1 can be deleted.
   */

  /* Delete columns of M1 that correspond to empty rows in M2. */
  for (uint64_t i = 0; i < nr; i++) {
    if (M1[i] == NULL)
      continue;
    uint32_t l = 1;
    for (uint32_t k = 1; k <= rowLength(M1, i); k++) {
      if (M2[rowCell(M1, i, k)] != NULL) {
        rowCell(M1, i, l) = rowCell(M1, i, k);
        l++;
      }
    }
    rowLength(M1, i) = l-1;
    if (l == 1) {
      free(M1[i]);
      M1[i] = NULL;
    }
  }

  /* Find the real row-size of M2 (== col-size of M1) */
  uint64_t nrp = count_non_empty_rows(M2, nr);
  printf("Size of internal vector between M1 and M2 is %" PRIu64 "\n", nrp);

  /* Move non-empty rows of M2 at the beginning, and change column
     numbering of M1 accodingly, so that we still have
        MM = M1 * M2
     Also update index_data (that recalls rel_sets) if necessary.
    */
  {
    uint64_t *col_renumber = (uint64_t *)malloc(nr*sizeof(uint64_t));
    uint64_t j = 0;
    for (uint64_t i = 0; i < nr; ++i) {
      col_renumber[i] = UINT64_MAX;
      if (M2[i] != NULL) {
        // M2
        ASSERT(j <= i);
        M2[j] = M2[i];
        if (j < i)
          M2[i] = NULL;
        // remember the renumbering
        col_renumber[i] = j;
        j++;
      }
    }

    for (uint64_t i = 0; i < nr; ++i) {
      if (M1[i] == NULL)
        continue;
      for (uint32_t k = 1; k <= rowLength(M1, i); k++) {
          uint64_t col = rowCell(M1, i, k);
          ASSERT(col_renumber[col] != UINT64_MAX);
          rowCell(M1, i, k) = col_renumber[col];
      }
    }
    free(col_renumber);
  }

  /* Change column numbering of M2 according to weights. */
  {
    // Compute the column weights of M2.
    uint32_t * M2_col_wt = (uint32_t*) malloc(nc * sizeof(uint32_t));
    ASSERT_ALWAYS(M2_col_wt != NULL);
    compute_col_weights(M2, M2_col_wt, nrp, nc);

    // Compute renumbering permutation
    uint64_t *col_renumber = renumber_decr_weights(M2_col_wt, nc);

    // Apply it in place to M2
    for (uint64_t i = 0; i < nrp; ++i) {
      if (M2[i] == NULL)
        continue;
      for (uint32_t k = 1; k <= rowLength(M2, i); k++) {
        uint64_t col = rowCell(M2, i, k);
        ASSERT(col_renumber[col] != UINT64_MAX);
        rowCell(M2, i, k) = col_renumber[col];
      }
    }
    free(col_renumber);
    free(M2_col_wt);
  }
  
  /* Move non-empty rows of M1 at the beginning. */
  {
    uint64_t j = 0;
    for (uint64_t i = 0; i < nr; ++i) {
      if (M1[i] != NULL) {
        ASSERT(j <= i);
        M1[j] = M1[i];
        if (j < i)
          M1[i] = NULL;
        j++;
      }
    }
  }

  char * basenameM1 = (char *) malloc(strlen(basename)+4);
  ASSERT_ALWAYS(basenameM1 != NULL);
  strcpy(basenameM1, basename);
  strcat(basenameM1, ".M1");

  char * basenameM2 = (char *) malloc(strlen(basename)+4);
  ASSERT_ALWAYS(basenameM2 != NULL);
  strcpy(basenameM2, basename);
  strcat(basenameM2, ".M2");

  write_matrix(basenameM1, M1, final_nr, nrp, 0, bin, zip);
  write_matrix(basenameM2, M2, nrp, final_nc, skip, bin, zip);
  free(basename);
  free(basenameM1);
  free(basenameM2);

  printf("Total weight of M1: %lu\n", total_weight(M1, final_nr, 0));
  printf("Total weight of M2: %lu\n", total_weight(M2, nrp, skip));

  if (indexname != NULL) {
    write_index(indexname, index_data, final_nr, final_nc);
    for (uint64_t i = 0; i < nr; ++i)
      free(index_data[i].rels);
    free(index_data);
  }

  for (uint64_t i = 0; i < nr; ++i) {
    free(M1[i]);
    free(M2[i]);
    free(MM[i]);
  }
  free(M1);
  free(M2);
  free(MM);
  free(MM_col_wt);

  param_list_clear(pl);
  print_timing_and_memory (stdout, cpu0, wct0);
  return 0;
}
