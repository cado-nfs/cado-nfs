/* replay --- replaying history of merges to build the sparse matrix

Copyright 2008-2019 Francois Morain, Emmanuel Thome, Paul Zimmermann,
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

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MINGW
#include <fcntl.h>   /* for _O_BINARY */
#endif
#include <string.h>
#include <inttypes.h>        // for PRIu64, PRIu32, PRIx64
#include <stdint.h>          // for uint64_t, uint32_t, UINT32_MAX
#include "purgedfile.h"      // for purgedfile_read_firstline
#include "typedefs.h"        // for index_t, ideal_merge_t, index_signed_t
#include "filter_config.h"
#include "filter_io.h"  // earlyparsed_relation_ptr
#include "fix-endianness.h" // fwrite32_little
#include "gzip.h"       // fopen_maybe_compressed
#include "misc.h"       // derived_filename
#include "params.h"     // param_list_parse_*
#include "sparse.h"
#include "stats.h"      // stats_data_t
#include "timing.h"     // seconds
#include "verbose.h"    // verbose_decl_usage
#include "portability.h" // strdup  // IWYU pragma: keep
#include "macros.h"
#include "merge_replay_matrix.h"

#define DEBUG 0

// newrows[i] contains a new row formed of old rows that correspond to
// true original relations (and not multirelations).
//
// After computing newrows, we deduce for all old rows the list of newrows
// containing it.

static unsigned long
flushSparse(const char *sparsename, typerow_t **sparsemat, index_t small_nrows,
            index_t small_ncols, const index_t * code, index_t skip, int bin)
{
#ifdef FOR_DL
  ASSERT_ALWAYS (skip == 0);
#endif
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
          .smat = "txt",
          .srw = "rw.txt",
          .scw = "cw.txt",
          .dmat = "dense.txt",
          .drw = "dense.rw.txt",
          .dcw = "dense.cw.txt",
        },
        {
          .ext = ".bin",
          .smat = "bin",
          .srw = "rw.bin",
          .scw = "cw.bin",
          .dmat = "dense.bin",
          .drw = "dense.rw.bin",
          .dcw = "dense.cw.bin",
        },
    }, * suf = &(suffixes[bin]);

    unsigned long W = 0;
    unsigned long DW = 0;
    char * zip = has_suffix(sparsename, ".gz") ? ".gz" : NULL;
    index_t * weights = malloc(small_ncols * sizeof(index_t));
    ASSERT_ALWAYS(weights != NULL);
    memset(weights, 0, small_ncols * sizeof(index_t));

    char wmode[3] = "w";
#ifdef HAVE_MINGW
    if (bin)
      strcpy (wmode, "wb");
#endif

    char * base = strdup(sparsename);
    if (zip) { base[strlen(base)-3]='\0'; }
    if (has_suffix(base, suf->ext)) { base[strlen(base)-4]='\0'; }

    char * smatname = NULL;
    FILE * smatfile = NULL;
    char * srwname  = NULL;
    FILE * srwfile  = NULL;
    char * scwname  = NULL;
    FILE * scwfile  = NULL;

    {
        smatname = derived_filename(base, suf->smat, zip);
        srwname  = derived_filename(base, suf->srw, zip);
        scwname = derived_filename(base, suf->scw, zip);
        smatfile = fopen_maybe_compressed(smatname, wmode);
        srwfile  = fopen_maybe_compressed(srwname, wmode);
        /* XXX sm-outside-matrix creates a square matrix here */
        if (!bin) fprintf(smatfile, "%" PRIu64 " %" PRIu64 "\n",
                          (uint64_t) small_nrows,
                          (uint64_t) small_ncols - skip);
    }

    char * dmatname = NULL;
    FILE * dmatfile = NULL;
    char * drwname  = NULL;
    FILE * drwfile  = NULL;
    char * dcwname  = NULL;
    FILE * dcwfile  = NULL;

    if (skip) {
        /* arrange so that we don't get file names like .sparse.dense */
        char * dbase = strdup(base);
        char * tmp = strstr(dbase, ".sparse");
        if (tmp) {
            memmove(tmp, tmp + 7, strlen(tmp + 7) + 1);
        }
        dmatname = derived_filename(dbase, suf->dmat, zip);
        drwname  = derived_filename(dbase, suf->drw, zip);
        dcwname = derived_filename(dbase, suf->dcw, zip);
        dmatfile = fopen_maybe_compressed(dmatname, wmode);
        drwfile  = fopen_maybe_compressed(drwname, wmode);
        if (!bin)
            fprintf(dmatfile, "%" PRIu64 " %" PRIu64 "\n",
                    (uint64_t) small_nrows, (uint64_t) skip);
        free(dbase);
    }

    for (index_t i = 0; i < small_nrows; i++){
//#ifdef FOR_DL
#if 0
        /* this is for sm-outside-matrix */
        if (i == small_ncols) {
            printf("Rotating file names\n");
            char * base2;
            int rc = asprintf(&base2, "%s.tail", base);
            ASSERT_ALWAYS(rc >= 0);
            fclose_maybe_compressed (smatfile, smatname);
            free(smatname);
            smatname = derived_filename(base2, suf->smat, zip);
            smatfile = fopen_maybe_compressed(smatname, wmode);
            free(base2);
            if (!bin) {
                fprintf(smatfile, "%d %d\n", small_nrows - small_ncols, small_ncols - skip);
            }
            fclose_maybe_compressed(srwfile, srwname); srwfile = NULL;
        }
#endif

	if(sparsemat[i] == NULL) {
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
	    for(unsigned int j = 1; j <= rowLength(sparsemat, i); j++){
		if (code[rowCell(sparsemat[i], j)]-1 < skip) {
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


#if 0 //#ifdef FOR_DL //for DL sort the index
      for (int k = 1; k < rowLength(sparsemat, i); k++)
        {
          for (int l = k+1; l <= rowLength(sparsemat, i); l++)
            {
              index_t x = code[rowCell(sparsemat[i], k)]-1;
              index_t y = code[rowCell(sparsemat[i], l)]-1;
              if (x > y)
                {
                  typerow_t tmp = sparsemat[i][k];
                  sparsemat[i][k] = sparsemat[i][l];
                  sparsemat[i][l] = tmp;
                }
            }
        }
#endif

	    for(unsigned int j = 1; j <= rowLength(sparsemat, i); j++){
#if DEBUG >= 1
		ASSERT(code[rowCell(sparsemat[i], j)] > 0);
#endif
		ASSERT_ALWAYS(code[rowCell(sparsemat[i], j)]-1 <= (index_t) UINT32_MAX);
		uint32_t x = code[rowCell(sparsemat[i], j)]-1;
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
			/* exponents are always int32_t */
                        uint32_t e = (uint32_t) sparsemat[i][j].e;
                        fwrite32_little(&e, 1, smatfile);
#endif
                    } else {
		        fprintf(smatfile, " %" PRIu32 "", x);
#ifdef FOR_DL
                        fprintf(smatfile, ":%d", sparsemat[i][j].e);
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

    {
        fclose_maybe_compressed (smatfile, smatname);
        if (srwfile) fclose_maybe_compressed (srwfile, srwname);
    }
    if (skip) {
        printf("%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
                DW, DW+W, dmatname,
                100.0 * (double) DW / (DW+W+(DW==0&&W==0)));
        fflush(stdout);

        fclose_maybe_compressed (dmatfile, dmatname);
        fclose_maybe_compressed (drwfile, drwname);
    }

    if (skip) {
        dcwfile = fopen_maybe_compressed(dcwname, wmode);
        for (unsigned int j = 0; j < skip; j++){
	    ASSERT_ALWAYS(weights[j] <= (index_t) UINT32_MAX);
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, dcwfile);
            } else {
                fprintf(dcwfile, "%" PRIu32 "\n", x);
            }
        }
        fclose_maybe_compressed(dcwfile, dcwname);
    }

    {
        scwfile = fopen_maybe_compressed(scwname, wmode);
        for(index_t j = skip; j < small_ncols; j++){
	    ASSERT_ALWAYS(weights[j] <= (index_t) UINT32_MAX);
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, scwfile);
            } else {
                fprintf(scwfile, "%" PRIu32 "\n", x);
            }
        }
        fclose_maybe_compressed(scwfile, scwname);
    }

    {
        free(smatname);
        free(srwname);
        free(scwname);
    }
    if (skip) {
        free(dmatname);
        free(drwname);
        free(dcwname);
    }

    free(base);
    free(weights);

    return W;
}


// on input, colweight[j] contains the weight; on exit, colweight[j]
// contains the new index for j. Heavier columns are in front of the new
// matrix.
static void
renumber (index_t small_ncols, index_t *colweight, index_t ncols,
          MAYBE_UNUSED const char *idealsfilename)
{
    index_t k, nb;
    index_t *tmp;

#ifdef FOR_DL
    FILE *renumberfile = fopen_maybe_compressed (idealsfilename, "w");
    if (renumberfile == NULL)
    {
      fprintf (stderr, "Error while opening file to save permutation of"
                       "ideals\n");
      exit(EXIT_FAILURE);
    }
#endif

    ASSERT_ALWAYS(small_ncols);
    tmp = (index_t *) malloc ((small_ncols << 1) * sizeof(index_t));
    ASSERT_ALWAYS(tmp != NULL);
    memset(tmp, 0, (small_ncols << 1) * sizeof(index_t));
    for (k = 0, nb = 0; k < ncols; k++)
      if (colweight[k] > 0)
	{
	  tmp[nb++] = colweight[k];
	  tmp[nb++] = k;
	}
    ASSERT_ALWAYS(nb == 2 * small_ncols);
#ifdef FOR_DL
    fprintf (renumberfile, "# %" PRIu64 "\n", (uint64_t) small_ncols);
#endif
    printf ("Sorting %" PRIu64 " columns by decreasing weight\n",
            (uint64_t) small_ncols);
    fflush (stdout);
    qsort (tmp, small_ncols, 2*sizeof(index_t), cmp_index2);
    memset (colweight, 0, ncols * sizeof(index_t));
    // useful for BW + skipping heavy part only...
    index_signed_t j;
    ASSERT_ALWAYS(nb >= 1);
    for (j = (index_signed_t) (nb - 1), k = 1; j >= 0; j -= 2)
      {
        colweight[tmp[j]] = k++; // always this +1 trick
#ifdef FOR_DL
        fprintf (renumberfile, "%" PRIu64 " %" PRIx64 "\n",
                 (uint64_t) colweight[tmp[j]]-1, (uint64_t) tmp[j]);
#endif
      }

#ifdef FOR_DL
    fclose_maybe_compressed (renumberfile, idealsfilename);
#endif
    free(tmp);
}

// A line is "i i1 ... ik [#j]".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also if i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
// If given, j is the index of the column used for pivoting (used in DL).
// If newrows=NULL, we only compute the index.
static void
doAllAdds(typerow_t **newrows, char *str, index_data_t index_data)
{
  index_t j;
  index_signed_t ind[MERGE_LEVEL_MAX], i0;
  int ni, destroy;

  ni = parse_hisfile_line (ind, str, &j);

  if (ind[0] < 0)
    {
      destroy = 0;
      i0 = -ind[0]-1;
    }
  else
    {
      destroy = 1;
      i0 = ind[0];
    }
#if DEBUG >= 1
    fprintf(stderr, "first i is %d in %s", i0, str);
#endif

  for (int k = 1; k < ni; k++)
    addRowsUpdateIndex(newrows, index_data, ind[k], i0, j);

  if(destroy)
    {
      //destroy initial row!
      if (newrows != NULL)
        {
          free(newrows[i0]);
          newrows[i0] = NULL;
        }

      if (index_data != NULL) // ie we want an index
      {
        index_data[i0].n = 0;
        free(index_data[i0].rels);
        index_data[i0].rels = NULL;
      }
    }
}


#define STRLENMAX 2048

// sparsemat is small_nrows x small_ncols
static void
toFlush (const char *sparsename, typerow_t **sparsemat, index_t *colweight,
         index_t ncols, index_t small_nrows, index_t small_ncols, int skip,
         int bin, const char *idealsfilename, unsigned int nsquare_matrices)

{
    unsigned long W;

    printf("Renumbering columns (including sorting w.r.t. weight)\n");
    fflush(stdout);
    renumber (small_ncols, colweight, ncols, idealsfilename);

    printf ("Sparse submatrix: nrows=%" PRIu64 " ncols=%" PRIu64 "\n",
            (uint64_t) small_nrows, (uint64_t) small_ncols - skip);
    ASSERT_ALWAYS(small_nrows >= small_ncols);

    if (nsquare_matrices)
    {
        ASSERT_ALWAYS(nsquare_matrices <= (small_nrows-small_ncols+1));
        printf ("Will output %u matri%s of size %" PRIu64 "x%" PRIu64 "\n",
                nsquare_matrices, nsquare_matrices == 1 ? "x" : "ces",
                (uint64_t) small_ncols, (uint64_t) small_ncols);
    }

    double tt = seconds();
    printf("Writing sparse representation to %s\n", sparsename);
    fflush(stdout);
    W = flushSparse(sparsename, sparsemat,
                    nsquare_matrices ? small_ncols : small_nrows,
                    small_ncols, colweight, skip, bin);
    printf("# Writing matrix took %.1lfs\n", seconds()-tt);
    printf("# Weight of the sparse submatrix: %lu\n", W);
    fflush(stdout);

    if (nsquare_matrices > 1)
    {
        int zip = has_suffix(sparsename, ".gz");
        const char *ext = bin ? (zip ? ".bin.gz" : ".bin")
                              : (zip ? ".txt.gz" : ".txt");
        char istr[] = { '1', '\0' };

        for (unsigned int i = 1; i < nsquare_matrices; i++)
        {
            double tt = seconds();
            char * smatname = derived_filename(sparsename, istr, ext);
            printf("Writing sparse representation to %s\n", smatname);
            fflush(stdout);
            W = flushSparse(smatname, sparsemat+i, small_ncols, small_ncols,
                            colweight, skip, bin);
            istr[0]++;
            printf("# Writing matrix took %.1lfs\n", seconds()-tt);
            printf("# Weight of the sparse submatrix: %lu\n", W);
            fflush(stdout);
            free(smatname);
        }
    }
}

/* If newrows = NULL, only compute the index */
static void
build_newrows_from_file(typerow_t **newrows, FILE *hisfile,
                        index_data_t index_data, index_t nrows, index_t Nmax)
{
    uint64_t addread = 0;
    char str[STRLENMAX];

    printf("Reading row additions\n");
    fflush(stdout);
    stats_data_t stats; /* struct for printing progress */
    /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
     * 2^23 primes after that */
    stats_init (stats, stdout, &addread, 23, "Read", "row additions", "", "lines");
    while(fgets(str, STRLENMAX, hisfile) && nrows >= Nmax)
    {
        if (str[0] == '#') continue;

      addread++;

      if (stats_test_progress(stats))
        stats_print_progress (stats, addread, 0, 0, 0);

      if(str[strlen(str)-1] != '\n')
      {
        fprintf(stderr, "Gasp: not a complete a line!");
        fprintf(stderr, " I stop reading and go to the next phase\n");
        break;
      }
      doAllAdds(newrows, str, index_data);
      /* if the line starts with a negative index, the number of rows remains
         unchanged, otherwise it decreases by one */
      nrows -= str[0] != '-';
    }
    stats_print_progress (stats, addread, 0, 0, 1);
}

typedef struct
{
  typerow_t **mat;
  index_t ncols;
  index_t col0;
  index_t colmax;
} replay_read_data_t;

void * fill_in_rows (void *context_data, earlyparsed_relation_ptr rel)
{
  replay_read_data_t *data = (replay_read_data_t *) context_data;
  typerow_t buf[UMAX(weight_t)];

  unsigned int nb = 0;
  for (unsigned int j = 0; j < rel->nb; j++)
  {
    index_t h = rel->primes[j].h;
    if (h < data->col0 || h >= data->colmax) continue;
    nb++;
#ifdef FOR_DL
    exponent_t e = rel->primes[j].e;
    buf[nb] = (ideal_merge_t) {.id = h, .e = e};
#else
    ASSERT_ALWAYS (rel->primes[j].e == 1);
    buf[nb] = h;
#endif
    ASSERT (h < data->ncols);
  }
#ifdef FOR_DL
  buf[0].id = nb;
#else
  buf[0] = nb;
#endif

  qsort (&(buf[1]), nb, sizeof(typerow_t), cmp_typerow_t);

  data->mat[rel->num] = mallocRow (nb + 1);
  compressRow (data->mat[rel->num], buf, nb);

  return NULL;
}

/* if for_msieve=1, generate the *.cyc file needed by msieve to construct
   its matrix, which is of the following (binary) format:
      small_nrows
      n1 i1 i2 ... in1
      n2 j1 j2 ... jn2
      ...
      nk ...
   where each value is stored as a 32-bit integer (no linebreak),
   small_nrows is the number of relation-sets of the matrix,
   n1 is the number of relations in the first relation-set,
   i1 is the index of the first relation in the first relation-set
   (should correspond to line i1+2 in *.purged.gz), and so on */

static void
read_purgedfile (typerow_t **mat, const char* filename, index_t nrows,
                 index_t ncols, index_t col0, index_t colmax, int for_msieve)
{
  index_t nread;
  if (for_msieve == 0)
  {
    printf("Reading sparse matrix from %s\n", filename);
    fflush(stdout);
    const char * fic[2] = {filename, NULL};
    replay_read_data_t tmp = (replay_read_data_t) {
        .mat= mat,
        .ncols = ncols,
        .col0 = col0,
        .colmax = colmax,
    };
    nread = filter_rels(fic, (filter_rels_callback_t) &fill_in_rows, &tmp,
                        EARLYPARSE_NEED_INDEX, NULL, NULL);
    ASSERT_ALWAYS (nread == nrows);
  }
  else /* for_msieve */
  {
    /* to generate the .cyc file for msieve, we only need to start from the
       identity matrix with newnrows relation-sets, where relation-set i
       contains only relation i. Thus we only need to read the first line of
       the purged file, to get the number of relations-sets. */

    for (index_t i = 0; i < nrows; i++)
    {
      mat[i] = (typerow_t *) malloc(2 * sizeof(typerow_t));
      setCell(mat[i], 1, i, 1);
      setCell(mat[i], 0, 1, 1);
    }
  }
}

static void
writeIndex(const char *indexname, index_data_t index_data, index_t small_nrows)
{
    FILE *indexfile = NULL;
    indexfile = fopen_maybe_compressed(indexname, "w");
    ASSERT_ALWAYS (indexfile != NULL);
    fprintf(indexfile, "%" PRIu64 "\n", (uint64_t) small_nrows);

    for (index_t i = 0; i < small_nrows; ++i) {
        ASSERT (index_data[i].n > 0);
        fprintf(indexfile, "%d", index_data[i].n);
        for (unsigned int j = 0; j < index_data[i].n; ++j) {
#ifdef FOR_DL
            fprintf(indexfile, " %" PRIx64 ":%d",
                    (uint64_t) index_data[i].rels[j].ind_row,
                    index_data[i].rels[j].e);
#else
            fprintf(indexfile, " %" PRIx64 "",
                    (uint64_t) index_data[i].rels[j].ind_row);
#endif
        }
        fprintf(indexfile, "\n");
    }
    fclose_maybe_compressed(indexfile, indexname);
}

void
generate_cyc (const char *outname, typerow_t **rows, index_t nrows)
{
  FILE *outfile;
  index_t t, u, i, k;

  outfile = fopen (outname, "w");
  ASSERT_ALWAYS(outfile != NULL);

  /* first write the number of relations */
  t = nrows;
  fwrite (&t, sizeof(index_t), 1, outfile);

  /* then for each relation-set write a 32-bit integer giving the number k
     of its element, followed by those k elements */
  for (i = 0; i < nrows; i++)
    {
      t = rowLength(rows, i);
      ASSERT_ALWAYS(t != 0);
      fwrite (&t, sizeof(uint32_t), 1, outfile);
      for (k = 1; k <= t; k++)
        {
          u = rowCell(rows[i], k);
          fwrite (&u, sizeof(uint32_t), 1, outfile);
        }
    }

  fclose (outfile);
}

static void
fasterVersion (typerow_t **newrows, const char *sparsename,
               const char *indexname, const char *hisname, index_t nrows,
               index_t ncols, int skip, int bin, const char *idealsfilename,
               int for_msieve, index_t Nmax, unsigned int nsquare_matrices)
{
  FILE *hisfile = NULL;
  index_t *colweight = NULL;
  index_t small_nrows = 0, small_ncols;
  index_data_t index_data = NULL;

  hisfile = fopen_maybe_compressed (hisname, "r");
  ASSERT_ALWAYS(hisfile != NULL);

  if (indexname != NULL)
  {
    /* At the beginning, the index_data consists of relsets that are just single
     * relations.*/
    index_data = (relset_t *) malloc (nrows * sizeof(relset_t));
    ASSERT_ALWAYS (index_data != NULL);
    for (index_t i = 0; i < nrows; ++i)
    {
      index_data[i].n = 1;
      index_data[i].rels = (multirel_t *) malloc (sizeof(multirel_t));
      ASSERT_ALWAYS (index_data[i].rels != NULL);
      index_data[i].rels[0].ind_row = i;
#ifdef FOR_DL
      index_data[i].rels[0].e = 1;
#endif
     }
  }
  else
    index_data = NULL;

  /* read merges in the *.merge.his file and replay them */
  build_newrows_from_file (newrows, hisfile, index_data, nrows, Nmax);

  if (sparsename != NULL)
    {
      /* crunch empty rows first to save memory and compute small_nrows */
      index_t j = 0;
      for (index_t i = 0; i < nrows; i++)
        if (newrows[i] != NULL)
          newrows[j++] = newrows[i]; /* we always have j <= i */
      small_nrows = j;
      CHECKED_REALLOC(newrows, small_nrows, typerow_t *);
    }

  /* if index was asked: crunch the empty rows as above, create the index and
   * free index_data before calling toFlush(), in order to decrease the total
   * memory usage */
  if (indexname != NULL)
  {
    index_t j = 0;
    for (index_t i = 0; i < nrows; i++)
    {
      if (index_data[i].n > 0)
        index_data[j++] = index_data[i];
      else
        free(index_data[i].rels);
    }
    if (sparsename != NULL)
      ASSERT_ALWAYS (j == small_nrows);
    else
      small_nrows = j;

    writeIndex (indexname, index_data, small_nrows);

    for (index_t i = 0; i < small_nrows; ++i)
      free (index_data[i].rels);
    free (index_data);
  }

  if (sparsename == NULL)
    goto close_and_exit;

  /* compute column weights */
  colweight = (index_t*) malloc (ncols * sizeof(index_t));
  ASSERT_ALWAYS (colweight != NULL);
  memset (colweight, 0, ncols * sizeof(index_t));
  for (index_t i = small_ncols = 0; i < small_nrows; i++)
    for(unsigned int k = 1; k <= rowLength(newrows, i); k++)
    {
      index_t j = rowCell(newrows[i], k);
      small_ncols += (colweight[j] == 0);
      colweight[j] ++;
    }
  /* small_ncols is the number of columns with non-empty weight,
     i.e., the number of columns of the final matrix */

#if defined FOR_DL && defined STAT_DL
  index_t count[11] = {0,0,0,0,0,0,0,0,0,0,0};
  index_t nonzero = 0;
  for (index_t i = 0; i < small_nrows ; i++)
    {
      for(unsigned int k = 1; k <= rowLength(newrows, i); k++)
        {
          if (abs(newrows[i][k].e) > 10)
            count[0]++;
          else
            count[abs(newrows[i][k].e)]++;
          nonzero++;
        }
    }
  fprintf (stderr, "# of non zero coeff: %lu\n", nonzero);
  for (int i = 1; i <= 10 ; i++)
    fprintf (stderr, "# of %d: %lu(%.2f%%)\n", i, count[i],
                     100 * (double) count[i]/nonzero);
  fprintf (stderr, "# of > 10: %lu(%.2f%%)\n", count[0],
                   100 * (double) count[0]/nonzero);
#endif

  if (for_msieve)
    {
      /* generate the <dat_file_name>.cyc file in "indexname" */
      if (skip != 0)
        {
          fprintf (stderr, "Error, skip should be 0 with --for_msieve\n");
          exit (1);
        }
      if (bin != 0)
        {
          fprintf (stderr, "Error, --binary incompatible with --for_msieve\n");
          exit (1);
        }
      generate_cyc (sparsename, newrows, small_nrows);
    }
  else if (sparsename != NULL)
    /* renumber columns after sorting them by decreasing weight */
    toFlush (sparsename, newrows, colweight, ncols, small_nrows, small_ncols,
             skip, bin, idealsfilename, nsquare_matrices);

  /* Free */
  free (colweight);
  for (index_t i = 0; i < small_nrows; i++)
    free (newrows[i]);
  free (newrows);

 close_and_exit:
  fclose_maybe_compressed (hisfile, hisname);
}

/*
static void
usage (const char *argv0)
{
  fprintf (stderr, "   --binary\n");
  fprintf (stderr, "   -ideals\n");
  exit (1);
}
*/

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "purged", "input purged file");
  param_list_decl_usage(pl, "his", "input history file");
  param_list_decl_usage(pl, "out", "basename for output matrices");
#ifndef FOR_DL
  param_list_decl_usage(pl, "skip", "number of heaviest columns that go to the "
                            "dense matrix (default " CADO_STRINGIZE(DEFAULT_MERGE_SKIP) ")");
#endif
  param_list_decl_usage(pl, "index", "file containing description of rows "
                                     "(relations-sets) of the matrix");
  param_list_decl_usage(pl, "ideals", "file containing correspondence between "
                                      "ideals and matrix columns");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "for_msieve", "output matrix in msieve format");
  param_list_decl_usage(pl, "Nmax", "stop at Nmax number of rows (default 0)");
#ifndef FOR_DL
  param_list_decl_usage(pl, "col0", "print only columns with index >= col0");
  param_list_decl_usage(pl, "colmax", "print only columns with index < colmax");
#endif
  param_list_decl_usage(pl, "nsquare_matrices", "if non-zero, output that many "
                                "square matrices by skipping rows (default 0)");

  verbose_decl_usage(pl);
}

static void
usage (param_list pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

// We start from M_purged which is nrows x ncols;
// we build M_small which is small_nrows x small_ncols.
// newrows[i] if != NULL, contains a list of the indices of the rows in
// M_purged that were added together to form this new row in M_small.
// TODO: replace this index by the index to rels directly to skip one
// indirection???
int main(int argc, char const * argv[])
{
  const char * argv0 = argv[0];
  uint64_t Nmax = 0;
  uint64_t nrows, ncols;
  typerow_t **newrows;
  int bin = -1, skip = DEFAULT_MERGE_SKIP, for_msieve = 0;
  double cpu0 = seconds ();
  double wct0 = wct_seconds ();

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    setvbuf(stderr, NULL, _IONBF, 0);
    setvbuf(stdout, NULL, _IONBF, 0);

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "for_msieve", &for_msieve);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

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
    const char * idealsfilename = param_list_lookup_string(pl, "ideals");
#ifndef FOR_DL
    param_list_parse_int(pl, "skip", &skip);
#endif
    param_list_parse_uint64(pl, "Nmax", &Nmax);
    const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    index_t col0 = 0;
    index_t colmax = UMAX(index_t);

#ifndef FOR_DL
    { uint64_t c; if (param_list_parse_uint64(pl, "col0", &c))   col0 = c; }
    { uint64_t c; if (param_list_parse_uint64(pl, "colmax", &c)) colmax = c; }
#endif

    unsigned int nsquare_matrices = 0;
    param_list_parse_uint(pl, "nsquare_matrices", &nsquare_matrices);
    if (nsquare_matrices)
    {
        printf ("# Will output %u square matri%s\n", nsquare_matrices,
                nsquare_matrices == 1 ? "x" : "ces" );
        ASSERT_ALWAYS (skip == 0); /* not compatible with skip */
        ASSERT_ALWAYS (nsquare_matrices < 10);
    }

    /* Some checks on command line arguments */
    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (purgedname == NULL)
    {
      fprintf(stderr, "Error, missing -purged command line argument\n");
      usage(pl, argv0);
    }
    if (hisname == NULL)
    {
      fprintf(stderr, "Error, missing -his command line argument\n");
      usage(pl, argv0);
    }
    if (sparsename == NULL && indexname == NULL)
      {
        fprintf(stderr, "Error, at least one of -out and -index is required\n");
        usage(pl, argv0);
      }
#ifdef FOR_DL
    if (idealsfilename == NULL)
    {
      fprintf(stderr, "Error, missing -ideals command line argument\n");
      usage(pl, argv0);
    }
    ASSERT_ALWAYS (skip == 0);
#endif
    if (sparsename != NULL)
      {
        if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
          {
            bin = 1;
            printf ("# Output matrices will be written in binary format\n");
          }
        else
          {
            bin = 0;
            printf ("# Output matrices will be written in text format\n");
          }
      }

    set_antebuffer_path (argv0, path_antebuffer);

  /* Read number of rows and cols on first line of purged file */
  purgedfile_read_firstline (purgedname, &nrows, &ncols);
  if (nrows >= 4294967296UL)
    {
      fprintf (stderr, "Error, cannot handle 2^32 rows or more after purge\n");
      fprintf (stderr, "change ind_row from uint32_t to uint64_t in sparse.h\n");
      exit (EXIT_FAILURE);
    }
  printf("Sparse matrix has %" PRIu64 " rows and %" PRIu64 " cols\n",
         nrows, ncols);
  fflush(stdout);

#if SIZEOF_INDEX == 4
  if (ncols >= UINT32_MAX) {
      fprintf(stderr, "You must recompile with -DSIZEOF_INDEX=8\n");
      exit(EXIT_FAILURE);
  }
#endif

  /* Allocate memory for rows of the matrix */
  if (sparsename != NULL)
    {
      newrows = (typerow_t **) malloc (nrows * sizeof(typerow_t *));
      ASSERT_ALWAYS(newrows != NULL);
    }
  else
    newrows = NULL;

  /* Read the matrix from purgedfile */
  if (sparsename != NULL)
    {
      read_purgedfile (newrows, purgedname, nrows, ncols, col0, colmax, for_msieve);
      printf("The biggest index appearing in a relation is %" PRIu64 "\n", ncols);
      fflush(stdout);
#if DEBUG >=1
      for (index_t i = 0; i < nrows; i++)
        {
          fprintf(stderr, "row[%" PRIu64 "] :", i);
          fprintRow(stderr, newrows[i]);
          fprintf(stderr, "\n");
        }
#endif
    }

  fasterVersion (newrows, sparsename, indexname, hisname, nrows, ncols, skip,
                 bin, idealsfilename, for_msieve, Nmax, nsquare_matrices);

  param_list_clear(pl);
  print_timing_and_memory (stdout, cpu0, wct0);
  return 0;
}
