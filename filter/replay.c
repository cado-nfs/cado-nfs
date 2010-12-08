#include "cado.h"
/* 
 * Program: replay
 * Author : F. Morain
 * Purpose: replaying history of merges to build the small matrix
 * 
 * Algorithm:
 *
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

#include "sparse.h"
#include "gzip.h"

#define DEBUG 0

#define TRACE_COL -1 // 231 // put to -1 if not...!

// newrows[i] contains a new row formed of old rows that correspond to
// true original relations (and not multirelations).
//
// After computing newrows, we deduce for all old rows the list of newrows
// containing it.

#if DEBUG >= 1
static void
printOldRows(int **oldrows, int nrows)
{
    int i;
    
    // where oldrows intervene
    for(i = 0; i < nrows; i++)
	if(oldrows[i][0] != 0){
	    fprintf(stderr, "old row_%d takes part to new row(s)", i);
	    fprintRow(stderr, oldrows[i]);
	    fprintf(stderr, "\n");
	}
}

static void
printBuf(FILE *file, int *buf, int ibuf)
{
    int i;

    fprintf(file, "buf[%d] =", ibuf);
    for(i = 0; i < ibuf; i++)
	fprintf(file, " %d", buf[i]);
}
#endif

// add buf[0..ibuf[ to row i of sparsemat
static void
addrel(int **sparsemat, int *colweight, int *buf, int ibuf, int i)
{
    int *rowi = sparsemat[i], *tmp, tmp_len, k1, k2, k;

    if(rowi == NULL){
	// first time, dump buf
	rowi = (int *)malloc((ibuf+1) * sizeof(int));
	rowi[0] = ibuf;
	memcpy(rowi+1, buf, ibuf * sizeof(int));
	sparsemat[i] = rowi;
	for(k = 1; k <= rowi[0]; k++){
#if TRACE_COL >= 0
	    if(rowi[k] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears 1st in R_%d\n",TRACE_COL,i);
#endif
	    colweight[rowi[k]]++;
	}
    }
    else{
#if DEBUG >= 1
	fprintRow(stderr, rowi); fprintf(stderr, "\n");
	printBuf(stderr, buf, ibuf); fprintf(stderr, "\n");
#endif
	// rowi is bound to disappear
	removeWeight(sparsemat, colweight, i);
	tmp_len = rowi[0] + ibuf + 1;
	tmp = (int *)malloc(tmp_len * sizeof(int));
	k1 = 1;
	k2 = 0;
	k = 1;
	while((k1 <= rowi[0]) && (k2 < ibuf)){
#if TRACE_COL >= 0
	    if(rowi[k1] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears in R_%d\n",TRACE_COL,i);
	    if(buf[k2] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears in buf_%d\n",TRACE_COL,i);
#endif
	    if(rowi[k1] < buf[k2]){
		tmp[k++] = rowi[k1++];
		colweight[rowi[k1-1]]++;
	    }
	    else if(rowi[k1] > buf[k2]){
		tmp[k++] = buf[k2++];
		colweight[buf[k2-1]]++;
	    }
	    else{
		k1++; k2++;
	    }
	}
	for(; k1 <= rowi[0]; k1++){
	    tmp[k++] = rowi[k1];
	    colweight[rowi[k1]]++;
	}
	for(; k2 < ibuf; k2++){
	    tmp[k++] = buf[k2];
	    colweight[buf[k2]]++;
	}
	tmp[0] = k-1;
	free(rowi);
	sparsemat[i] = tmp;
#if DEBUG >= 1
	fprintRow(stderr, tmp); fprintf(stderr, "\n");
	//	if(tmp[0] != (tmp_len-1))
	//	    printf("#W# shorter length\n"); // who cares, really?
#endif
	//	fprintf(stderr, "W: 231 -> %d\n", colweight[231]);
    }
}

static void
makeSparse(int **sparsemat, int *colweight, purgedfile_stream ps,
           int jmin, int jmax, int **oldrows, int verbose, int nslices)
{
    fprintf(stderr, "Reading and treating relations from purged file\n");

    for(int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
	if (verbose && purgedfile_stream_disp_progress_now_p(ps))
	    fprintf(stderr, "Treating old rel #%d at %2.2lf\n",
                    ps->rrows,ps->dt);

        // take only the primes within the requested interval.
        int ibuf = 0;
        for(int k = 0; k < ps->nc; k++) {
            int j = ps->cols[k];
            if (j >= jmin && j < jmax)
                ps->cols[ibuf++] = j;
        }
        ps->nc = ibuf;

	if(ibuf > 0){
	    qsort(ps->cols, ps->nc, sizeof(int), cmp);
	    // now, we "add" relation [a, b] in all new relations in which it
	    // participates
	    for(int k = 1; k <= oldrows[i][0]; k++)
		addrel(sparsemat, colweight, ps->cols, ps->nc, oldrows[i][k]);
	}
	if(nslices == 0){
	    // free space just in case
	    free(oldrows[i]);
	    oldrows[i] = NULL;
	}
    }
}

static unsigned long
flushSparse(const char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code, int skip, int bin)
{
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
    uint32_t * weights = malloc(small_ncols * sizeof(uint32_t));
    memset(weights, 0, small_ncols * sizeof(uint32_t));

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
        smatfile = gzip_open(smatname, "w");
        srwfile  = gzip_open(srwname, "w");
        if (!bin)
            fprintf(smatfile, "%d %d\n", small_nrows, small_ncols - skip);
    }

    char * dmatname = NULL;
    FILE * dmatfile = NULL;
    char * drwname  = NULL;
    FILE * drwfile  = NULL;
    char * dcwname  = NULL;
    FILE * dcwfile  = NULL;

    if (skip) {
        dmatname = derived_filename(base, suf->dmat, zip);
        drwname  = derived_filename(base, suf->drw, zip);
        dmatfile = gzip_open(dmatname, "w");
        drwfile  = gzip_open(drwname, "w");
        if (!bin)
            fprintf(dmatfile, "%d %d\n", small_nrows, skip);
    }

    for(int i = 0; i < small_nrows; i++){
	if(sparsemat[i] == NULL) {
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
	    for(int j = 1; j <= sparsemat[i][0]; j++){
		if (code[sparsemat[i][j]]-1 < skip) {
                    dw++;
                    DW++;
                } else {
                    sw++;
                    W++;
                }
            }
            if (bin) {
                fwrite32_little(&sw, 1, smatfile);
                fwrite32_little(&sw, 1, srwfile);
                if (skip) fwrite32_little(&dw, 1, dmatfile);
                if (skip) fwrite32_little(&dw, 1, drwfile);
            } else {
                fprintf(smatfile, "%"PRIu32"", sw);
                fprintf(srwfile, "%"PRIu32"\n", sw);
                if (skip) fprintf(dmatfile, "%"PRIu32"", dw);
                if (skip) fprintf(drwfile, "%"PRIu32"\n", dw);
            }
	    for(int j = 1; j <= sparsemat[i][0]; j++){
#if DEBUG >= 1
		ASSERT(code[sparsemat[i][j]] > 0);
#endif
		uint32_t x = code[sparsemat[i][j]]-1;
                weights[x]++;
                if ((int) x < skip) {
                    ASSERT_ALWAYS(skip);
                    if (bin) {
                        fwrite32_little(&x, 1, dmatfile);
                    } else {
                        fprintf(dmatfile, " %"PRIu32"", x);
                    }
                } else {
                    x-=skip;
                    if (bin) {
                        fwrite32_little(&x, 1, smatfile);
                    } else {
                        fprintf(smatfile, " %"PRIu32"", x);
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
        gzip_close (smatfile, smatname);
        gzip_close (srwfile, srwname);
        free(smatname);
        free(srwname);
    }
    if (skip) {
        fprintf(stderr, "%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
                DW, DW+W, dmatname,
                100.0 * (double) DW / (DW+W));

        gzip_close (dmatfile, dmatname);
        gzip_close (drwfile, drwname);
        free(dmatname);
        free(drwname);
    }

    if (skip) {
        dcwname = derived_filename(base, suf->dcw, zip);
        dcwfile = gzip_open(dcwname, "w");
        for(int j = 0; j < skip; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, dcwfile);
            } else {
                fprintf(dcwfile, "%"PRIu32"\n", x);
            }
        }
        gzip_close(dcwfile, dcwname);
        free(dcwname);
    }

    {
        scwname = derived_filename(base, suf->scw, zip);
        scwfile = gzip_open(scwname, "w");
        for(int j = skip; j < small_ncols; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, scwfile);
            } else {
                fprintf(scwfile, "%"PRIu32"\n", x);
            }
        }
        gzip_close(scwfile, scwname);
        free(scwname);
    }

    free(base);
    free(weights);

    return W;
}

// dump of newrows in indexname.
static void
makeIndexFile(const char *indexname, int nrows, int **newrows, int small_nrows, int small_ncols)
{
    FILE *indexfile;
    int i, j;

    indexfile = gzip_open(indexname, "w");
    fprintf(indexfile, "%d %d\n", small_nrows, small_ncols);
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    fprintf(indexfile, "%d", newrows[i][0]);
	    for(j = 1; j <= newrows[i][0]; j++){
		fprintf(indexfile, " ");
		fprintf(indexfile, PURGE_INT_FORMAT, newrows[i][j]);
	    }
	    fprintf(indexfile, "\n");
	}
    gzip_close(indexfile, indexname);
}

// on input, colweight[j] contains the weight; on exit, colweight[j]
// contains the new index for j.
static void
renumber(const char * sosname, int *small_ncols, int *colweight, int ncols)
{
    int j, k, nb, *tmp;

    tmp = (int *)malloc((ncols<<1) * sizeof(int));
    memset(tmp, 0, (ncols<<1) * sizeof(int));
    for(j = 0, nb = 0; j < ncols; j++)
	if(colweight[j] > 0){
#if DEBUG >= 1
	    fprintf(stderr, "J %d %d\n", j, colweight[j]);
#endif
	    tmp[nb++] = colweight[j];
	    tmp[nb++] = j;
	}
    *small_ncols = nb>>1;
    qsort(tmp, nb>>1, 2*sizeof(int), cmp);
    memset(colweight, 0, ncols * sizeof(int));
#if 0
    // useful for Gauss only...
    fprintf(stderr, "Sorting in INcreasing weight order of j\n");
    for(j = 0, k = 1; j < nb; j += 2)
	colweight[tmp[j+1]] = k++; // always this +1 trick
#else
    // useful for BW + skipping heavy part only...
    fprintf (stderr, "Sorting columns by decreasing weight\n");
    for(j = nb-1, k = 1; j >= 0; j -= 2)
	colweight[tmp[j]] = k++; // always this +1 trick
    if (sosname) {
        fprintf(stderr, "Saving ideal renumbering info to %s\n", sosname);
        FILE * f = fopen(sosname, "w");
        ASSERT_ALWAYS(f);
        for(j = nb-1, k = 0; j >= 0; j -= 2) {
            // column number k in the output matrix corresponds to ideal
            // number tmp[j] (whose meaning can be fetched from the sos
            // file produced by purge)
            uint32_t z = tmp[j];
            int r = fwrite32_little(&z, 1, f);
            ASSERT_ALWAYS(r == 1);
        }
        fclose(f);
    }

#endif
    free(tmp);
}

// A line is "i i1 ... ik".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
static void
doAllAdds(int **newrows, char *str)
{
    char *t = str;
    int i, ii, destroy = 1;

    if(*t == '-'){
	destroy = 0;
	t++;
    }
    for(i = 0; (*t != ' ') && (*t != '\n'); t++)
	i = 10 * i + (*t - '0');
    if(!destroy)
	i--; // what a trick, man!
#if DEBUG >= 1
    fprintf(stderr, "first i is %d\n", i);
#endif
    if(*t != '\n'){
	++t;
	ii = 0;
	while(1){
	    if((*t == '\n') || (*t == ' ')){
#if DEBUG >= 1
		fprintf(stderr, "next ii is %d\n", ii);
#endif
		addRows(newrows, ii, i, -1);
		ii = 0;
	    }
	    else
		ii = 10 * ii + (*t - '0');
	    if(*t == '\n')
		break;
	    t++;
	}
    }
    if(destroy){
	// destroy initial row!
	free(newrows[i]);
	newrows[i] = NULL;
    }
}

#define STRLENMAX 2048

static int
oneFile(const char *sparsename, const char * sosname, int **sparsemat, int *colweight, purgedfile_stream_ptr ps, int **oldrows, int ncols, int small_nrows, int verbose, int skip, int bin)
{
    unsigned long W;
    int small_ncols;

    makeSparse(sparsemat, colweight, ps, 0, ncols, oldrows, verbose,0);

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber(sosname, &small_ncols, colweight, ncols);

    fprintf(stderr, "small_nrows=%d small_ncols=%d\n",small_nrows,small_ncols);

    double tt = seconds();
    fprintf(stderr, "Writing sparse representation to file\n");
    W = flushSparse(sparsename, sparsemat, small_nrows, small_ncols, colweight, skip, bin);
    fprintf(stderr, "#T# writing sparse: %2.2lf\n", seconds()-tt);
    fprintf(stderr, "# Weight(M_small) = %lu\n", W);

    return small_ncols;
}

static int
manyFiles(const char *sparsename, int **sparsemat, int *colweight, purgedfile_stream_ptr ps, int **oldrows, int ncols, int small_nrows, int verbose, int skip, int bin, int nslices)
{
    char *name;
    int small_ncols = 1; // always the +1 trick
    int slice, i, jmin, jmax, jstep, j, *tabnc;
    unsigned long *Wslice;

    // to add ".00" or to add ".infos"
    size_t namelen = strlen(sparsename)+7;
    name = (char *)malloc(namelen * sizeof(char));
    jmin = jmax = 0;
    // tabnc[slice] = number of columns in slice
    tabnc = (int *)malloc(nslices * sizeof(int));
    Wslice = (unsigned long *)malloc(nslices * sizeof(unsigned long));
    // we want nslices always
    if((ncols % nslices) == 0)
	jstep = ncols/nslices;
    else
	jstep = ncols/(nslices-1);
    ASSERT_ALWAYS(skip <= jstep);
    for(slice = 0; slice < nslices; slice++){
	// we operate on [jmin..jmax[
	jmin = jmax;
	jmax += jstep;
	if(slice == (nslices-1))
	    jmax = ncols;
	snprintf(name, namelen, "%s.%02d", sparsename, slice);
	fprintf(stderr, "Dealing with M[%d..%d[ -> %s\n", jmin, jmax, name);
        purgedfile_stream_rewind(ps);
	makeSparse(sparsemat,colweight,ps,jmin,jmax,oldrows,verbose,
		   nslices);
	// colweight[jmin..jmax[ contains the new weights of the
	// corresponding columns
	tabnc[slice] = 0;
	for(j = jmin; j < jmax; j++)
	    if(colweight[j] > 0){
		colweight[j] = small_ncols++;
		tabnc[slice]++;
	    }
	// warning: take care to the +1 trick
	Wslice[slice] = flushSparse(name, sparsemat, small_nrows, small_ncols,
				    colweight, bin, jmin == 0 ? skip : 0);
	// do not forget to clean the rows in sparsemat to be ready for
	// next time if any
	if(slice < nslices-1){
	    for(i = 0; i < small_nrows; i++)
		if(sparsemat[i] != NULL){
		    free(sparsemat[i]);
		    sparsemat[i] = NULL;
		}
	}
    }
    small_ncols--; // undoing the +1 trick
    // infos4Manu(name,sparsename,small_nrows,small_ncols,nslices,tabnc,Wslice);
    free(name);
    free(tabnc);
    free(Wslice);
    return small_ncols;
}

static void
build_newrows_from_file(int **newrows, int nrows, FILE *hisfile, uint64_t bwcostmin)
{
    uint64_t bwcost;
    unsigned long addread = 0;
    int i;
    char str[STRLENMAX];

    for(i = 0; i < nrows; i++){
	newrows[i] = (int *)malloc(2 * sizeof(int));
	newrows[i][0] = 1;
	newrows[i][1] = i;
    }
    fprintf(stderr, "Reading row additions\n");
	double tt = wct_seconds();
    while(fgets(str, STRLENMAX, hisfile)){
	addread++;
	if((addread % 100000) == 0)
	    fprintf(stderr, "%lu lines read at %2.2lf\n", addread, wct_seconds() - tt);
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Gasp: not a complete a line!");
	    fprintf(stderr, " I stop reading and go to the next phase\n");
	    break;
	}
	if(strncmp(str, "BWCOST", 6) != 0)
	    doAllAdds(newrows, str);
	else{
	    if(strncmp(str, "BWCOSTMIN", 9) != 0){
		sscanf(str+8, "%"PRIu64"", &bwcost);
		//	fprintf(stderr, "Read bwcost=%"PRIu64"\n", bwcost);
		if((bwcostmin != 0) && (bwcost == bwcostmin)){
		    // what a damn tricky feature!!!!!!!
		    fprintf(stderr, "Activating tricky stopping feature");
		    fprintf(stderr, " since I reached %"PRIu64"\n", bwcostmin);
		    break;
		}
	    }
	}
    }
}

static void
read_newrows_from_file(int **newrows, int nrows, FILE *file)
{
    int small_nrows, small_ncols, i, nc, k, *tmp;
    int rc;

    rc = fscanf(file, "%d %d", &small_nrows, &small_ncols);
    ASSERT_ALWAYS(rc == 2);

    for(i = 0; i < small_nrows; i++){
	rc = fscanf(file, "%d", &nc);
        ASSERT_ALWAYS(rc == 1);

	tmp = (int *)malloc((1+nc) * sizeof(int));
	tmp[0] = nc;
	for(k = 0; k < nc; k++) {
	    rc = fscanf(file, PURGE_INT_FORMAT, tmp+k+1);
            ASSERT_ALWAYS(rc == 1);
        }

	newrows[i] = tmp;
    }
    for( ; i < nrows; i++)
	newrows[i] = NULL;
}

// We start from M_purged which is nrows x ncols;
// we build M_small which is small_nrows x small_ncols.
// newrows[i] if != NULL, contains a list of the indices of the rows in
// M_purged that were added together to form this new row in M_small.
// TODO: replace this index by the index to rels directly to skip one
// indirection???
int
main(int argc, char *argv[])
{
    FILE *hisfile, *fromfile;
    uint64_t bwcostmin = 0;
    int nrows, ncols, nslices = 0;
    int **newrows, i, j, nb, *nbrels, **oldrows, *colweight;
    int ind, small_nrows, small_ncols, **sparsemat;
    int verbose = 0, writeindex;
    char str[STRLENMAX];
    int bin=0;
    char * rp;
    int skip=0;

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    // printing the arguments as everybody does these days
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list pl;

    param_list_init(pl);
    argv++,argc--;
    param_list_configure_knob(pl, "--verbose", &verbose);
    param_list_configure_knob(pl, "--binary", &bin);
    param_list_configure_alias(pl, "--verbose", "-v");

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        abort();
    }
    const char * purgedname = param_list_lookup_string(pl, "purged");
    const char * hisname = param_list_lookup_string(pl, "his");
    const char * sparsename = param_list_lookup_string(pl, "out");
    const char * indexname = param_list_lookup_string(pl, "index");
    const char * fromname = param_list_lookup_string(pl, "from");
    const char * sosname = param_list_lookup_string(pl, "sos");
    param_list_parse_int(pl, "binary", &bin);
    param_list_parse_int(pl, "nslices", &nslices);
    param_list_parse_int(pl, "skip", &skip);
    param_list_parse_uint64(pl, "bwcostmin", &bwcostmin);
    if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
        bin=1;

    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, purgedname);

    hisfile = fopen(hisname, "r");
    ASSERT(hisfile != NULL);
    rp = fgets(str, STRLENMAX, hisfile);
    ASSERT_ALWAYS(rp);

    // read parameters that should be the same as in purgedfile!
    sscanf(str, "%d %d", &nrows, &ncols);
    ASSERT_ALWAYS(nrows == ps->nrows);
    ASSERT_ALWAYS(ncols == ps->ncols);

    fprintf(stderr, "Original matrix has size %d x %d\n", nrows, ncols);
    newrows = (int **)malloc(nrows * sizeof(int *));

    if(fromname == NULL){
	writeindex = 1;
	build_newrows_from_file(newrows, nrows, hisfile, bwcostmin);
    } else {
	writeindex = 0;
	fromfile = fopen(fromname, "r");
	ASSERT(fromfile != NULL);
	read_newrows_from_file(newrows, nrows, fromfile);
    }
    fclose(hisfile);

    nbrels = (int *) malloc(nrows * sizeof(int));
    memset (nbrels, 0, nrows * sizeof(int));
    // nbrels[oldi] contains the number of new relations in which old
    // relation of index oldi takes part
    small_nrows = 0;
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    small_nrows++;
#if DEBUG >= 1
	    fprintf(stderr, "New row %d:", small_nrows-1);
#endif
	    for(j = 1; j <= newrows[i][0]; j++){
#if DEBUG >= 1
		fprintf(stderr, " %d", newrows[i][j]);
#endif
		nbrels[newrows[i][j]] += 1;
	    }
#if DEBUG >= 1
	    fprintf(stderr, "\n");
#endif
	}
    fprintf(stderr, "Allocating oldrows\n");
    oldrows = (int **)malloc(nrows * sizeof(int *));
    for(i = 0; i < nrows; i++){
	oldrows[i] = (int *)malloc((nbrels[i]+1) * sizeof(int));
	oldrows[i][0] = 0;
    }
    free (nbrels);
    fprintf(stderr, "Filling oldrows\n");
    for(i = 0, nb = 0; i < nrows; i++)
        if(newrows[i] != NULL){
	    // this is row of index nb in the new matrix
	    for(j = 1; j <= newrows[i][0]; j++){
		ind = newrows[i][j];
		oldrows[ind][0]++;
		oldrows[ind][oldrows[ind][0]] = nb;
	    }
	    nb++;
	}
#if DEBUG >= 1
    printOldRows(oldrows, nrows);
#endif
    colweight = (int *)malloc(ncols * sizeof(int *));
    memset(colweight, 0, ncols * sizeof(int *));
    fprintf(stderr, "Building sparse representation\n");
    sparsemat = (int **)malloc(small_nrows * sizeof(int *));
    for(i = 0; i < small_nrows; i++)
	sparsemat[i] = NULL;
    if(nslices == 0) {
	small_ncols = oneFile(sparsename, sosname, sparsemat, colweight, ps,
			      oldrows, ncols, small_nrows, verbose, skip, bin);
    } else {
	small_ncols = manyFiles(sparsename, sparsemat, colweight, ps,
				oldrows, ncols, small_nrows, verbose, skip, bin,
				nslices);
	exit(0);
    }
    if(writeindex){
	// this part depends on newrows only
	double tt = wct_seconds();
	fprintf(stderr, "Writing index file\n");
	makeIndexFile(indexname, nrows, newrows, small_nrows, small_ncols);
	fprintf(stderr, "#T# writing index file: %2.2lf\n", wct_seconds()-tt);
    }

    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);

    for(i = 0; i < small_nrows; i++)
	if(sparsemat[i] != NULL)
	    free(sparsemat[i]);
    free(sparsemat);

    for(i = 0; i < nrows; i++){
	if(newrows[i] != NULL)
	    free(newrows[i]);
	if(oldrows[i] != NULL)
	    free(oldrows[i]);
    }
    param_list_clear(pl);
    free(newrows);
    free(oldrows);
    free(colweight);
    return 0;
}
