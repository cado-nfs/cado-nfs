#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include "macros.h"

static void
printVec (int *vec, int ncols)
{
    int i;

    fprintf(stderr, "vec=");
    for(i = 0; i < ncols; i++)
	fprintf(stderr, "%d", vec[i] & 1);
    fprintf(stderr, "\n");
}

static int
checkSparse (FILE *matfile, FILE *kerfile, int ncols, int nlimbs, int *vec,
             int compact, int verbose)
{
    long a;
    unsigned long w, b;
    int i, j, ret = 0, nc, cc, jj;
    char c;

    memset(vec, 0, ncols * sizeof(int));
    // get next dependance relation
    rewind(matfile);
    ret = fscanf(matfile, "%d %d", &i, &j);
    ASSERT_ALWAYS(ret == 2);
    while((c = getc(matfile)) != '\n');
    for(i = 0; i < nlimbs; ++i){
	ret = fscanf(kerfile, "%lx", &w);
	if(ret == -1)
	    break;
	ASSERT_ALWAYS (ret == 1);
	if(verbose)
	    fprintf(stderr, "w=%lx\n", w);
	for(j = 0; j < GMP_NUMB_BITS; ++j){
	    if(w & 1UL){
		if(verbose)
		    fprintf(stderr, "+R_%d\n", (i * GMP_NUMB_BITS)+j);
		if(compact == 0) {
		    ret = fscanf(matfile, "%ld %lu %d", &a, &b, &nc);
                    ASSERT_ALWAYS(ret == 3);
                } else {
		    ret = fscanf(matfile, "%d", &nc);
                    ASSERT_ALWAYS(ret == 1);
                }
		for(jj = 0; jj < nc; jj++){
		    ret = fscanf(matfile, "%d", &cc);
                    ASSERT_ALWAYS(ret == 1);
		    if(verbose >= 2)
			fprintf(stderr, "vec[%d]++\n", cc);
		    if(cc >= ncols){
			fprintf(stderr, "GASP: cc=%d > ncols=%d", cc, ncols);
			fprintf(stderr, " at line %d\n", (i*GMP_NUMB_BITS)+j);
		    }
		    vec[cc]++;
		}
		if(verbose >= 2)
		    printVec(vec, ncols);
		while((c = getc(matfile)) != '\n')
		    if(feof(matfile))
			break;
	    }
	    else{
		while((c = getc(matfile)) != '\n')
		    if(feof(matfile))
			break;
	    }
	    w >>= 1;
	}
    }
    return ret;
}

static void
checkVector (int *vec, int ncols, int skip)
{
    int ok, i;

    ok = 1;
    for(i = skip; i < ncols; i++)
	if(vec[i] & 1){
	    ok = 0;
	    break;
	}
    if(ok)
	fprintf(stderr, "y");
    else
	fprintf(stderr, "[n (%d)]", i);
}

static void
checkSparseAll (char *matname, char *kername, int compact, int skip, int verbose)
{
    FILE *matfile, *kerfile;
    int nrows, ncols, nlimbs, *vec, ndep, ret;

    matfile = fopen(matname, "r");
    kerfile = fopen(kername, "r");
    ret = fscanf(matfile, "%d %d", &nrows, &ncols);
    ASSERT_ALWAYS(ret == 2);
    nlimbs = (nrows / GMP_NUMB_BITS) + 1;
    vec = (int *)malloc(ncols * sizeof(int));
    ndep = 0;
    while(1){
	if(verbose)
	    fprintf(stderr, "Checking dep %d\n", ndep++);
        ret = checkSparse(matfile,kerfile,ncols,nlimbs,vec,compact,verbose);
	if(ret == -1)
	    break;
	if(verbose >= 2)
	    printVec(vec, ncols);
	checkVector(vec, ncols, skip);
    }
    fprintf(stderr, "\n");
    free(vec);
    fclose(matfile);
    fclose(kerfile);
}

static int
checkWithIndex (FILE *purgedfile, FILE *indexfile, FILE *kerfile, int verbose,
                int nrows, int ncols, int small_nrows, int *vec)
{
    unsigned long w;
    int i, j, k, ret, nlimbs, ind, nrel, r, nc;
    char *small_row_used, *rel_used;

    nlimbs = (small_nrows / GMP_NUMB_BITS) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    memset(small_row_used, 0, small_nrows * sizeof(char));
    for(i = 0; i < nlimbs; ++i){
        ret = fscanf(kerfile, "%lx", &w);
        if(ret == -1){
            free(small_row_used);
            return ret;
        }
        ASSERT_ALWAYS (ret == 1);
        if(verbose)
            fprintf(stderr, "w=%lx\n", w);
        for(j = 0; j < GMP_NUMB_BITS; ++j){
            if(w & 1UL){
                ind = (i * GMP_NUMB_BITS)+j;
                if(verbose)
                    fprintf(stderr, "+R_%d\n", ind);
                small_row_used[ind] = 1;
            }
            w >>= 1;
        }
    }
    // now map to the rels of the big matrix
    rel_used = (char *)malloc(nrows * sizeof(char));
    memset(rel_used, 0, nrows * sizeof(char));
    rewind(indexfile);
    ret = fscanf(indexfile, "%d %d", &i, &j); // skip first line
    ASSERT_ALWAYS(ret == 2);
    for(i = 0; i < small_nrows; i++){
        ret = fscanf(indexfile, "%d", &nrel);
        ASSERT_ALWAYS(ret == 1);
        for(j = 0; j < nrel; j++){
            ret = fscanf(indexfile, "%d", &r);
            ASSERT_ALWAYS(ret == 2);
            if(small_row_used[i])
                rel_used[r] ^= 1;
        }
    }
    // now really read the big matrix in
    rewind(purgedfile);
    ret = fscanf(purgedfile, "%d %d", &i, &j);
    ASSERT_ALWAYS(ret == 2);
    memset(vec, 0, ncols * sizeof(int));
    for(i = 0; i < nrows; i++){
        ret = fscanf(purgedfile, "%d %d", &j, &nc);
        ASSERT_ALWAYS(ret == 2);
        for(j = 0; j < nc; j++){
            ret = fscanf(purgedfile, "%d", &k);
            ASSERT_ALWAYS(ret == 1);
            if(rel_used[i])
		vec[k]++;
	}
    }
    return 1;
}

static void
checkWithIndexAll (char *purgedname, char *indexname, char *kername,
		   int skip, int verbose)
{
    FILE *purgedfile = fopen(purgedname, "r");
    FILE *indexfile = fopen(indexname, "r");
    FILE *kerfile = fopen(kername, "r");
    int nrows, ncols, *vec, small_nrows, small_ncols, ret, ndep;

    ret = fscanf(purgedfile, "%d %d", &nrows, &ncols);
    ASSERT_ALWAYS(ret == 2);
    ret = fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    ASSERT_ALWAYS(ret == 2);
    vec = (int *)malloc(ncols * sizeof(int));
    ndep = 0;
    while(1){
        if(verbose)
            fprintf(stderr, "Checking dep %d\n", ndep++);
        ret = checkWithIndex(purgedfile,indexfile,kerfile,verbose,nrows,ncols,small_nrows, vec);
	if(ret == -1)
	    break;
	checkVector(vec, ncols, skip);
    }
    fprintf(stderr, "\n");
    free(vec);
    fclose(kerfile);
    fclose(purgedfile);
    fclose(indexfile);
}

static void
usage (void)
{
  fprintf (stderr, "Usage: checkdep [-verbose] [-skip n] -purged file -index file -ker file\n");
  fprintf (stderr, "       checkdep [-verbose] [-compact] [-skip n] -mat file -ker file\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
    char *matname = NULL, *indexname = NULL, *purgedname = NULL;
    char *kername = NULL;
    int verbose = 0, compact = 0, skip = 0;

    while(argc > 1 && argv[1][0] == '-'){
        if(argc > 2 && strcmp (argv[1], "-mat") == 0){
	    matname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if(argc > 2 && strcmp (argv[1], "-ker") == 0){
	    kername = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if(argc > 2 && strcmp (argv[1], "-index") == 0){
	    indexname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if(argc > 2 && strcmp (argv[1], "-purged") == 0){
	    purgedname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if(argc > 2 && strcmp (argv[1], "-skip") == 0){
	    skip = atoi(argv[2]);
            argc -= 2;
            argv += 2;
        }
	else if(argc > 1 && strcmp (argv[1], "-compact") == 0){
	    compact = 1;
	    argc -= 1;
	    argv += 1;
	}
	else if(argc > 1 && strcmp (argv[1], "-verbose") == 0){
	    verbose = 1;
	    argc -= 1;
	    argv += 1;
	}
    }

    if (indexname != NULL)
      {
        if (purgedname == NULL || kername == NULL)
          usage ();
	checkWithIndexAll (purgedname, indexname, kername, skip, verbose);
      }
    else
      {
        if (matname == NULL || kername == NULL)
          usage ();
	checkSparseAll (matname, kername, compact, skip, verbose);
      }
    return 0;
}
