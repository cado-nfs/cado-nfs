/* This takes a range of row numbers in the post-merge relation matrix,
 * and checks the corresponding factorization into rational primes and
 * algebraic prime ideals.
 */

/* Usage string:
 *
  $bindir/filter/explain -sosp $wdir/$name.sos-purged.bin -sosr $wdir/$name.sos-replay.bin -purged $wdir/$name.purged -index $wdir/$name.index -matrix $wdir/$name.small.bin -skip 32 -start 0 -end 10

 * Do not forget the skip parameter -- otherwise the output is completely
 * useless.
 *
 * The file name.sos-replay.bin  is created by replay.
 * The file sos-purged.bin is created by purge in text format, the
 * following snippet converts it to binary (yeah, purge should do this).
#include <stdio.h>
main()
{
    unsigned int i,p,r;
    for( ; scanf("%u %x %x", &i,&p,&r) == 3 ; ) {
        fwrite(&p, sizeof(unsigned int), 1, stdout);
        fwrite(&r, sizeof(unsigned int), 1, stdout);
    }
}
 *
 * -start defaults to 0
 * -end defaults to nrows
 */

#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "utils.h"

/* copy my_malloc_int interface from purge {{{ */
#define BLOCK_SIZE 1000000 /* memory blocks are allocated of that # of int's */

/* mmi_list is a list of blocks, each one of BLOCK_SIZE ints */
int **mmi_list = NULL;
int mmi_current = -1; /* index of current block */
unsigned long mmi_used = BLOCK_SIZE; /* usage of current block */

/* converts an index to a pointer */
static int* mmi_translate(int v)
{
    int block = v / BLOCK_SIZE;
    return mmi_list[block] + (v % BLOCK_SIZE);
}

/* return the index of an array of n ints */
static int
mmi_get (unsigned long n)
{
  int index;

  if (mmi_used + n > BLOCK_SIZE) {
      /* new block */
      mmi_current ++;
      mmi_list = (int**) realloc (mmi_list, (mmi_current + 1) * sizeof (int*));
      mmi_list[mmi_current] = (int*) malloc (BLOCK_SIZE * sizeof (int));
      mmi_used = 0;
    }
  index = mmi_current * BLOCK_SIZE + mmi_used;
  mmi_used += n;
  return index;
}

static void
mmi_clear (void)
{
  while (mmi_current >= 0)
    free (mmi_list[mmi_current--]);
  free (mmi_list);
  mmi_list = NULL;
  mmi_used = BLOCK_SIZE;
}
/* }}} */

typedef int (*sortfunc_t) (const void*, const void*);

int uint32cmp(const uint32_t * a, const uint32_t * b)
{
    return (*a>*b)-(*b>*a);
}


struct ideal {
    uint32_t p;
    uint32_t r;
};
struct ideal * ideal_decode;
int intcmp2(const int * a, const int * b)
{
    if (ideal_decode[*a].r == -2u && ideal_decode[*b].r != -2u)
        return -1;
    if (ideal_decode[*a].r != -2u && ideal_decode[*b].r == -2u)
        return 1;
    int r = uint32cmp(&ideal_decode[*a].p, &ideal_decode[*b].p);
    if (r) return r;
    r = uint32cmp(&ideal_decode[*a].r, &ideal_decode[*b].r);
    return r;
}


int
main (int argc, char **argv)
{
    param_list pl;
    param_list_init(pl);

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        abort(); 
    }

    const char * sosp = param_list_lookup_string(pl, "sosp");
    const char * sosr = param_list_lookup_string(pl, "sosr");
    const char * purged = param_list_lookup_string(pl, "purged");
    const char * index = param_list_lookup_string(pl, "index");
    const char * matrix = param_list_lookup_string(pl, "matrix");
    const char * outname = param_list_lookup_string(pl, "out");

    int skip=0;
    param_list_parse_int(pl, "skip", &skip);

    ASSERT_ALWAYS(sosp && sosr && purged && index && matrix);

    /* just avoid warnings. */
    param_list_lookup_string(pl, "start");
    param_list_lookup_string(pl, "end");
    if (param_list_warn_unused(pl))
        exit(1);

    FILE * outfile = stdout;
    int pipe = 0;
    if (outname) {
        outfile = fopen_compressed_w(outname,&pipe,NULL);
    }

    /* Prepare the list of (a,b) pairs */
    unsigned int n_ab;
    unsigned int n_ideals;

    size_t nb_stored_ints_purged = 0;
    struct abpair {
        int64_t a;
        uint64_t b;
        int mmi_index;
    };
    struct abpair * abs;
    {
        purgedfile_stream ps;
        purgedfile_stream_init(ps);
        purgedfile_stream_openfile(ps, purged);
        n_ab = ps->nrows;
        n_ideals = ps->ncols;
        fprintf(stderr, "Allocating room for %u abs: %.2f MB\n",
                n_ab,
                (ps->nrows * sizeof(struct abpair)) /1048576.0);
        abs = malloc(ps->nrows * sizeof(struct abpair));
        ASSERT_ALWAYS(abs);
        for(int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
            abs[i].a = ps->a;
            abs[i].b = ps->b;
            abs[i].mmi_index = mmi_get(ps->nc + 1);
            int * ptr = mmi_translate(abs[i].mmi_index);
            *ptr++ = ps->nc;
            for(int k = 0 ; k < ps->nc ; k++) {
                ptr[k] = ps->cols[k];
            }
            nb_stored_ints_purged += (ps->nc + 1);
        }
        fprintf(stderr, "Storing (a,b) pairs factorizations: %.2f MB\n",
                nb_stored_ints_purged * sizeof(int) /1048576.0);
        purgedfile_stream_closefile(ps);
        purgedfile_stream_clear(ps);
    }


    /* preparing the correspondence between ideals appearing in .purged,
     * and real ideals
     */
    {
        fprintf(stderr, "Allocating room for %u ideal descriptons: %.2f MB\n",
                n_ideals,
                (n_ideals * sizeof(struct ideal)) /1048576.0);
        ideal_decode = malloc(n_ideals * sizeof(struct ideal));
        ASSERT_ALWAYS(ideal_decode);
        FILE * f = fopen(sosp, "r"); ASSERT_ALWAYS(f);
        int r = fread(ideal_decode, sizeof(struct ideal), n_ideals, f);
        ASSERT_ALWAYS(n_ideals == (unsigned int) r);
        fclose(f);
    }

    /* preparing to decode ideals which appear in the matrix */
    unsigned int n_ideals_inmat;
    uint32_t * ideals_inmat;
    {
        struct stat sbuf[1];
        if (stat(sosr, sbuf) < 0) { perror(sosr); exit(1); }
        n_ideals_inmat = sbuf->st_size / sizeof(uint32_t);
        fprintf(stderr, "Allocating room for %u ideal indirections: %.2f MB\n",
                n_ideals_inmat,
                (n_ideals_inmat * sizeof(uint32_t)) /1048576.0);
        ideals_inmat = malloc(n_ideals_inmat * sizeof(uint32_t));
        ASSERT_ALWAYS(ideals_inmat);
        FILE * f = fopen(sosr, "r"); ASSERT_ALWAYS(f);
        int r = fread(ideals_inmat, sizeof(uint32_t), n_ideals_inmat, f);
        ASSERT_ALWAYS(n_ideals_inmat == (unsigned int) r);
        fclose(f);
    }

    /* read relation-sets and matrix rows simultaneously */
    {
        FILE * f = fopen(index, "r"); ASSERT_ALWAYS(f);
        unsigned int nrows;
        unsigned int ncols;
        int r = fscanf(f, "%u %u",&nrows, &ncols);
        unsigned int i0 = 0;
        unsigned int i1 = nrows;
        param_list_parse_uint(pl, "start", &i0);
        param_list_parse_uint(pl, "end", &i1);
        ASSERT_ALWAYS(r == 2);
        ASSERT_ALWAYS(ncols == n_ideals_inmat);

        FILE * g = fopen(matrix, "r"); ASSERT_ALWAYS(f);

        fprintf(stderr, "Dumping relset info for relsets [%d..%d[ (complete range is [%d..%d[\n",i0,i1,0,nrows);
        for(unsigned int i = 0 ; i < i0 ; i++) {
            int nc;
            r = fscanf(f, "%u", &nc);
            ASSERT_ALWAYS(r == 1);
            for(int j = 0 ; j < nc ; j++) {
                unsigned int k;
                r = fscanf(f, "%x", &k);
                ASSERT_ALWAYS(r == 1);
            }
            uint32_t z;
            int r = fread(&z, sizeof(uint32_t), 1, g);
            ASSERT_ALWAYS(r == 1);
            for(int j = 0 ; j < nc ; j++) {
                r = fread(&z, sizeof(uint32_t), 1, g);
                ASSERT_ALWAYS(r == 1);
            }
        }
        int * fac = NULL;
        int fac_alloc = 0;
        int fac_size = 0;
        for(unsigned int i = i0 ; i < i1 ; i++) {
            /* Use the index file to decode the relation-set */
            int nc;
            fac_size = 0;
            r = fscanf(f, "%u", &nc);
            ASSERT_ALWAYS(r == 1);
            fprintf(outfile, "relset %u\n", i);
            for(int j = 0 ; j < nc ; j++) {
                unsigned int k;
                r = fscanf(f, "%x", &k);
                ASSERT_ALWAYS(r == 1);
                ASSERT_ALWAYS(k < n_ab);
                fprintf(outfile, "%"PRId64" %"PRIu64"\n", abs[k].a, abs[k].b);
                int * ff = mmi_translate(abs[k].mmi_index);
                if (fac_size + *ff > fac_alloc) {
                    fac_alloc += *ff + 64 + fac_alloc / 2;
                    fac = realloc(fac, fac_alloc * sizeof(int));
                }
                memcpy(fac + fac_size, ff + 1, *ff * sizeof(int));
                fac_size += *ff;
            }
            fprintf(outfile, "expected\n");
            qsort(fac, fac_size, sizeof(int), (sortfunc_t)&intcmp2);
            for(int s = 0; s < fac_size ;) {
                int idl = fac[s];
                int v = 0;
                for( ; s + v < fac_size && fac[s+v] == idl ; v++) ;
                ASSERT_ALWAYS((unsigned int) idl < n_ideals);
                fprintf(outfile, "%d %d %d\n", v, ideal_decode[idl].p, ideal_decode[idl].r);
                s+=v;
            }
            fprintf(outfile, "recorded\n");
            uint32_t z;
            int r = fread(&z, sizeof(uint32_t), 1, g);
            ASSERT_ALWAYS(r == 1);
            nc = z;
            fac_size = 0;

            if (fac_size + nc > fac_alloc) {
                fac_alloc += nc + 64 + fac_alloc / 2;
                fac = realloc(fac, fac_alloc * sizeof(int));
            }

            for(int j = 0 ; j < nc ; j++) {
                r = fread(&z, sizeof(uint32_t), 1, g);
                z += skip;
                ASSERT_ALWAYS(r == 1);
                ASSERT_ALWAYS(z < n_ideals_inmat);
                int idl = ideals_inmat[z];
                fac[fac_size++]=idl;
            }
            qsort(fac, fac_size, sizeof(int), (sortfunc_t)&intcmp2);
            for(int s = 0; s < fac_size ;) {
                int idl = fac[s];
                int v = 0;
                for( ; s + v < fac_size && fac[s+v] == idl ; v++) ;
                ASSERT_ALWAYS((unsigned int) idl < n_ideals);
                fprintf(outfile, "%d %d %d\n", v, ideal_decode[idl].p, ideal_decode[idl].r);
                s+=v;
            }
            fprintf(outfile, "\n");
        }
        free(fac);
        fclose(f);
        fclose(g);
    }

    if (outname) {
        if (pipe)
            pclose(outfile);
        else
            fclose(outfile);
    }

    mmi_clear();
    param_list_clear(pl);
    return 0;
}
