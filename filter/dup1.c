/* dup1: 1st duplicate pass, split relation files into NSLICES slices
   (adapted from check).

   Usage:
   dup1 [-bz] -out <dir> file1 ... filen

   Files file1 ... filen are split into NSLICES slices in
   <dir>/0/filej ... <dir>/31/filej.

   If option -bz is given, then the output is compressed with bzip2
   instead of gzip.
   Input can be in gzipped or bzipped format.
*/

#include "cado.h"

#define NSLICES_LOG 2
#define NSLICES (1 << NSLICES_LOG)

#define CA 314159265358979323UL
#define CB 271828182845904523UL

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <limits.h> /* for CHAR_BIT */
#include <unistd.h>
#include <assert.h>

#include "macros.h"
#include "utils.h"
#include "relation.h"

/* Only (a,b) are parsed on input. This flags control whether we copy the
 * rest of the relation data to the output file, or if we content
 * ourselves with smaller .ab files */
static int only_ab = 0;

/* output relations in dirname/0/name, ..., dirname/31/name */
int
split_relfile (relation_stream_ptr rs, const char *name, const char *dirname, const char * outfmt, int *do_slice)
{
    FILE * f_in;
    int p_in;
    const char * suffix_in;

    int ok;
    char * oname[NSLICES];
    FILE * ofile[NSLICES];
    int p_out[NSLICES];

    unsigned long count[NSLICES] = {0,};
    uint64_t h;

    f_in = fopen_compressed_r(name, &p_in, &suffix_in);
    ASSERT_ALWAYS(f_in != NULL);

    char * newname = strdup(name);
    const char * suffix_out = outfmt;
    if (!suffix_out) suffix_out = suffix_in;
    // remove recognized suffix.
    ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
    newname[strlen(newname)-strlen(suffix_in)]='\0';
    for(int i = 0 ; i < NSLICES ; i++) {
        int rc = asprintf(&oname[i],
                only_ab ? "%s/%d/%s.ab%s" : "%s/%d/%s%s",
                dirname, i, path_basename(newname), suffix_out);
        ASSERT_ALWAYS(rc >= 0);
        ofile[i] = fopen_compressed_w(oname[i], &p_out[i], NULL);
        ASSERT_ALWAYS(ofile[i] != NULL);
    }
    free(newname);

    relation_stream_bind(rs, f_in);
    rs->parse_only_ab = 1;

    size_t pos0 = rs->pos;
    for (int lnum = 0 ; ; lnum++) {
        char line[RELATION_MAX_BYTES];
        size_t rpos = rs->pos - pos0;

        if (relation_stream_get(rs, line) < 0)
            break;

	ok = 1;

        h = CA * (uint64_t) rs->rel.a + CB * rs->rel.b;
        /* Using the low bit of h is not a good idea, since then
           odd values of i are twice more likely. The second low bit
           also gives a small bias with RSA768 (but not for random
           coprime a, b). We use here the NSLICES_LOG high bits.
        */
        int i = h >> (64 - NSLICES_LOG);

        /* print relation */
        if (do_slice[i]) {
            if (only_ab) {
                struct {
                    int64_t a;
                    uint64_t b;
                    size_t pos;
                } foo = { .a = rs->rel.a, .b=rs->rel.b, .pos=rpos, };
                fwrite(&foo, sizeof(foo), 1, ofile[i]);
            } else {
                if (fputs(line, ofile[i]) < 0) {
                    fprintf (stderr, "Problem writing to file %d\n", i);
                    exit(1);
                }
            }
        }
        count[i] ++;

        if (relation_stream_disp_progress_now_p(rs)) {
            fprintf (stderr,
                    "# split %d relations in %.1fs"
                    " -- %.1f MB/s -- %.1f rels/s\n",
                    rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
        }
    }

    relation_stream_unbind(rs);

    for (int i = 0; i < NSLICES; i++) {
      if (do_slice[i])
        {
          if (p_out[i]) pclose (ofile[i]); else fclose(ofile[i]);
          fprintf (stderr, "%d:%lu ", i, count[i]);
        }
      free(oname[i]);
    }
    fprintf (stderr, " %s\n", name);



    if (p_in) pclose(f_in); else fclose(f_in);
    return 0;
}

void usage()
{
    fprintf(stderr, "Usage: ./dup1 [-outfmt <fmt> | -bz] [-only <i>] -out <output_dir> [files...]\n");
}

int
main (int argc, char * argv[])
{
    int had_error = 0;

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int bz = 0;
    param_list_configure_knob(pl, "bz", &bz);
    param_list_configure_knob(pl, "ab", &only_ab);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort(); 
    }

    const char * dirname = param_list_lookup_string(pl, "out");
    int only_slice = -1;
    param_list_parse_int(pl, "only", &only_slice);
    const char * outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");

    if (param_list_warn_unused(pl)) {
        exit(1);
    }

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (!dirname)
        usage();
    if (bz) {
        if (outfmt) {
            fprintf(stderr, "-bz and -outfmt are exclusive");
            usage();
        } else {
            outfmt = ".bz2";
        }
    }
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "output compression format unsupported\n");
        usage();
    }

    int do_slice[NSLICES];

    if (only_slice == -1) { /* split all slices */
        for (int i = 0; i < NSLICES; i++)
          do_slice[i] = 1;
    } else { /* split only slide i */
        for (int i = 0; i < NSLICES; i++)
          do_slice[i] = i == only_slice;
    }

    if ((filelist != NULL) + (argc != 0) != 1) {
        fprintf(stderr, "Provide either -filelist or freeform file names\n");
        usage();
    }

    char ** files = filelist ? filelist_from_file(basepath, filelist) : argv;

    relation_stream rs;
    relation_stream_init(rs);
    for (char ** fp = files ; *fp ; fp++) {
        had_error |= split_relfile (rs, *fp, dirname, outfmt, do_slice);
    }
    relation_stream_trigger_disp_progress(rs);
    fprintf (stderr,
            "# split %d relations in %.1fs"
            " -- %.1f MB/s -- %.1f rels/s\n",
            rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    relation_stream_clear(rs);

    if (filelist) filelist_clear(files);

    param_list_clear(pl);

    if (had_error)
      return 1;
    else
      return 0;
}
