/* dup2: 2nd pass

   Usage: dup2 [-out <dir>] [-rm] [-filelist <fl>] -K <K> file1 ... filen

   Puts non-duplicate files (among file1 ... filen only) into:
   <dir>/file1 ... <dir>/filen

   If -out <dir> is missing, no output is done.

   Input files can be given on command line, or via a filelist file.
   (it is possible to do both)

   If -rm is given, the input files are removed after having been treated.
   This is also the case if -out . is given. 

   By default, the output will be bzipped/gzipped according to the status
   of the input file.

   Allocates a hash-table of 2^k 32-bit entries.

   Algorithm: for each (a,b) pair, we compute h(a,b) = (CA*a+CB*b) % 2^64.

   h has 64 bits. We know bits 2..6 are identical for a given slice, thus
   we remove them, and we obtain a 59-bit value h'.

   Let h' = j * 2^k + i.

   We store j % 2^32 at the next empty cell after index i in the hash-table.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>     // for unlink
#include <inttypes.h>
#include <ctype.h>  // for isspace
#include "utils.h"
#include "relation.h"

#include "macros.h"

#define CA 271828182845904523UL
#define CB 577215664901532889UL

unsigned long dupl = 0;
unsigned long nodu = 0;
double cost = 0.0;

/* sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
   and check for hash collisions */
unsigned long sanity_size;
int64_t  *sanity_a;
uint64_t *sanity_b;
unsigned long sanity_checked = 0;
unsigned long sanity_collisions = 0;
/* end sanity check */

int rm = 0;

/* infile is the input file
   if dirname is NULL, no output is done */
unsigned long
remove_dup_in_files (char ** files, const char *dirname, const char * outfmt, uint32_t * H, unsigned long K)
{
    FILE * f_in;
    int p_in;
    const char * suffix_in;
    FILE * f_out;
    int p_out;
    const char * suffix_out;

    relation_stream rs;
    relation_stream_init(rs);

    rs->parse_only_ab = 1;

    for( ; *files ; files++) {
        const char * name = *files;

        f_in = fopen_compressed_r(name, &p_in, &suffix_in);
        ASSERT_ALWAYS(f_in != NULL);

        suffix_out = outfmt ? outfmt : suffix_in;
        char * newname = strdup(name);
        ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
        newname[strlen(newname)-strlen(suffix_in)]='\0';
        int only_ab = has_suffix(newname, ".ab");
        char * oname = NULL;
        char * oname_tmp = NULL;
        if (dirname) {
            int rc;
            rc = asprintf(&oname_tmp, "%s/pre-%s%s", dirname, path_basename(newname), suffix_out);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&oname, "%s/%s%s", dirname, path_basename(newname), suffix_out);
            ASSERT_ALWAYS(rc >= 0);

            f_out = fopen_compressed_w(oname_tmp, &p_out, NULL);
            ASSERT_ALWAYS(f_out != NULL);
        } else {
            f_out = NULL;
        }
        free(newname);

        relation_stream_bind(rs, f_in);

        uint64_t h;
        uint32_t i, j;
        static double factor = 1.0;
        double full_table;

        for (;;) {
            char line[RELATION_MAX_BYTES];
            struct {
                int64_t a;
                uint64_t b;
                size_t pos;
            } desc;

            if (only_ab) {
                int readbytes = fread(&desc, sizeof(desc), 1, f_in);
                if (readbytes == 0)
                    break;
                rs->pos += readbytes * sizeof(desc);
                rs->nrels++;
                rs->rel.a = desc.a;
                rs->rel.b = desc.b;
            } else {
                if (relation_stream_get(rs, line) < 0)
                    break;
            }

            h = CA * (uint64_t) rs->rel.a + CB * rs->rel.b;
            /* dup1 uses the high 5 bits of h to identify the slices 
               but now we use a different hash function, so we can keep these
               5 bits
               */
            i = h % K;
            j = (uint32_t) (h >> 32);
            /* Note: in the case where K > 2^32, i and j share some bits.
             * The high bits of i are in j. These bits correspond therefore to
             * far away positions in the tables, and keeping them in j can only
             * help.
             * FIXME:
             * TODO: that's wrong!!! it would be better do take i from high
             * bits instead!
             */
            while (H[i] != 0 && H[i] != j) {
                i++;
                if (UNLIKELY(i == K))
                    i = 0;
                cost++;
            }
            if (H[i] == j) {
                dupl++;
                continue;		/* probably duplicate */
            }

            if (i < sanity_size) {
                sanity_checked++;
                if (sanity_a[i] == 0) {
                    sanity_a[i] = rs->rel.a;
                    sanity_b[i] = rs->rel.b;
                } else if (sanity_a[i] != rs->rel.a) {
                    sanity_collisions++;
                    fprintf(stderr,
                            "Collision between (%" PRId64 ",%" PRIu64 ") and (%"
                            PRId64 ",%" PRIu64 ")\n", sanity_a[i], sanity_b[i], rs->rel.a,
                            rs->rel.b);
                }
            }

            nodu++;
            if (cost >= factor * (double) rs->nrels) {
                full_table = 100.0 * (double) nodu / (double) K;
                fprintf(stderr, "Warning, hash table is %1.0f%% full\n",
                        full_table);
                if (full_table >= 99) {
                    fprintf(stderr, "Error, hash table is full\n");
                    exit(1);
                }
                factor += 1.0;
            }
            /* now H[i] = 0 */
            H[i] = j;
            if (f_out) {
                if (only_ab) {
                    fwrite(&desc, sizeof(desc), 1, f_out);
                } else {
                    fputs(line, f_out);
                }
            }

            if (relation_stream_disp_progress_now_p(rs)) {
                fprintf(stderr,
                        "Read %d relations, %lu duplicates (%1.2f%%)"
                        " in %.1f s -- %.1f MB/s, %.1f rels/s\n",
                        rs->nrels, dupl, 100.0 * (double) dupl / (double) rs->nrels,
                        rs->dt, rs->mb_s, rs->rels_s);
            }
        }
        relation_stream_unbind(rs);
        if (p_out) pclose(f_out); else fclose(f_out);
        if (p_in) pclose(f_in); else fclose(f_in);


        if (dirname) {
            if (rm || (strcmp(dirname, ".") == 0)) {
                fprintf(stderr, "Removing old file %s\n", name);
                int ret = unlink(name);
                if (ret) {
                    perror("Problem removing file");
                    fprintf(stderr, "Let's hope that it's ok to continue!\n");
                }
            }
            fprintf(stderr, "%s/{%s => %s}\n",
                    dirname, path_basename(oname_tmp), path_basename(oname));
            int ret = rename(oname_tmp, oname);
            if (ret) {
                perror("Problem renaming result file");
                fprintf(stderr, "Let's hope that it's ok to continue!\n");
            }
            free(oname);
            free(oname_tmp);
        }
    }
    relation_stream_trigger_disp_progress(rs);
    fprintf(stderr,
            "Read %d relations, %lu duplicates (%1.2f%%)"
            " in %.1f s -- %.1f MB/s, %.1f rels/s\n",
            rs->nrels, dupl, 100.0 * (double) dupl / (double) rs->nrels,
            rs->dt, rs->mb_s, rs->rels_s);
    unsigned long rread = rs->nrels;
    relation_stream_clear(rs);
    return rread;
}

void usage()
{
    fprintf (stderr, "Usage: dup2 [-rm] [-out <dir>] [-filelist <fl>] -K <K> file1 ... filen\n");
    exit (1);
}

int main (int argc, char *argv[])
{
    unsigned long K = 0;
    uint32_t *H;

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int bz = 0;
    param_list_configure_knob(pl, "bz", &bz);
    param_list_configure_knob(pl, "rm", &rm);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort(); 
    }

    const char * dirname = param_list_lookup_string(pl, "out");
    const char * outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");

    param_list_parse_ulong(pl, "K", &K);

    if (param_list_warn_unused(pl)) {
        usage();
    }

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (K == 0)
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


  /* sanity check: since we allocate two 64-bit words for each, instead of
     one 32-bit word for the hash table, taking K/100 will use 2.5% extra
     memory */
  sanity_size = 1 + (K / 100);
  fprintf (stderr, "[checking true duplicates on sample of %lu cells]\n",
           sanity_size);
  sanity_a = (int64_t*)  malloc (sanity_size * sizeof (int64_t));
  memset (sanity_a, 0, sanity_size * sizeof (int64_t));
  sanity_b = (uint64_t*) malloc (sanity_size * sizeof (uint64_t));

  H = (uint32_t*) malloc (K * sizeof (uint32_t));
  if (H == NULL)
    {
      fprintf (stderr, "Error, cannot allocate hash table\n");
      exit (1);
    }
  memset (H, 0, K * sizeof (uint32_t));
  fprintf (stderr, "Allocated hash table of %lu entries (%luMb)\n", K,
           (K * sizeof (uint32_t)) >> 20);


  if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Provide either -filelist or freeform file names\n");
      usage();
  }

  char ** files = filelist ? filelist_from_file(basepath, filelist) : argv;
  unsigned long rread = remove_dup_in_files (files, dirname, outfmt, H, K);
  if (filelist) filelist_clear(files);


  fprintf (stderr, "Read %lu relations, %lu duplicates (%1.2f%%)\n",
           rread, dupl, 100.0 * (double) dupl / (double) rread);

  free (H);

  fprintf (stderr, "     %lu remaining relations (hash table %1.2f%% full)\n",
           nodu, 100.0 * (double) nodu / (double) K);
  fprintf (stderr, "     Hash-table cost %1.2f per relation\n",
           1.0 + cost / (double) rread);

  free (sanity_a);
  free (sanity_b);
  fprintf (stderr, "[found %lu true duplicates on sample of %lu relations]\n",
           sanity_collisions, sanity_checked);

  param_list_clear(pl);
  return 0;
}
