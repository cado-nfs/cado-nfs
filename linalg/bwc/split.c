#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include "bwc_config.h"
#include "portability.h"
#include "macros.h"
#include "filenames.h"
#include "bw-common.h"
#include "params.h"
// #include "balancing.h"
#include "misc.h"


/* Maximum number of possible splits */
#define MAXSPLITS 16

/* splits for the different sites */
int splits[MAXSPLITS + 1];

int split_y = 0;
int split_f = 0;
int force = 0;

void usage()
{
    fprintf(stderr, "Usage: ./splits <options> [--split-y|--split-f] splits=0,<n1>,<n2>,...\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg n\n");
    fprintf(stderr, "Note: data files must be found in wdir !\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    FILE ** files;
    // balancing bal;

    param_list pl;
    param_list_init(pl);
    param_list_configure_switch(pl, "--split-y", &split_y);
    param_list_configure_switch(pl, "--split-f", &split_f);
    param_list_configure_switch(pl, "--force", &force);
    bw_common_init(bw, pl, &argc, &argv);
    int nsplits;
    nsplits = param_list_parse_int_list(pl, "splits", splits, MAXSPLITS, ",");

    /* one byte of data contains 8 elements of the field we are
     * considering
     */
    int scale[2] = {1, 8};      /* default for binary */

    param_list_parse_int_and_int(pl, "binary-ratio", scale, "/");

    /*
    const char * balancing_filename = param_list_lookup_string(pl, "balancing");
    if (!balancing_filename) {
        fprintf(stderr, "Required argument `balancing' is missing\n");
        usage();
    }
    balancing_read_header(bal, balancing_filename);
    */

    char * ifile = NULL;
    char * ofile_fmt = NULL;

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "ifile")) != NULL)
        ifile = strdup(tmp);
    if ((tmp = param_list_lookup_string(pl, "ofile-fmt")) != NULL)
        ofile_fmt = strdup(tmp);

    ASSERT_ALWAYS(!(!ifile ^ !ofile_fmt));

    if (param_list_warn_unused(pl)) usage();
    if (split_f + split_y + (ifile || ofile_fmt) != 1) {
        fprintf(stderr, "Please select one of --split-y or --split-f\n");
        usage();
    }
    if (nsplits <= 0) {
        fprintf(stderr, "Please indicate the splitting points\n");
        usage();
    }

    files = malloc(nsplits * sizeof(FILE *));

    for(int i = 0 ; i < nsplits ; i++) {
        ASSERT_ALWAYS((i == 0) == (splits[i] == 0));
        ASSERT_ALWAYS((i == 0) || (splits[i-1] < splits[i]));
    }
    if (splits[nsplits-1] != bw->n) {
        fprintf(stderr, "last split does not coincide with configured n\n");
        exit(1);
    }

    for(int i = 0 ; i < nsplits ; i++) {
        ASSERT_ALWAYS(splits[i] % scale[1] == 0);
    }

    nsplits--;

    int rc;

    /* prepare the file names */
    if (split_y) {
        rc = asprintf(&ifile, "%s.%u", Y_FILE_BASE, 0);
        rc = asprintf(&ofile_fmt, "%s.%u", V_FILE_BASE_PATTERN, 0);
    } else if (split_f) {
        ifile = strdup(LINGEN_F_FILE);
        ofile_fmt = strdup(F_FILE_SLICE_PATTERN);
    } else {
        ASSERT_ALWAYS(ifile);
        ASSERT_ALWAYS(ofile_fmt);
    }

    struct stat sbuf[1];
    if (stat(ifile, sbuf) < 0) {
        perror(ifile);
        exit(1);
    }

    if ((sbuf->st_size) % (splits[nsplits]*scale[0]/scale[1]) != 0) {
        fprintf(stderr, 
                "Size of %s (%ld bytes) is not a multiple of %d bytes\n",
                ifile, (long) sbuf->st_size,
                splits[nsplits]*scale[0]/scale[1]);
    }

    FILE * f = fopen(ifile, "rb");
    if (f == NULL) {
        fprintf(stderr,"%s: %s\n", ifile, strerror(errno));
        exit(1);
    }

#ifndef HAVE_MINGW
    /* mingw does not have link() */
    if (nsplits == 1) {
        char * fname;
        int i = 0;
        /* Special case ; a hard link is enough */
        rc = asprintf(&fname, ofile_fmt, splits[i], splits[i+1]);
        rc = stat(fname, sbuf);
        if (rc == 0 && !force) {
            fprintf(stderr,"%s already exists\n", fname);
            exit(1);
        }
        if (rc == 0 && force) { unlink(fname); }
        if (rc < 0 && errno != ENOENT) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        rc = link(ifile, fname);
        if (rc < 0) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        free(fname);
        free(ifile);
        free(ofile_fmt);
        return 0;
    }
#endif

    for(int i = 0 ; i < nsplits ; i++) {
        char * fname;
        rc = asprintf(&fname, ofile_fmt, splits[i], splits[i+1]);
        rc = stat(fname, sbuf);
        if (rc == 0 && !force) {
            fprintf(stderr,"%s already exists\n", fname);
            exit(1);
        }
        if (rc == 0 && force) { unlink(fname); }
        if (rc < 0 && errno != ENOENT) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        files[i] = fopen(fname, "wb");
        if (files[i] == NULL) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        free(fname);
    }

    void * ptr = malloc(splits[nsplits]*scale[0]/scale[1]);

    for(;;) {
        rc = fread(ptr, 1, splits[nsplits]*scale[0]/scale[1], f);
        if (rc != splits[nsplits]*scale[0]/scale[1] && rc != 0) {
            fprintf(stderr, "Unexpected short read\n");
            exit(1);
        }
        if (rc == 0)
            break;

        char * q = ptr;

        for(int i = 0 ; i < nsplits ; i++) {
            int d = splits[i+1]*scale[0]/scale[1] - splits[i]*scale[0]/scale[1];
            rc = fwrite(q, 1, d, files[i]);
            if (rc != d) {
                fprintf(stderr, "short write\n");
                exit(1);
            }
            q += d;
        }
    }

    param_list_clear(pl);
    for(int i = 0 ; i < nsplits ; i++) {
        fclose(files[i]);
    }
    fclose(f);

    free(ifile);
    free(ofile_fmt);
}

