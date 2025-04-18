/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015, 2016
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

/* Program to update the thresholds table given a tuning result file. */

/* How to use this program:

   1) Use the option -o to provide the location of thresholds.h file
      which will be overwritten.
   2) feed as stdin the result of either the tunetoom or tunefft program.

*/

#define _DEFAULT_SOURCE /* _BSD_SOURCE is deprecated */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>             /* int64_t */
#include <inttypes.h>           /* SCNd64 */
#include <limits.h>		/* for LONG_MAX */
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/utsname.h>        /* for uname */

/* We need some constants from gf2x */
#include "gf2x/gf2x-impl.h"

#include "replace.h"

#ifndef MAX
#define MAX(a,b)        ((a)<(b) ? (b) : (a))
#endif

const char * thresholds_file = "thresholds.h";

struct hash_define replacements[100];
unsigned long nrepl=0;

/* We build our private table here */
short best_tab[GF2X_TOOM_TUNING_LIMIT];
short best_utab[GF2X_TOOM_TUNING_LIMIT];
unsigned long blim = 0;
unsigned long ulim = 0;

struct fft_tuning_pair {
    int64_t sz;
    int method;
};

struct fft_tuning_pair fft_tbl[1000];
unsigned long nstairs = 0;

char * toom_info_string = "\"placeholder\"";
char * fft_info_string = "\"placeholder\"";

#define MINI_GF2X_MUL_TOOM_THRESHOLD 	17
#define MINI_GF2X_MUL_TOOMW_THRESHOLD	8
#define MINI_GF2X_MUL_TOOM4_THRESHOLD	30
#define MINI_GF2X_MUL_TOOMU_THRESHOLD	33

#define BESTMIN (GF2X_MUL_KARA_THRESHOLD-1)
#define BESTMINU (GF2X_MUL_TOOMU_THRESHOLD-1)

int startswith(char ** pptr, const char * key)
{
    unsigned long l = strlen(key);
    if (strncmp(*pptr, key, l) != 0)
        return 0;
    *pptr += l;
    for( ; **pptr && isspace(**pptr) ; (*pptr)++) ;
    return 1;
}

void read_table()
{
    memset(best_tab, -1, sizeof(best_tab));
    memset(best_utab, -1, sizeof(best_utab));
    nstairs = 0;

    for(;;) {
        char buf[1024];
        int r;
        int64_t sz;
        int m;
        char token[16];
        char * ptr;

        ptr = fgets(buf, sizeof(buf), stdin);
        if (feof(stdin))
            break;
        if (buf[0] == '#')
            continue;
        ptr = buf;
        if (startswith(&ptr, "info")) {
            if (startswith(&ptr, "-fft")) {
                fft_info_string = strdup(ptr);
            } else if (startswith(&ptr, "-toom")) {
                toom_info_string = strdup(ptr);
            } else {
                abort();
            }
            continue;
        }
        r = sscanf(buf, "%15s %" SCNd64 " %d", token, &sz, &m);
        if (r != 3) {
            fprintf(stderr, "Bad input line: %s", buf);
            exit(1);
        }

        if (strcmp(token, "fft") == 0) {
            fft_tbl[nstairs].sz = sz;
            fft_tbl[nstairs].method = m;
            nstairs++;
            continue;
        }

        if (strcmp(token, "toom") == 0) {
            best_tab[sz - 1] = m;
            if (sz > (int64_t) blim) {
                blim = sz;
            }
            continue;
        }

        if (strcmp(token, "utoom") == 0) {
            best_utab[sz - 1] = m;
            if (sz > (int64_t) ulim) {
                ulim = sz;
            }
            continue;
        }
    }
}

typedef int (*sortfunc_t) (const void *, const void *);
int compare_fft_tuning(const struct fft_tuning_pair * a, const struct fft_tuning_pair * b)
{
    return a->sz - b->sz;
}

void fill_holes(short * tab, long lim, short dflt)
{
    long i,j;
    for(i = 0 ; i < lim ; ) {
        for( ; i < lim && tab[i] >= 0 ; i++);
        if (i == lim)
            break;
        for(j = i ; j < lim && tab[j] < 0 ; j++);
        short m0 = i ? tab[i-1] : 0;
        short m1 = (j < lim) ? tab[j] : dflt;
        long k;
        k = j-1;
        for( ; i <= k ; i++,k--) {
            tab[k]=m1;
            tab[i]=m0;
        }
    }
}

void prepare_and_push_hash_define_tbl(const char * name, short * tbl, size_t n)
{
    size_t sz = 3 + 3*n+(n+19)/20*4 + 10;
    size_t h = 0;
    char * table_string = malloc(sz);
    size_t i;

    h += snprintf(table_string + h, sz-h, "{");
    for (i = 1 ; i <= n ; i++) {
        if ((i-1) % 20 == 0) {
            h += snprintf(table_string + h, sz-h, "\t\\\n\t");
        }
        h += snprintf(table_string + h, sz-h, "%d, ", tbl[i-1]);
    }
    /* The last carriage return is added by the replacement engine anyway
     */
    h += snprintf(table_string + h, sz-h, "}");
    if (h >= sz) abort();

    set_hash_define(replacements + nrepl++, name, table_string);

    free(table_string);
}

void prepare_and_push_hash_define_fft_tbl(const char *name,
					  struct fft_tuning_pair *tbl,
					  size_t n)
{
    size_t sz = 10 + 30 * n;
    size_t h = 0;
    char *table_string = malloc(sz);
    size_t i;

    h += snprintf(table_string + h, sz - h, "{");
    for (i = 0; i < n; i++) {
	if (i % 4 == 0) {
	    h += snprintf(table_string + h, sz - h, "\t\\\n\t");
	}
        if (tbl[i].sz >> 32) {
            // make sure that the constant looks somewhat right (even
            // though it's a stupid idea, using a thresholds file for a
            // 64bit cpu on a 32bit box shouldn't break the meaning of
            // the data).
            h += snprintf(table_string + h, sz - h,
                    "{ INT64_C(%"PRId64"), %d }, ",
                    tbl[i].sz, tbl[i].method);
        } else {
            h += snprintf(table_string + h, sz - h,
                    "{ %"PRId64", %d }, ",
                    tbl[i].sz, tbl[i].method);
        }
    }
    h += snprintf(table_string + h, sz - h, "}\n");
    if (h >= sz)
	abort();

    set_hash_define(replacements + nrepl++, name, table_string);

    free(table_string);
}


void get_toom_thresholds()
{
    /* Get the values for the balanced thresholds:
     * GF2X_MUL_TOOM_THRESHOLD
     * GF2X_MUL_TOOMW_THRESHOLD
     * GF2X_MUL_TOOM4_THRESHOLD
     * GF2X_MUL_TOOM4_ALWAYS_THRESHOLD
     */
    unsigned long t3 = UINT_MAX;
    unsigned long tw = UINT_MAX;
    unsigned long t4 = UINT_MAX;
    unsigned long t4a = UINT_MAX;

    unsigned long i;

    for (i = BESTMIN + 1; i && i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC3) {
	    t3 = i;
	    break;
	}
    }
    if (t3 == UINT_MAX) { t3 = MAX(blim, 17); }

    for (i = BESTMIN + 1; i && i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC3W) {
	    tw = i;
	    break;
	}
    }
    if (tw == UINT_MAX) { tw = MAX(blim, 8); }

    for (i = BESTMIN + 1; i && i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC4) {
	    t4 = i;
	    break;
	}
    }
    if (t4 == UINT_MAX) { t4 = MAX(blim, 30); }

    if (blim && t4 < blim) {
	for (t4a = blim; t4a >= MINI_GF2X_MUL_TOOM4_THRESHOLD; t4a--) {
	    if (best_tab[t4a - 1] != GF2X_SELECT_TC4) {
		t4a++;
		break;
	    }
	}
    } else {
	t4a = blim;
    }

    /* Now do some sanity checks */
    if (!(tw <= t3)) {
	fprintf(stderr,
		"Note: we normally expect GF2X_MUL_TOOMW_THRESHOLD(%lu) to be below GF2X_MUL_TOOM_THRESHOLD(%lu). It is not the case here, but it is not an error either.\n",
		tw, t3);
    }

    /* GF2X_MUL_TOOM4_ALWAYS_THRESHOLD must be >= 30 */
    if (t4a < 30)
      {
	fprintf (stderr,
                 "GF2X_MUL_TOOM4_ALWAYS_THRESHOLD(%lu) must be >= 30\n",
                 t4a);
        exit (1);
      }

    set_hash_define_int(replacements + nrepl++, "GF2X_MUL_TOOM_THRESHOLD", t3);
    set_hash_define_int(replacements + nrepl++, "GF2X_MUL_TOOMW_THRESHOLD", tw);
    set_hash_define_int(replacements + nrepl++, "GF2X_MUL_TOOM4_THRESHOLD", t4);
    set_hash_define_int(replacements + nrepl++,
            "GF2X_MUL_TOOM4_ALWAYS_THRESHOLD", t4a);
    blim = t4a;
}

void get_utoom_thresholds()
{
    /**************************************************/
    /* Get the values for the unbalanced thresholds:
     * GF2X_MUL_TOOM3U_THRESHOLD
     * GF2X_MUL_TOOM3U_ALWAYS_THRESHOLD
     */
    unsigned long tu = UINT_MAX;
    unsigned long tua = UINT_MAX;
    unsigned long i;

    for (i = BESTMINU + 1; i && i <= ulim; i++) {
	if (best_utab[i - 1] == GF2X_SELECT_UNB_TC3U) {
	    tu = i;
	    break;
	}
    }
    if (tu == UINT_MAX) { tu = MAX(blim, 33); }

    if (ulim && tu < ulim) {
	for (tua = ulim; tua >= MINI_GF2X_MUL_TOOMU_THRESHOLD; tua--) {
	    if (best_utab[tua - 1] != GF2X_SELECT_UNB_TC3U) {
		tua++;
		break;
	    }
	}
    } else {
	tua = ulim;
    }
    set_hash_define_int(replacements + nrepl++, "GF2X_MUL_TOOMU_THRESHOLD", tu);
    set_hash_define_int(replacements + nrepl++, "GF2X_MUL_TOOMU_ALWAYS_THRESHOLD", tua);
    ulim = tua;
}

void digest_table()
{
    unsigned long i;

    if (blim || ulim) {
        set_hash_define(replacements + nrepl++,
                "GF2X_TOOM_TUNING_INFO", toom_info_string);

        if (blim) {
            printf("toom_table computed up to size %lu\n", blim);
            fill_holes(best_tab, blim, GF2X_SELECT_TC4);
            get_toom_thresholds();
            prepare_and_push_hash_define_tbl("GF2X_BEST_TOOM_TABLE", best_tab, blim);
        }

        if (ulim) {
            printf("utoom_table computed up to size %lu\n", ulim);
            fill_holes(best_utab, ulim, GF2X_SELECT_UNB_TC3U);
            get_utoom_thresholds();
            prepare_and_push_hash_define_tbl("GF2X_BEST_UTOOM_TABLE",best_utab,ulim);
        }
    }

    if (nstairs) {
        set_hash_define(replacements + nrepl++,
                "GF2X_FFT_TUNING_INFO", fft_info_string);
        qsort(fft_tbl, nstairs, sizeof(struct fft_tuning_pair),
                (sortfunc_t) &compare_fft_tuning);

        /* coalesce fft results */
        unsigned long j = 0;
        for(i = 0 ; i < nstairs ; ) {
            unsigned long k;
            for(k = i ; k < nstairs ; k++) {
                if (fft_tbl[k].method != fft_tbl[i].method)
                    break;
            }
            memcpy(&(fft_tbl[j]),&(fft_tbl[i]),sizeof(struct fft_tuning_pair));
            j++;
            i = k;
        }
        nstairs = j;
        prepare_and_push_hash_define_fft_tbl("GF2X_MUL_FFT_TABLE", fft_tbl, nstairs);
    }

    replace(replacements, nrepl, thresholds_file);
    for(i = 0 ; i < nrepl; i++ ) {
	free(replacements[i].identifier);
	free(replacements[i].string);
    }
}

void usage(int rc)
{
    FILE * f = rc ? stderr : stdout;

    fprintf(f, "Usage: modify-thresholds [-o <thresholds file>]\n");
    exit(rc);
}

int main(int argc, char *argv[])
{
    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        if (strcmp(argv[0], "--help") == 0) {
            usage(0);
        }
        if (strcmp(argv[0], "--output") == 0 || strcmp(argv[0], "-o") == 0) {
            argc--,argv++;
            if (! argc) {
                usage(1);
            }
            thresholds_file = argv[0];
            continue;
        }
        usage(1);
    }

    read_table();
    digest_table();

    fflush(stdout);

    return 0;
}

/* vim: set sw=4 sta et: */
