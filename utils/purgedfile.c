#include "cado.h"
#include <string.h>
#include <errno.h>
#include "purgedfile.h"
#include "timing.h"
#include "macros.h"
#include "gzip.h"

void purgedfile_stream_init(purgedfile_stream_ptr ps)
{
    memset(ps, 0, sizeof(purgedfile_stream));
    ps->t0 = ps->t1 = wct_seconds();
}

void purgedfile_stream_clear(purgedfile_stream_ptr ps)
{
    purgedfile_stream_closefile(ps);
    free(ps->cols); ps->cols = NULL;
}

void purgedfile_stream_openfile(purgedfile_stream_ptr ps, const char * fname)
{
    ASSERT_ALWAYS(ps->source == NULL);
    ps->source = fopen_compressed_r(fname, &(ps->pipe), NULL);
    if (ps->source == NULL) {
        fprintf(stderr, "opening %s: %s\n", fname, strerror(errno));
        exit(1);
    }
    int rc = fscanf (ps->source, "%d %d\n", &ps->nrows, &ps->ncols);
    ASSERT_ALWAYS(rc == 2);
}

void purgedfile_stream_closefile(purgedfile_stream_ptr ps)
{
    ASSERT_ALWAYS(ps->pipe != -1);
    if (ps->source) {
        if (ps->pipe) pclose(ps->source); else fclose(ps->source);
    }
    ps->source = NULL;
    ps->pipe = 0;
}

#if 0
void purgedfile_stream_bind(purgedfile_stream_ptr ps, FILE * f)
{
    ps->source = f;
    ps->pipe = -1;
}

void purgedfile_stream_unbind(purgedfile_stream_ptr ps)
{
    ps->source = NULL;
    ps->pipe = 0;
}
#endif


int purgedfile_stream_disp_progress_now_p(purgedfile_stream_ptr ps)
{
    if (ps->rrows % 100)
        return 0;
    double t = wct_seconds();
    if (ps->rrows % 1000000 == 0 || t >= ps->t1 + 1) {
        ps->dt = t - ps->t0;
        ps->mb_s = ps->dt > 0.01 ? (ps->pos/ps->dt * 1.0e-6) : 0;
        ps->rows_s = ps->dt > 0.01 ? ps->rrows/ps->dt : 0;
        ps->t1 = t;
        return 1;
    }
    return 0;
}

void purgedfile_stream_trigger_disp_progress(purgedfile_stream_ptr ps)
{
    double t = wct_seconds();
    ps->dt = t - ps->t0;
    ps->mb_s = ps->dt > 0.01 ? (ps->pos/ps->dt * 1.0e-6) : 0;
    ps->rows_s = ps->dt > 0.01 ? ps->rrows/ps->dt : 0;
    ps->t1 = t;
}

void purgedfile_stream_rewind(purgedfile_stream_ptr ps)
{
    ASSERT_ALWAYS(ps->source);
    int rc = fseek(ps->source, 0, SEEK_SET);
    ASSERT_ALWAYS(rc == 0);
    /* re-parse header */
    rc = fscanf (ps->source, "%d %d\n", &ps->nrows, &ps->ncols);
    ASSERT_ALWAYS(rc == 2);
}


static void purgedfile_stream_provision_for_primes(purgedfile_stream_ptr ps, int nc)
{
    if (nc > 0 && ps->nc_alloc < nc) {
        ps->nc_alloc = nc + nc / 2;
        ps->cols = (int *) realloc(ps->cols, ps->nc_alloc * sizeof(int));
    }
}

#define STORE(p, v)     do { if (p) *p++=(v); else (v);  nread++; } while (0)

int purgedfile_stream_get(purgedfile_stream_ptr ps, char * line)
{
    FILE * f = ps->source;

    static const int ugly[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
        -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    int c;
    char * p;
    
another_line:
    ps->lnum++;

    p = line;
    int nread = 0;

    STORE(p, c=fgetc(f));
    if (c == EOF) return -1;

    /* parse comments, if any */
    if (c == '#') {
        for( ; c != EOF && c != '\n' ; ) STORE(p, c=fgetc(f));
        goto another_line;
    }

    /* parse nodup_index */
    {
        int v, * w = &ps->nodup_index;
        for(*w = 0 ; (v=ugly[(unsigned char) c]) >= 0 ; ) {
            *w=*w*10+v;
            STORE(p, c=fgetc(f));
        }
    }
    ASSERT_ALWAYS(c == ' ');
    STORE(p, c=fgetc(f));

    /* parse a */
    {
        int64_t * w = &ps->a;
        int v, s = 1;
        if (c == '-') { s=-1; STORE(p, c=fgetc(f)); }
        for(*w = 0 ; (v=ugly[(unsigned char) c]) >= 0 ; ) {
            *w=*w*10+v;
            STORE(p, c=fgetc(f));
        }
        *w*=s;
    }
    ASSERT_ALWAYS(c == ' ');
    STORE(p, c=fgetc(f));

    /* parse b */
    {
        uint64_t * w = &ps->b;
        int v;
        for(*w = 0 ; (v=ugly[(unsigned char) c]) >= 0 ; ) {
            *w=*w*10+v;
            STORE(p, c=fgetc(f));
        }
    }
    ASSERT_ALWAYS(c == ' ');
    STORE(p, c=fgetc(f));

    if (!ps->parse_only_ab) {

        /* parse nc */
        {
            int v, * w = &ps->nc;
            for(*w = 0 ; (v=ugly[(unsigned char) c]) >= 0 ; ) {
                *w=*w*10+v;
                STORE(p, c=fgetc(f));
            }
        }
        /* don't advance the file pointer right now, do it from within
         * the loop instead.
         */

        purgedfile_stream_provision_for_primes(ps, ps->nc);

        for(int i = 0 ; i < ps->nc ; i++) {
            /* parse prime */
            ASSERT_ALWAYS(c == ' ');
            STORE(p, c=fgetc(f));
            {
                int v, * w = &ps->cols[i];
                for(*w = 0 ; (v=ugly[(unsigned char) c]) >= 0 ; ) {
                    *w=*w*16+v;
                    STORE(p, c=fgetc(f));
                }
            }
        }
        ASSERT_ALWAYS(c == '\n');
    }

    /* skip rest of line -- a no-op if we've been told to parse
     * everything. */
    for( ; c != EOF && c != '\n' ; ) STORE(p, c=fgetc(f));

    if (p) *p++='\0';   // don't update nread for this one.

    ps->pos += nread;
    ps->rrows++;

    return nread;
}
#undef  STORE

