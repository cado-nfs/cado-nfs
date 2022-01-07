#include "cado.h" // IWYU pragma: keep
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <gmp.h>

#include "macros.h"
#include "mf.h"
#include "fix-endianness.h" // fwrite32_little
#include "misc.h"       // derived_filename

/* This our favorite way of expanding buffer. We always clear to zero */
static inline void expand(struct mf_io_file * m, uint64_t nn)
{
    m->p = realloc(m->p, nn * sizeof(uint32_t));
    ASSERT_ALWAYS(m->p);
    if (nn > m->alloc)
        memset(m->p + m->alloc, 0, (nn - m->alloc) * sizeof(uint32_t));
    m->alloc = nn;
}

static void abort_unexpected_eof()
{
    fprintf(stderr, "Unexpected EOF -- (input file truncated ?).\n");
    exit(1);
}

char * build_mat_auxfile(const char * prefix, const char * what, const char * ext)
{
    /* This function gets called only if we have automatically decided on
     * the prefix, which means that normally, the decision has been taken
     * in accordance to the type of the matrix file.
     *
     * We check for a match between the elcted prefix and the matrix
     * filename. If found, that match is chopped.
     */
    return derived_filename(prefix, what, ext);
}

int matrix_autodetect_input(struct mf_io_file * m_in, const char * mfile)
{
    const char * file_ext[2] = { ".txt", ".bin" };

    if (mfile == NULL)
        return -1;

    /* try to auto-detect */
    if (has_suffix(mfile, file_ext[0])) {
        m_in->ascii = 1;
    } else if (has_suffix(mfile, file_ext[1])) {
        m_in->ascii = 0;
    } else {
        /* for *regular* files, where it's possible to read somewhat
         * ahead of time, try to auto-detect based on the contents */
        ASSERT_ALWAYS(m_in->f);
        struct stat sbuf[1];
        int rc = fstat(fileno(m_in->f), sbuf);
        DIE_ERRNO_DIAG(rc < 0, "fstat(%s)", mfile);
        if (!S_ISREG(sbuf->st_mode)) {
            // guard against tricks like /dev/fd/ to unseekable fd's.
            return -1;
        }
        char test[1024];
        int n = fread(test, 1, 1024, m_in->f);
        DIE_ERRNO_DIAG(n < 1024 && !feof(m_in->f), "fread(%s)", mfile);
        int k;
        for(k = 0 ; k < n && (isdigit((int)(unsigned char)test[k]) || isspace((int)(unsigned char)test[k])) ; k++);
        if (k < n) {
            // assume binary.
            m_in->ascii = 0;
        } else {
            m_in->ascii = 1;
        }
        rc = fseek(m_in->f, 0L, SEEK_SET);
        DIE_ERRNO_DIAG(rc < 0, "rewind(%s)", mfile);
        fprintf(stderr, "auto-detected %s as %s based on contents\n",
                mfile, m_in->ascii ? "ascii" : "binary");
    }
    return 0;
}

void matrix_read_pass(
        struct mf_io_file * m_in,
        struct mf_io_file * m_out,
        struct mf_io_file * rw_out,
        struct mf_io_file * cw_out,
        unsigned int rskip,
        unsigned int cskip,
        int progress,
        /* number of 32-bit words used to store (signed) coefficients */
        int withcoeffs
        )
{
    /* This temp buffer is used only if skipping columns and outputting a
     * matrix file -- and even then, we can avoid it if the matrix rows are
     * known to be sorted (although we don't do so at the moment). */
    struct mf_io_file temp[1];
    memset(temp, 0, sizeof(struct mf_io_file));

    /* This gathers several possible ways of reading the matrix. It
     * encompasses the task of both mf_b_fq and mf_a2b_fq. It can also be used
     * to write anf mf_b2a tool.
     *
     * This supersedes the read_easy functions which used to reside in
     * ../readmat-easy.c -- see build.c for an example of use (the
     * function there could be exported somewhere).
     */
    ASSERT_ALWAYS(m_in);
    ASSERT_ALWAYS(m_in->p == NULL || m_in->f == NULL);
    ASSERT_ALWAYS(m_out == NULL || m_out->p == NULL || m_out->f == NULL);

    int rc;
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;

    /* In some cases we have to seek in the input file. If turns out to
     * be impossible, then avoid. */
    int can_seek = 1;

    unsigned int exp_nr=0;
    unsigned int exp_nc=0;
    int drop_cw_p = 0;
    if (m_in->f && m_in->ascii) {
        rc = fscanf(m_in->f, "%u %u", &exp_nr, &exp_nc);
        if (rc == EOF) abort_unexpected_eof();
        if (cw_out && cw_out->p == NULL)  {
            ASSERT_ALWAYS(exp_nc >= cskip);
#ifdef __COVERITY__
            __coverity_mark_pointee_as_sanitized(&exp_nc, ALLOCATION);
#endif
            expand(cw_out, exp_nc - cskip);
            drop_cw_p=1;
        }
        // NOTE: here we have the chance to know the number of columns we
        // expect. So we'll be able to write the cw file with the proper
        // number of columns. On the contrary, when reading from binary
        // data, we can't compute cw precisely. Therefore, and I admit
        // it's quite odd, we have a behaviour which diverges depending
        // on whether we're ascii or not.
        if (cw_out)
            cw_out->size = exp_nc;
    } else {
        if (cw_out && cw_out->p == NULL)  {
            expand(cw_out, 1000 * 1000);
            drop_cw_p=1;
        }
    }
    unsigned int i;

    uint32_t * ptr = NULL;
    if (!m_in->f) {
        ptr = m_in->p;
    }

    mpz_t large_coeff;
    mpz_init(large_coeff);
    uint32_t * large_coeff_buffer = NULL;
    if (withcoeffs > 1) {
        large_coeff_buffer = malloc(withcoeffs * sizeof(uint32_t));
    }

    size_t readbytes = 0;
    int delta;

    for(i = 0 ; ; i++) {
        uint32_t j;
        /* {{{ read row weight w */
        uint32_t w;
        if (m_in->f) {
            if (m_in->ascii) {
                delta=0;
                rc = fscanf(m_in->f, "%" SCNu32"%n", &w, &delta);
                readbytes += delta;
                if (i >= exp_nr) {
                    if (rc == EOF)
                        break;
                    if (i == exp_nr) {
                        fprintf(stderr,
                                "Warning: reading extra rows at end of file\n");
                    }
                }
                if (rc == EOF) abort_unexpected_eof();
            } else {
                rc = fread32_little(&w, 1, m_in->f);
                if (rc == 0 && feof(m_in->f)) break;
                readbytes += rc?4:0;
            }
        } else {
            if (((uint64_t) (ptr - m_in->p) >= m_in->size))
                break;
            w = *ptr++;
        }
        /* }}} */

        uint32_t ww = w + withcoeffs * w;

        /* reading the row weight to the rw file is deferred to after the
         * point where we've cleaned the columns to be skipped. */

        if ((!m_out && !cw_out) || (i < rskip)) {
            /* {{{ skip fast over. */
            if (m_in->f) {
                if (!m_in->ascii) {
                    if (can_seek) {
                        rc = fseek(m_in->f, ww * sizeof(uint32_t), SEEK_CUR);
                        if (rc != 0)
                            can_seek = 0;
                    }
                    if (!can_seek) {
                        for(j = 0 ; j < ww ; j++) {
                            uint32_t c;
                            rc = fread32_little(&c, 1, m_in->f);
                            if (!rc) abort_unexpected_eof();
                            readbytes += rc?4:0;
                        }
                    }
                } else {
                    /* TODO: do better */
                    for(j = 0 ; j < w ; j++) {
                        uint32_t c;
                        delta=0;
                        rc = fscanf(m_in->f, "%" SCNu32"%n", &c, &delta);
                        readbytes += delta;
                        if (rc == EOF) abort_unexpected_eof();
                        if (withcoeffs == 1) {
                            /* XXX Can probably be more tolerant on how
                             * coeffs are presented */
                            int sp = fgetc(m_in->f);
                            ASSERT_ALWAYS(sp == ' ' || sp == ':');
                            delta=0;
                            rc = fscanf(m_in->f, "%" SCNd32"%n", (int32_t *) &c, &delta);
                            readbytes += delta;
                            if (rc == EOF) abort_unexpected_eof();
                        } else if (withcoeffs > 1) {
                            int sp = fgetc(m_in->f);
                            ASSERT_ALWAYS(sp == ' ' || sp == ':');
                            readbytes++;
                            /* suppressed */
                            delta = 0;
                            rc = gmp_fscanf(m_in->f, "%*Zd%n", &delta);
                            if (rc == EOF) abort_unexpected_eof();
                            readbytes += delta;
                        }
                    }
                }
            } else {
                ptr += w;
            }
            /* }}} */
            goto row_done;
        }

        uint32_t * ptr_w = NULL;

        if (m_out) { /* {{{ start row for matrix file */
            if (m_out->f) {
                if (!cskip) {
                    if (m_out->ascii) {
                        fprintf(m_out->f, "%" PRIu32, w);
                    } else {
                        rc = fwrite32_little(&w, 1, m_out->f);
                    }
                }
                /* Otherwise we have to defer */
            } else {
                if (m_out->size + ww >= m_out->alloc) {
                    expand(m_out, (1 << 20) + (m_out->size + ww) * 1.05);
                }
                ptr_w = &(m_out->p[m_out->size]);
                *ptr_w = w;
                m_out->size++;
            }
        } /* }}} */

        uint32_t real_weight = 0;

        if (m_out && cskip && temp->alloc < ww) {
            expand(temp, ww + ww/2);
        }

        for(j = 0 ; j < w ; j++) {
            /* {{{ read this column index */
            uint32_t c;
            int32_t coeff=0;
            if (m_in->f) {
                if (m_in->ascii) {
                    delta = 0;
                    rc = fscanf(m_in->f, "%" SCNu32"%n", &c, &delta);
                    if (rc == EOF) abort_unexpected_eof();
                    readbytes += delta;
                    if (withcoeffs == 1) {
                        int sp = fgetc(m_in->f);
                        ASSERT_ALWAYS(sp == ' ' || sp == ':');
                        readbytes++;
                        delta = 0;
                        rc = fscanf(m_in->f, "%" SCNd32"%n", (int32_t*) &coeff, &delta);
                        if (rc == EOF) abort_unexpected_eof();
                        readbytes += delta;
                    } else if (withcoeffs > 1) {
                        int sp = fgetc(m_in->f);
                        ASSERT_ALWAYS(sp == ' ' || sp == ':');
                        readbytes++;
                        delta = 0;
                        rc = gmp_fscanf(m_in->f, "%Zd%n", large_coeff,&delta);
                        if (rc == EOF) abort_unexpected_eof();
                        readbytes += delta;
                        uint32_t * q = large_coeff_buffer;
                        memset(q, 0, withcoeffs * sizeof(uint32_t));
                        size_t countp;
                        mpz_export(q, &countp, -1, sizeof(uint32_t), -1, 0, large_coeff);
                        ASSERT_ALWAYS(countp <= (size_t) withcoeffs);
                    }
                } else {
                    rc = fread32_little(&c, 1, m_in->f);
                    if (!rc) abort_unexpected_eof();
                    readbytes += rc?4:0;
                    if (withcoeffs == 1) {
                        rc = fread32_little((uint32_t*) &coeff, 1, m_in->f);
                        if (!rc) abort_unexpected_eof();
                        readbytes += rc?4:0;;
                    } else {
                        rc = fread(large_coeff_buffer, sizeof(uint32_t), withcoeffs, m_in->f);
                        if (rc != withcoeffs) abort_unexpected_eof();
                        readbytes += rc?4:0;;
                    }
                }
            } else {
                c = *ptr++;
                if (withcoeffs == 1) {
                    coeff = *ptr++;
                } else if (withcoeffs > 1) {
                    memcpy(large_coeff_buffer, ptr, withcoeffs * sizeof(uint32_t));
                    ptr += withcoeffs;
                }
            }
            /* }}} */
            if (cskip) {
                if (c < cskip) {
                    continue;   /* next coefficient */
                } else {
                    if (m_out) {
                        temp->p[(1 + withcoeffs) * real_weight] = c - cskip;
                    }
                    if (withcoeffs == 1) {
                        temp->p[2 * real_weight + 1] = coeff;
                    } else if (withcoeffs > 1) {
                        uint32_t * q = temp->p + (1+withcoeffs)*real_weight + 1; 
                        memcpy(q, large_coeff_buffer, withcoeffs * sizeof(uint32_t));
                    }
                    real_weight++;
                }
            } else {
                real_weight++;
            }

            if (m_out) { /* {{{ write to matrix file */
                if (m_out->f) {
                    if (!cskip) {
                        if (m_out->ascii) {
                            fprintf(m_out->f, " %" PRIu32, c);
                            if (withcoeffs == 1) {
                                fprintf(m_out->f, ":%" PRId32, coeff);
                            } else if (withcoeffs > 1) {
                                mpz_import(large_coeff, withcoeffs, sizeof(uint32_t), -1, -1, 0, large_coeff_buffer);
                                gmp_fprintf(m_out->f, ":%Zd", large_coeff);
                            }
                        } else {
                            rc = fwrite32_little(&c, 1, m_out->f);
                            if (withcoeffs == 1) {
                                rc = fwrite32_little((uint32_t*)&coeff, 1, m_out->f);
                            } else if (withcoeffs > 1) {
                                uint32_t * q = large_coeff_buffer;
                                rc = fwrite(q, sizeof(uint32_t), withcoeffs, m_out->f);
                                if (rc != withcoeffs) abort_unexpected_eof();

                            }
                        }
                    }
                    // otherwise it's too early, size hasn't been written yet.
                } else {
                    m_out->p[m_out->size++] = c-cskip;
                    if (withcoeffs == 1) {
                        m_out->p[m_out->size++] = coeff;
                    } else if (withcoeffs > 1) {
                        memcpy(m_out->p + m_out->size, large_coeff_buffer, withcoeffs * sizeof(uint32_t));
                        m_out->size += withcoeffs;
                    }
                }
            }/*}}}*/

            if (exp_nc && c >= exp_nc) {
                fprintf(stderr, "Warning: column index %" PRIu32 " exceeds header value %u\n", c, exp_nc);
            }
            if (cw_out) {
                if (c-cskip >= cw_out->alloc) {
                    expand(cw_out, c-cskip + cw_out->alloc / 2);
                }
                if (c-cskip >= cw_out->size) {
                    cw_out->size = c-cskip + 1;
                }
                cw_out->p[c-cskip]++;
            }
        }
        w = real_weight;
        ww = w + withcoeffs * w;

        if (m_out && cskip) {/* {{{ start row for matrix file */
            if (m_out->f) {
                if (m_out->ascii) {
                    fprintf(m_out->f, "%" PRIu32, w);
                    if (withcoeffs == 1) {
                        for(j = 0 ; j < ww ; ) {
                            fprintf(m_out->f, " %" PRIu32, temp->p[j++]);
                            fprintf(m_out->f, ":%" PRId32, temp->p[j++]);
                        }
                    } else if (withcoeffs > 1) {
                        for(j = 0 ; j < ww ; ) {
                            fprintf(m_out->f, " %" PRIu32, temp->p[j++]);
                            uint32_t * q = temp->p;
                            mpz_import(large_coeff, withcoeffs, sizeof(uint32_t), -1, -1, 0, q);
                            gmp_fprintf(m_out->f, ":%Zd", large_coeff);
                            j += withcoeffs;
                        }
                    } else {
                        for(j = 0 ; j < w ; j++) {
                            fprintf(m_out->f, " %" PRIu32, temp->p[j]);
                        }
                    }
                } else {
                    rc = fwrite32_little(&w, 1, m_out->f);
                    rc = fwrite32_little(temp->p, ww, m_out->f);
                }
            } else {
                // fix this up. It's the only thing we have to do.
                *ptr_w = w;
            }
        } /* }}} */

        /* row done */
        if (m_out && m_out->f && m_out->ascii)
            fprintf(m_out->f, "\n");

        if (rw_out) { /* {{{ write row weight to rw file */
            rw_out->size++;
            if (rw_out->p) rw_out->p[i]=w;
            if (rw_out->f) {
                if (rw_out->ascii) {
                    fprintf(rw_out->f, "%" PRIu32 "\n", w);
                } else {
                    rc = fwrite32_little(&w, 1, rw_out->f);
                }
            }
        }/*}}}*/

row_done:
        if (progress) {/*{{{*/
            time_t t = time(NULL);
            if (t >= t1) {
                t1 = t + 1;
                int dt = t - t0;
                size_t sz;
                long o = -1;
                if (m_in->f && (o = ftell(m_in->f)) >= 0) {
                    sz = o;
                } else {
                    sz = (ptr - m_in->p) * sizeof(uint32_t);
                }
                fprintf(stderr, "read %zu MB in %d s (%.1f MB/s) %zu %zu  \r",
                        readbytes >> 20,
                        dt,
                        1.0e-6 * (double) readbytes / dt, readbytes, sz);
            }
        }/*}}}*/
    }

    mpz_clear(large_coeff);
    free(large_coeff_buffer);

    if (progress)
        fprintf(stderr, "\n");

    /* write cw if requested *//*{{{*/
    if (cw_out) {
        if (cw_out->f) {
            if (cw_out->ascii) {
                uint32_t j;
                for(j = 0 ; j < cw_out->size ; j++) {
                    fprintf(cw_out->f, "%" PRIu32 "\n", cw_out->p[j]);
                }
            } else {
                rc = fwrite32_little(cw_out->p, cw_out->size, cw_out->f);
            }
        }
    }
    if (drop_cw_p) {
        free(cw_out->p);
        cw_out->p = NULL;
    }/*}}}*/

    if (exp_nr && i != exp_nr) {
        fprintf(stderr,
                "Warning: counted %u rows,"
                " not according with header value %u\n", i, exp_nr);
    }
}

