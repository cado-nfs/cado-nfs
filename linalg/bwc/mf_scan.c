#include "cado.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "macros.h"
#include "portability.h"
#include "utils.h"
#include "mf.h"


void usage(int rc)
{
    fprintf(stderr,
    "This program make one reading pass through a matrix, producing one or\n"
    "several output files. Depending on the arguments given, the behaviour\n"
    "differs.\n"
    "\n"
    "read an ascii matrix, produce a binary one as well as the companion\n"
    ".rw and .cw files:\n"
    "\n"
    "./mf-scan --ascii-in mfile=<infile> --binary-out ofile=<outfile> --freq\n"
    "\n"
    "read a binary matrix, produce an ascii one, and no frequency files.\n"
    "(NOTE though that the produced ascii file lacks a header !)\n"
    "\n"
    "./mf-scan --binary-in mfile=<infile> --ascii-out ofile=<outfile> --nofreq\n"
    "\n"
    "produce only the frequency files corresponding to an existing binary\n"
    "matrix.\n"
    "\n"
    "./mf-scan --binary-in mfile=<infile>\n"
    "\n"
    "\n"
            "Recognized options"
                " (<option_name>=<value>, or --<option_name> <value>:\n"
            " mfile       matrix file (can also be given freeform)\n"
            " rwfile      row weight file\n"
            " cwfile      col weight file\n"
            " ofile       output file name ; no output if unspecified\n"
            "Recognized flags:"
            " --ascii-in   mfile in ascii\n"
            " --ascii-out  ofile in ascii\n"
            " --binary-in  mfile in binary\n"
            " --binary-out ofile in binary\n"
            " --quiet      no progress info\n"
            " --nofreq     no .rw and .cw files\n"
            " --freq       output .rw and .cw files\n"
            " --ascii-freq .rw and .cw files in ascii\n"
            " --binary-freq .rw and .cw files in binary\n"
    );
    exit(rc);
}

const char * file_ext[2] = { ".txt", ".bin" };

static void autodetection_error(const char * what)
{
    fprintf(stderr, "Error: cannot auto-detect format for %s. Please use the options.\n", what);
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;

    const char * rwfile = NULL;
    const char * cwfile = NULL;
    const char * mfile = NULL;
    const char * ofile = NULL;

    unsigned int wild =  0;

    /* many flags -- several are mutually incompatible */
    unsigned int rskip = 0;
    unsigned int cskip = 0;
    int quiet =  0;

    int ascii_in = 0;
    int ascii_out = 0;
    int binary_in = 0;
    int binary_out = 0;

    int freq = 0;
    int nofreq = 0;
    int ascii_freq = 0;
    int binary_freq = 0;

    int withcoeffs = 0;

    param_list_init(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "--quiet", &quiet);
    param_list_configure_switch(pl, "--ascii-in", &ascii_in);
    param_list_configure_switch(pl, "--binary-in", &binary_in);
    param_list_configure_switch(pl, "--ascii-out", &ascii_out);
    param_list_configure_switch(pl, "--binary-out", &binary_out);
    param_list_configure_switch(pl, "--ascii-freq", &ascii_freq);
    param_list_configure_switch(pl, "--binary-freq", &binary_freq);
    param_list_configure_switch(pl, "--nofreq", &nofreq);
    param_list_configure_switch(pl, "--freq", &freq);
    param_list_configure_switch(pl, "--withcoeffs", &withcoeffs);

    for(;argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild == 0) {
            mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", argv[0]);
        usage(1);
    }

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        cwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "ofile")) != NULL) {
        ofile = tmp;
    }

    param_list_parse_uint(pl, "rskip", &rskip);
    param_list_parse_uint(pl, "cskip", &cskip);

    if (ascii_in && binary_in) usage(1);
    if (ascii_out && binary_out) usage(1);
    if (ascii_freq && binary_freq) usage(1);
    if (ascii_freq || binary_freq) freq=1;
    if (rwfile || cwfile) freq=1;
    if (freq && nofreq) usage(1);
    if ((ascii_out || binary_out) && !ofile) usage(1);

    const char * freq_default_prefix = NULL;

    /* Start with the input */

    struct mf_io_file m_in[1];
    memset(m_in, 0, sizeof(struct mf_io_file));

    if (!mfile) usage(1);
    /* special case to mean stdin */
    if (strcmp(mfile, "-") == 0) {
        mfile = NULL;
        m_in->f = stdin;
    } else {
#ifdef  HAVE_MINGW
        if (ascii_in) {
            m_in->f = fopen(mfile, "r");
        } else if (binary_in) {
            m_in->f = fopen(mfile, "rb");
        } else {
            fprintf(stderr, "Under MinGW, please specify explicitly --binary-in or --ascii-in\n");
            exit(1);
        }
#else
        /* Then we probably don't care */
        m_in->f = fopen(mfile, "rb");
#endif
        if (m_in->f == NULL) { perror(mfile); exit(1); }
    }
    if (mfile) freq_default_prefix=mfile;

    if (ascii_in || binary_in) {
        m_in->ascii = ascii_in;
    } else if (matrix_autodetect_input(m_in, mfile) < 0) {
        autodetection_error("input");
    }

    /* Now the output, if there is one. */
    struct mf_io_file m_out[1];
    memset(m_out, 0, sizeof(struct mf_io_file));
    
    if (ofile) {
        if (strcmp(ofile, "-") == 0) {
            ofile = NULL;
            m_out->f = stdout;
        } else {
#ifdef  HAVE_MINGW
            if (ascii_out) {
                m_out->f = fopen(ofile, "w");
            } else if (binary_out) {
                m_out->f = fopen(ofile, "wb");
            } else {
                fprintf(stderr, "Under MinGW, please specify explicitly --binary-out or --ascii-out\n");
                exit(1);
            }
#else
            /* Then we probably don't care */
            m_out->f = fopen(ofile, "wb");
#endif
        }
        if (ascii_out || binary_out) {
            m_out->ascii = ascii_out;
        } else if (ofile) {
            /* try to auto-detect */
            if (has_suffix(ofile, file_ext[0])) m_out->ascii = 1;
            else if (has_suffix(ofile, file_ext[1])) { m_out->ascii = 0; }
            else autodetection_error("output");
        } else {
            autodetection_error("output");
        }
        if (ofile) freq_default_prefix=ofile;
    }

    struct mf_io_file rw[1];
    struct mf_io_file cw[1];
    memset(rw, 0, sizeof(struct mf_io_file));
    memset(cw, 0, sizeof(struct mf_io_file));
    /* nofreq is the default */
    if (freq) {
        if (!rwfile || !cwfile) {
            if (!freq_default_prefix) {
                fprintf(stderr, "Cannot figure out names for .rw and .cw files\n");
                exit(1);
            }
        }
        if (ascii_freq || binary_freq) rw->ascii=ascii_freq;
        else if (rwfile && has_suffix(rwfile, file_ext[0])) rw->ascii = 1;
        else if (rwfile && has_suffix(rwfile, file_ext[1])) rw->ascii = 0;
        else rw->ascii=ofile ? m_out->ascii : m_in->ascii;

        if (ascii_freq || binary_freq) cw->ascii=ascii_freq;
        else if (cwfile && has_suffix(cwfile, file_ext[0])) cw->ascii = 1;
        else if (cwfile && has_suffix(cwfile, file_ext[1])) cw->ascii = 0;
        else cw->ascii=ofile ? m_out->ascii : m_in->ascii;

        if (!rwfile) {
            char * leakme;
            rwfile = leakme = build_mat_auxfile(freq_default_prefix, "rw", file_ext[!rw->ascii]);
        }

        if (!cwfile) {
            char * leakme;
            cwfile = leakme = build_mat_auxfile(freq_default_prefix, "cw", file_ext[!cw->ascii]);
        }
    }

    if (rwfile) rw->f = fopen(rwfile, "wb");
    if (cwfile) cw->f = fopen(cwfile, "wb");

    matrix_read_pass(m_in,
            m_out->f ? m_out : NULL,
            rw->f ? rw : NULL,
            cw->f ? cw : NULL,
            rskip,
            cskip,
            !quiet,
            withcoeffs);

    if (m_out->f && m_out->ascii) {
        fprintf(stderr, "\nWarning: output matrix is headerless\n");
    }


    if (rwfile) fclose(rw->f);
    if (cwfile) fclose(cw->f);
    if (mfile) fclose(m_in->f);
    if (ofile) fclose(m_out->f);

    param_list_clear(pl);
    return 0;
}
