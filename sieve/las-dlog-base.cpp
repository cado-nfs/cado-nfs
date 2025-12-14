#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdint>
#include <cstdio>

#include "las-dlog-base.hpp"
#include "las-multiobj-globals.hpp"
#include "macros.h"
#include "typedefs.h"
#include "params.h"
#include "portability.h"
#include "verbose.h"

void las_dlog_base::declare_usage(cxx_param_list & pl)
{
    if (!dlp_descent) return;
    param_list_decl_usage(pl, "renumber", "renumber table (for the descent)");
    param_list_decl_usage(pl, "log", "log table, as built by reconstructlog");
    /* These belong to las-siever-config of course. But we do a lookup
     * from here as well.
     */
    param_list_decl_usage(pl, "lpb0", "set large prime bound on side 0 to 2^lpb0");
    param_list_decl_usage(pl, "lpb1", "set large prime bound on side 1 to 2^lpb1");
}

las_dlog_base::las_dlog_base(cxx_cado_poly const & cpoly, cxx_param_list & pl)
    : cpoly(cpoly)
    , renumber_table(cpoly)
{
    if (!dlp_descent) return;
    renumberfilename = NULL;
    logfilename = NULL;
    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "renumber")) != NULL) {
        renumberfilename = strdup(tmp);
    }
    if ((tmp = param_list_lookup_string(pl, "log")) != NULL) {
        logfilename = strdup(tmp);
    }
    if (!logfilename != !renumberfilename) {
        fprintf(stderr, "In descent mode, want either renumber+log, or none\n");
        exit(EXIT_FAILURE);
    }
    if (!param_list_parse_ulong(pl, "lpb0", &(lpb[0]))) {
        fprintf(stderr, "In descent mode, want lpb0 for the final descent\n");
        exit(EXIT_FAILURE);
    }
    if (!param_list_parse_ulong(pl, "lpb1", &(lpb[1]))) {
        fprintf(stderr, "In descent mode, want lpb1 for the final descent\n");
        exit(EXIT_FAILURE);
    }
    read();
}

las_dlog_base::~las_dlog_base()
{
    if (!dlp_descent) return;
    free(renumberfilename);
    free(logfilename);
}

bool las_dlog_base::is_known(int side, p_r_values_t p, p_r_values_t r) const
{
    ASSERT_ALWAYS(dlp_descent);
    // if p is above large prime bound,  its log is not known.
    if (lpb[side] >= 64)
        return false;
    if (p >> lpb[side]) {
        return false;
    }
    if (renumberfilename) {
        /* For now we assume that we know the log of all bad ideals */
        /* If we want to be able to do a complete lookup for bad ideals,
         * then we need to use
         *      renumber_table.indices_from_p_a_b
         * which needs a,b as well.
         */
        if (renumber_table.is_bad ({ p, r, side }))
            return true;
        index_t const h = renumber_table.index_from_p_r({ p, r, side });
        return known_logs[h];
    }
    return true;
}


void las_dlog_base::read()
{
    ASSERT_ALWAYS(dlp_descent);
    if (!renumberfilename) {
        verbose_fmt_print(0, 1, "# Descent: no access to renumber table given, using lpb({}/{}) to decide what are the supposedly known logs\n",
                lpb[0], lpb[1]);
        return;
    }

    verbose_fmt_print(0, 1, "# Descent: will get list of known logs from {}, using also {} for mapping\n", logfilename, renumberfilename);

    renumber_table.read_from_file(renumberfilename, 1);

    for(int side = 0 ; side < renumber_table.get_nb_polys() ; side++) {
        if (lpb[side] != renumber_table.get_lpb(side)) {
            fmt::print(stderr, "lpb{}={} different from lpb{}={} stored in renumber table, probably a bug\n", side, lpb[side], side, renumber_table.get_lpb(side));
            exit(EXIT_FAILURE);
        }
    }

    uint64_t const nprimes = renumber_table.size();
    known_logs.assign(nprimes + 32, false);
    /* 32 is because the SM columns are here, too ! We would like to
     * avoid reallocation, so let's be generous (anyway we'll
     * reallocate if needed)
     */

    /* format of the log table: there are FIVE different line types.
     *
     * [index] added column [log]
     * [index] bad ideals [log]
     * [index] [p] [side] rat [log]
     * [index] [p] [side] [r] [log]
     * [index] SM col [i] [log]
     * where by default side=0 means the rational side, side=1 the algebraic one
     * r is the root of f(x) mod p, where f is the algebraic polynomial
     * log is the virtual logarithm
     * (see https://sympa.inria.fr/sympa/arc/cado-nfs/2019-02/msg00010.html)
     *
     * Here we care only about the index anyway. By the way, this index
     * is written in hex.
     */
    FILE * f = fopen(logfilename, "r");
    ASSERT_ALWAYS(f != NULL);
    size_t nlogs=0;
    for(int lnum = 0 ; ; lnum++) {
        char line[1024];
        char * x = fgets(line, sizeof(line), f);
        if (x == NULL) break;
        for( ; *x && isspace(*x) ; x++) ;
        if (*x == '#') continue;
        if (!*x) continue;
        errno=0;
        unsigned long const z = strtoul(x, &x, 16);
        if (errno) {
            fprintf(stderr, "Parse error at line %d in %s: %s\n", lnum, logfilename, line);
            break;
        }
        if (z >= known_logs.size()) {
            /* happens for SM columns: we have more than the size of the
             * renumber table ! */
            known_logs.resize(z + 1);
        }
	nlogs+=!known_logs[z];
        known_logs[z] = true;
    }
    fclose(f);
    verbose_fmt_print(0, 1, "# Got {} known logs from {}\n", nlogs, logfilename);
}


