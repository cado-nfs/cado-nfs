#include "cado.h" // IWYU pragma: keep
#include <cmath>       // for ldexp, sqrt
#include <cstdlib>     // for malloc, free, atoi, calloc
#include <cstring>     // for strcmp, strlen, strncpy
#include <regex.h>      // for regmatch_t, regcomp, regexec, regfree, REG_EX...
#include "facul_strategies.hpp"
#include "pm1.h"        // for pm1_plan_t, pm1_clear_plan, pm1_make_plan
#include "pp1.h"        // for pp1_plan_t, pp1_clear_plan, pp1_make_plan
#include "macros.h"
#include "verbose.h"

//#define USE_LEGACY_DEFAULT_STRATEGY 1

/*{{{*/
static int nb_curves90 (const unsigned int lpb);
#if 0
static int nb_curves95 (const unsigned int lpb);
static int nb_curves99 (const unsigned int lpb);
#endif

/* don't use nb_curves90, nb_curves95 and nb_curves99, only use nb_curves */
int nb_curves (const unsigned int lpb, const unsigned int mfb)
{
  /* for up to 2 large primes, we use a very large number of curves (say 100),
     since the more curves we need, the more likely the cofactor is smooth. */
  if (mfb <= 2 * lpb)
    return 100;
  return nb_curves90 (lpb);
}

static int
nb_curves_with_fbb (const unsigned long fbb,
		    const unsigned int lpb, const unsigned int mfb)
{
  /* if 2^mfb < fbb^2, we can have only one large prime */
  return (ldexp (1.0, mfb) < (double) fbb * (double) fbb)
    ? 0 : nb_curves (lpb, mfb);
}

static int
nb_curves90 (const unsigned int lpb)
{
  /* The following table, computed with do_table(10,33,ntries=10000) from the
     facul.sage file, ensures a probability of at least about 90% to find a
     factor below 2^lpb with n = T[lpb]. */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10-19 */
	     /* lpb=20 */ 1, /* 0:0.878100, 1:0.969600 */
	     /* lpb=21 */ 1, /* 1:0.940400 */
	     /* lpb=22 */ 1, /* 1:0.907400 */
	     /* lpb=23 */ 2, /* 1:0.859100, 2:0.903900 */
	     /* lpb=24 */ 4, /* 3:0.897000, 4:0.929200 */
	     /* lpb=25 */ 5, /* 4:0.884100, 5:0.916400 */
	     /* lpb=26 */ 6, /* 5:0.868600, 6:0.904200 */
	     /* lpb=27 */ 8, /* 7:0.873500, 8:0.901700 */
             /* lpb=28 */ 11, /* 10:0.896600, 11:0.918000 */
	     /* lpb=29 */ 13, /* 12:0.881700, 13:0.905600 */
	     /* lpb=30 */ 16, /* 15:0.889400, 16:0.909400 */
	     /* lpb=31 */ 18, /* 17:0.883500, 18:0.901300 */
	     /* lpb=32 */ 21, /* 20:0.884500, 21:0.903700 */
	     /* lpb=33 */ 25, /* 24:0.892800, 25:0.910000 */
        /* The extra ones below are computed with do_table(10,64,ntries=200) */
	     /* lpb=34 */ 24, /* 23:0.895000, 24:0.905000 */
	     /* lpb=35 */ 29, /* 28:0.895000, 29:0.905000 */
	     /* lpb=36 */ 30, /* 29:0.890000, 30:0.905000 */
	     /* lpb=37 */ 35, /* 34:0.890000, 35:0.900000 */
	     /* lpb=38 */ 37, /* 36:0.890000, 37:0.900000 */
	     /* lpb=39 */ 42, /* 41:0.890000, 42:0.900000 */
	     /* lpb=40 */ 45, /* 44:0.890000, 45:0.900000 */
	     /* lpb=41 */ 50, /* 49:0.885000, 50:0.905000 */
	     /* lpb=42 */ 56, /* 55:0.875000, 56:0.900000 */
	     /* lpb=43 */ 60, /* 59:0.895000, 60:0.900000 */
	     /* lpb=44 */ 68, /* 67:0.890000, 68:0.910000 */
	     /* lpb=45 */ 73, /* 72:0.895000, 73:0.905000 */
	     /* lpb=46 */ 82, /* 81:0.895000, 82:0.900000 */
	     /* lpb=47 */ 90, /* 89:0.895000, 90:0.900000 */
	     /* lpb=48 */ 93, /* 92:0.895000, 93:0.900000 */
	     /* lpb=49 */ 93, /* 93:0.900000 */
	     /* lpb=50 */ 111, /* 110:0.885000, 111:0.900000 */
	     /* lpb=51 */ 117, /* 116:0.890000, 117:0.905000 */
	     /* those below were computed with do_table(10,64,ntries=100) */
	     /* lpb=52 */ 130, /* 129:0.890000, 130:0.900000 */
	     /* lpb=53 */ 137, /* 136:0.890000, 137:0.900000 */
	     /* lpb=54 */ 142, /* 141:0.890000, 142:0.900000 */
	     /* lpb=55 */ 152, /* 151:0.890000, 152:0.920000 */
	     /* lpb=56 */ 167, /* 166:0.890000, 167:0.900000 */
	     /* lpb=57 */ 167, /* 167:0.910000 */
	     /* lpb=58 */ 186, /* 185:0.890000, 186:0.900000 */
	     /* lpb=59 */ 202, /* 201:0.880000, 202:0.900000 */
	     /* lpb=60 */ 227, /* 226:0.890000, 227:0.910000 */
	     /* lpb=61 */ 244, /* 243:0.890000, 244:0.910000 */
	     /* lpb=62 */ 244, /* 244:0.920000 */
	     /* lpb=63 */ 248, /* 247:0.880000, 248:0.900000 */
	     /* lpb=64 */ 294, /* 293:0.890000, 294:0.900000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}

#if 0
static int
nb_curves95 (const unsigned int lpb)
{
    /* same, but with target probability 95% */
    /* do_table(10,64,ntries=500,target_prob=0.95)
     */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
	     /* lpb=10 */ 0, /* 0:1.000000 */
	     /* lpb=11 */ 0, /* 0:1.000000 */
	     /* lpb=12 */ 0, /* 0:1.000000 */
	     /* lpb=13 */ 0, /* 0:1.000000 */
	     /* lpb=14 */ 0, /* 0:1.000000 */
	     /* lpb=15 */ 0, /* 0:0.998000 */
	     /* lpb=16 */ 0, /* 0:0.986000 */
	     /* lpb=17 */ 0, /* 0:0.972000 */
	     /* lpb=18 */ 0, /* 0:0.964000 */
	     /* lpb=19 */ 1, /* 0:0.926000, 1:0.986000 */
	     /* lpb=20 */ 1, /* 1:0.986000 */
	     /* lpb=21 */ 2, /* 1:0.946000, 2:0.976000 */
	     /* lpb=22 */ 3, /* 2:0.940000, 3:0.966000 */
	     /* lpb=23 */ 3, /* 3:0.952000 */
	     /* lpb=24 */ 5, /* 4:0.926000, 5:0.962000 */
	     /* lpb=25 */ 7, /* 6:0.934000, 7:0.956000 */
	     /* lpb=26 */ 8, /* 7:0.934000, 8:0.956000 */
	     /* lpb=27 */ 10, /* 9:0.936000, 10:0.956000 */
	     /* lpb=28 */ 13, /* 12:0.938000, 13:0.954000 */
	     /* lpb=29 */ 16, /* 15:0.946000, 16:0.956000 */
	     /* lpb=30 */ 18, /* 17:0.936000, 18:0.950000 */
	     /* lpb=31 */ 22, /* 21:0.934000, 22:0.958000 */
	     /* lpb=32 */ 26, /* 25:0.940000, 26:0.950000 */
	     /* lpb=33 */ 29, /* 28:0.942000, 29:0.952000 */
	     /* lpb=34 */ 33, /* 32:0.948000, 33:0.956000 */
	     /* lpb=35 */ 37, /* 36:0.938000, 37:0.952000 */
	     /* lpb=36 */ 42, /* 41:0.946000, 42:0.952000 */
	     /* lpb=37 */ 45, /* 44:0.946000, 45:0.958000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}

static int
nb_curves99 (const unsigned int lpb)
{
    /* same, but with target probability 99% */
    /* do_table(10,64,ntries=100,target_prob=0.99)
     */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
	     /* lpb=10 */ 0, /* 0:1.000000 */
	     /* lpb=11 */ 0, /* 0:1.000000 */
	     /* lpb=12 */ 0, /* 0:1.000000 */
	     /* lpb=13 */ 0, /* 0:1.000000 */
	     /* lpb=14 */ 0, /* 0:1.000000 */
	     /* lpb=15 */ 0, /* 0:1.000000 */
	     /* lpb=16 */ 0, /* 0:1.000000 */
	     /* lpb=17 */ 1, /* 0:0.980000, 1:1.000000 */
	     /* lpb=18 */ 1, /* 1:0.990000 */
	     /* lpb=19 */ 2, /* 1:0.980000, 2:1.000000 */
	     /* lpb=20 */ 2, /* 2:0.990000 */
	     /* lpb=21 */ 2, /* 2:0.990000 */
	     /* lpb=22 */ 4, /* 3:0.980000, 4:1.000000 */
	     /* lpb=23 */ 4, /* 4:0.990000 */
	     /* lpb=24 */ 5, /* 4:0.950000, 5:0.990000 */
	     /* lpb=25 */ 8, /* 7:0.980000, 8:0.990000 */
	     /* lpb=26 */ 12, /* 11:0.980000, 12:1.000000 */
	     /* lpb=27 */ 14, /* 13:0.980000, 14:1.000000 */
	     /* lpb=28 */ 14, /* 14:1.000000 */
	     /* lpb=29 */ 14, /* 14:0.990000 */
	     /* lpb=30 */ 17, /* 16:0.970000, 17:0.990000 */
	     /* lpb=31 */ 25, /* 24:0.980000, 25:1.000000 */
	     /* lpb=32 */ 25, /* 25:0.990000 */
	     /* lpb=33 */ 25, /* 25:0.990000 */
	     /* lpb=34 */ 32, /* 31:0.980000, 32:0.990000 */
	     /* lpb=35 */ 36, /* 35:0.980000, 36:0.990000 */
	     /* lpb=36 */ 39, /* 38:0.980000, 39:0.990000 */
	     /* lpb=37 */ 43, /* 42:0.980000, 43:0.990000 */
	     /* lpb=38 */ 49, /* 48:0.980000, 49:0.990000 */
	     /* lpb=39 */ 52, /* 51:0.970000, 52:0.990000 */
	     /* lpb=40 */ 56, /* 55:0.980000, 56:0.990000 */
	     /* lpb=41 */ 60, /* 59:0.980000, 60:0.990000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}
#endif
/*}}}*/

struct strategy_file_parser {/*{{{*/
    typedef std::array<unsigned int, 2> key_type;;
    typedef std::vector<facul_method::parameters_with_side> value_type;
private:
    class regexp_define_t {/*{{{*/
        regex_t re;
        public:
        typedef std::string T;
        regexp_define_t() {
            const char * re_txt =
                "^[[:space:]]*"
                "define[[:space:]]+"
                "([[:alnum:]_]+)";
            regcomp (&re, re_txt, REG_ICASE|REG_EXTENDED);
        }
        ~regexp_define_t() {
            regfree(&re);
        }
        bool operator()(T & name, const char * & str) const
        {
            constexpr const int nmatch = 2;
            regmatch_t p[nmatch];
            if (regexec (&re, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            name = std::string(str + p[1].rm_so, str + p[1].rm_eo);
            str += p[0].rm_eo;
            return true;
        }
    };/*}}}*/
    regexp_define_t regexp_define;

    class regexp_use_t {/*{{{*/
        regex_t re;
        public:
        typedef std::string T;
        regexp_use_t() {
            const char * re_txt =
                "^[[:space:]]*"
                "use[[:space:]]+"
                "([[:alnum:]_]+)";
            regcomp (&re, re_txt, REG_ICASE|REG_EXTENDED);
        }
        ~regexp_use_t() {
            regfree(&re);
        }
        bool operator()(T & name, const char * & str) const
        {
            constexpr const int nmatch = 2;
            regmatch_t p[nmatch];
            if (regexec (&re, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            name = std::string(str + p[1].rm_so, str + p[1].rm_eo);
            str += p[0].rm_eo;
            return true;
        }
    };/*}}}*/
    regexp_use_t regexp_use;

    class regexp_index_t {/*{{{*/
        regex_t re;
        public:
        typedef key_type T;
        regexp_index_t() {
            const char * re_txt =
                "^\\[?[[:space:]]*"
                "r0=([[:digit:]]+)"
                ",[[:space:]]*"
                "r1=([[:digit:]]+)"
                "\\]?[[:space:]]*";
            regcomp (&re, re_txt, REG_ICASE|REG_EXTENDED);
        }
        ~regexp_index_t() {
            regfree(&re);
        }
        bool operator()(T & index, const char * & str) const
        {
            constexpr const int nmatch = 3;
            regmatch_t p[nmatch];
            if (regexec (&re, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            index[0] = std::stoi(std::string(str + p[1].rm_so, str + p[1].rm_eo));
            index[1] = std::stoi(std::string(str + p[2].rm_so, str + p[2].rm_eo));
            str += p[0].rm_eo;
            return true;
        }
    };/*}}}*/
    regexp_index_t regexp_index;

    class regexp_timing_comment_t {/*{{{*/
        regex_t re;
        public:
        struct T {
            double p;
            double t;
        };
        regexp_timing_comment_t() {
            const char *re_txt =
                "^:[[:space:]]*"
                "\\([[:space:]]*"
                "p[[:space:]]*=[[:space:]]*([\\.[:digit:]]+)"
                "[[:space:]]*,[[:space:]]*"
                "t[[:space:]]*=[[:space:]]*([\\.[:digit:]]+)"
                "[[:space:]]*\\)";
            regcomp (&re, re_txt, REG_ICASE|REG_EXTENDED);
        }
        ~regexp_timing_comment_t() {
            regfree(&re);
        }
        bool operator()(T & comment, const char * & str) const
        {
            constexpr const int nmatch = 3;
            regmatch_t p[nmatch];
            if (regexec (&re, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            comment.p = std::stod(std::string(str + p[1].rm_so, str + p[1].rm_eo));
            comment.t = std::stod(std::string(str + p[2].rm_so, str + p[2].rm_eo));
            str += p[0].rm_eo;
            return true;
        }
    };/*}}}*/
    regexp_timing_comment_t regexp_timing_comment;

    class regexp_fm_t {/*{{{*/
        regex_t preg_fm;
        public:
        typedef facul_method::parameters_with_side T;
        regexp_fm_t() {
            // regular expression for the strategy
            const char *str_preg_fm =
                "^\\[?[[:space:]]*"
                "S([[:alnum:]]+):[[:space:]]*"  /* side, like "S0: " */
                "([[:alnum:]_-]+)"               /* method, like "PP1-65" or "ECM-M12" */
                ",[[:space:]]*([[:digit:]]+)"   /* B1, an integer */
                ",[[:space:]]*([[:digit:]]+)"   /* B2, an integer */
                "[[:space:]]*"
                "\\]?[[:space:]]*";
            regcomp (&preg_fm, str_preg_fm, REG_ICASE|REG_EXTENDED);
        }
        ~regexp_fm_t() {
            regfree(&preg_fm);
        }
        bool operator()(T & res, const char * & str) const
        {
            constexpr const int nmatch = 5;
            regmatch_t p[nmatch];
            if (regexec (&preg_fm, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            if (p[0].rm_so == p[0].rm_eo)
                return false;
            res.side = std::stoi(std::string(str + p[1].rm_so, str + p[1].rm_eo));
            res.B1 = std::stoi(std::string(str + p[3].rm_so, str + p[3].rm_eo));
            res.B2 = std::stoi(std::string(str + p[4].rm_so, str + p[4].rm_eo));
            std::string mtoken(str + p[2].rm_so, str + p[2].rm_eo);
            if (mtoken == "PM1")
                res.method = PM1_METHOD;
            else if (mtoken == "PP1-27")
                res.method = PP1_27_METHOD;
            else if (mtoken == "PP1-65")
                res.method = PP1_65_METHOD;
            else {
                res.method = EC_METHOD;
                if (mtoken == "ECM-B12")
                    res.parameterization = BRENT12;
                else if (mtoken == "ECM-M12")
                    res.parameterization = MONTY12;
                else if (mtoken == "ECM-M16")
                    res.parameterization = MONTY16;
                else if (mtoken == "ECM-TM12")
                    res.parameterization = MONTYTWED12;
                else if (mtoken == "ECM-TM16")
                    res.parameterization = MONTYTWED16;
                else {
                    fprintf (stderr,
                            "error : unknown method '%s' in strategy file\n",
                            mtoken.c_str());
                    return false;
                }
            }
            str += p[0].rm_eo;
            return true;
        }
    };/*}}}*/
    regexp_fm_t regexp_fm;

    public:

    /* This only returns the vector of descriptions. The methods are not
     * instantiated yet
     */
    facul_strategies::strategy_file
        operator()(std::array<unsigned int, 2> const & mfb, FILE * file)
    {
        verbose_output_print(0, 2, "# Read the cofactorization strategy file\n");
        // first, read linearly.
        std::vector<std::pair<key_type, value_type>> pre_parse;
        std::map<std::string, value_type> macros;

        value_type * current = nullptr;
        
        fseek (file, 0, SEEK_SET);
        for(char line[10000]; fgets (line, sizeof(line), file) != NULL ; )
        {
            const char * str = line;

            for( ; *str != '\0' ; ) {
                if (isspace(*str)) { str++; continue; }

                // This removes comments.
                if (*str == '#') break;

                key_type index_st;
                value_type::value_type fm;
                std::string macro;

                if (regexp_index(index_st, str)) {
                    regexp_timing_comment_t::T p_and_t;
                    regexp_timing_comment(p_and_t, str);

                    pre_parse.emplace_back(index_st, value_type());
                    current = &pre_parse.back().second;

                    continue;
                } else if (regexp_define(macro, str)) {
                    auto it = macros.emplace(macro, value_type());
                    if (!it.second) {
                        fprintf(stderr, "# macro %s redefined in strategies file\n", macro.c_str());
                        exit(EXIT_FAILURE);
                    }
                    current = &it.first->second;
                } else if (regexp_use(macro, str)) {
                    if (current == nullptr) {
                        fprintf(stderr, "# dangling methods in strategies file\n");
                        exit(EXIT_FAILURE);
                    }
                    auto it = macros.find(macro);
                    if (it == macros.end()) {
                        fprintf(stderr, "# macro %s unknown in strategies file\n", macro.c_str());
                        exit(EXIT_FAILURE);
                    }
                    auto const & M = it->second;
                    if (current == &M) {
                        fprintf(stderr, "# circular use of macro %s in strategies file\n", macro.c_str());
                        exit(EXIT_FAILURE);
                    }
                    current->insert(current->end(), M.begin(), M.end());
                } else if (regexp_fm(fm, str)) {
                    if (current == nullptr) {
                        fprintf(stderr, "# dangling methods in strategies file\n");
                        exit(EXIT_FAILURE);
                    }
                    current->push_back(fm);
                } else {
                    fprintf(stderr, "# parse error in strategies file\n");
                    fprintf(stderr, "# cannot parse: %s\n", str);
                    exit(EXIT_FAILURE);
                }

            }
        }

        std::map<key_type, value_type> parsed_file;

        for(auto & c : pre_parse) {
            key_type index_st = c.first;

            if (index_st[0] > mfb[0])
                continue;
            if (index_st[1] > mfb[1])
                continue;
            if (c.second.empty())
                continue;

            std::map<unsigned int, std::array<unsigned long, 2>>
                param_sequence {
                    { BRENT12, { 0, 0 } },
                    { MONTY12, { 0, 0 } },
                    { MONTY16, { 0, 0 } },
                    { MONTYTWED12, { 0, 0 } },
                    /* we have no parameter_from_sequence for MONTY16
                     * anyway. FIXME */
                    { MONTYTWED16, { 0, 0 } },
                };

            for(auto & fm : c.second) {
                if (fm.method == EC_METHOD) {
                    fm.parameter = ec_valid_parameter_from_sequence(fm.parameterization, param_sequence[fm.parameterization][fm.side]++);

                    /* it seems that we always have this on */
                    fm.extra_primes = 1;
                }
            }

            parsed_file.insert(c);
        }

        return parsed_file;
    }
};/*}}}*/

/* TODO: refactor this closer to facul_method */
void facul_strategies::print(FILE * file) const/*{{{*/
{
    if (file == NULL)
        return;
    // print info lpb ...
    fprintf (file,
            "(lpb = [%u,%u], as...=[%lf, %lf], BBB = [%lf, %lf])\n",
            lpb[0], lpb[1],
            BB[0], BB[1],
            BBB[0], BBB[1]);
    fprintf (file, "mfb = [%d, %d]\n", mfb[0], mfb[1]);
    for (unsigned int r = 0; r <= mfb[0]; r++) {
        for (unsigned int a = 0; a <= mfb[1]; a++) {
            fprintf (file, "r0=%u, r1=%u", r, a);

            for(facul_method_side const & ms : (*this)(r, a)) {
                int side = ms.side;
                facul_method const & fm = * ms.method;
                switch(fm.method) {
                    case PM1_METHOD:
                        fprintf (file, " S%d PM1,%d,%d", side,
                                ((pm1_plan_t*) fm.plan)->B1,
                                ((pm1_plan_t*) fm.plan)->stage2.B2);
                        break;
                    case PP1_27_METHOD:
                        fprintf (file, " S%d PP1-27,%d,%d", side,
                                ((pp1_plan_t*) fm.plan)->B1,
                                ((pp1_plan_t*) fm.plan)->stage2.B2);
                        break;
                    case PP1_65_METHOD:
                        fprintf (file, " S%d PP1-65,%d,%d", side,
                                ((pp1_plan_t*) fm.plan)->B1,
                                ((pp1_plan_t*) fm.plan)->stage2.B2);
                        break;
                    case EC_METHOD:
                        switch(((ecm_plan_t*) fm.plan)->parameterization) {
                            case BRENT12:
                                fprintf (file, " S%d ECM-B12,%d,%d", side,
                                        ((ecm_plan_t*) fm.plan)->B1,
                                        ((ecm_plan_t*) fm.plan)->stage2.B2);
                                break;
                            case MONTY12:
                                fprintf (file, " S%d ECM-M12,%d,%d", side,
                                        ((ecm_plan_t*) fm.plan)->B1,
                                        ((ecm_plan_t*) fm.plan)->stage2.B2);
                                break;
                            case MONTY16:
                                fprintf (file, " S%d ECM-M16,%d,%d", side,
                                        ((ecm_plan_t*) fm.plan)->B1,
                                        ((ecm_plan_t*) fm.plan)->stage2.B2);
                                break;
                            case MONTYTWED12:
                                fprintf (file, " S%d ECM-TM12,%d,%d", side,
                                        ((ecm_plan_t*) fm.plan)->B1,
                                        ((ecm_plan_t*) fm.plan)->stage2.B2);
                                break;
                            case MONTYTWED16:
                                fprintf (file, " S%d ECM-TM16,%d,%d", side,
                                        ((ecm_plan_t*) fm.plan)->B1,
                                        ((ecm_plan_t*) fm.plan)->stage2.B2);
                                break;
                        }
                        break;
                }
                fprintf (file, "\n");
            }
        }
    }
}/*}}}*/

facul_strategies_base::facul_strategies_base (
        std::array<unsigned long, 2> const & lim,
        std::array<unsigned int, 2> const & lpb,
        std::array<unsigned int, 2> const & mfb,
        bool perfectly_sieved)
    : B(lim), lpb(lpb), mfb(mfb)
{
    auto pBB = BB.begin();
    auto pBBB = BBB.begin();
    for(double b: B) {
        *pBB++ = perfectly_sieved ? b*b : 0;
        *pBBB++ = b*b*b;
    }
}

void facul_strategies::precompute_method(facul_method::parameters const & mp, int verbose)
{
    /* How do I emplace with a guarantee of not calling the ctor
     * needlessly? map::emplace is actually allowed to construct
     * unconditionally.
     */
    auto it = precomputed_methods.emplace(mp, facul_method());
    if (it.second)
        it.first->second = facul_method(mp, verbose);
}

/* Create our strategy book from a file. */
facul_strategies::facul_strategies(
        std::array<unsigned long, 2> const & lim,
        std::array<unsigned int, 2> const & lpb,
        std::array<unsigned int, 2> const & mfb,
        bool perfectly_sieved,
        FILE * file,
        const int verbose)
    : facul_strategies(lim, lpb, mfb, perfectly_sieved, strategy_file_parser()(mfb, file), verbose)
{}
 
/* Create our strategy book from an in-memory version of a file. */
facul_strategies::facul_strategies(
        std::array<unsigned long, 2> const & lim,
        std::array<unsigned int, 2> const & lpb,
        std::array<unsigned int, 2> const & mfb,
        bool perfectly_sieved,
        facul_strategies::strategy_file const & parsed_file,
        const int verbose)
    : facul_strategies_base(lim, lpb, mfb, perfectly_sieved)
{
    /* Add all methods to the cache */
    for(auto const & v : parsed_file) {
        auto const & chain_parameters = v.second;
        for(facul_method::parameters const & mp: chain_parameters) {
            precompute_method(mp, verbose);
        }
    }

    /* Add all method chains to the cache */
    for(auto const & v : parsed_file) {
        auto const & chain_parameters = v.second;
        auto it = precomputed_strategies.find(chain_parameters);
        if (it != precomputed_strategies.end())
            continue;
        auto & chain = precomputed_strategies[chain_parameters];
        for(auto const & ms: chain_parameters)
            chain.emplace_back(&precomputed_methods[ms], ms.side);
        facul_method_side::fix_is_last(chain);
    }

    direct_access.assign((1 + mfb[0]) * (1 + mfb[1]), nullptr);

    /* fill the direct_access array */
    for(auto const & v : parsed_file) {
        auto const & i = v.first;
        auto const & chain_parameters = v.second;
        auto const * w = &precomputed_strategies.at(chain_parameters);
        direct_access_get(i[0], i[1]) = w;
    }
}

std::vector<facul_method_side> const & facul_strategies::operator()(unsigned int r, unsigned int a) const
{
    if (direct_access.empty()) {
        return uniform_strategy[r < a];
    } else {
        static std::vector<facul_method_side> placeholder;
        auto it = direct_access_get(r, a);
        if (it != nullptr)
            return *it;
        else
            return placeholder;
    }
}

std::vector<facul_method_side> const * & facul_strategies::direct_access_get(unsigned int r, unsigned int a)
{
    ASSERT_ALWAYS(r <= mfb[0]);
    ASSERT_ALWAYS(a <= mfb[1]);
    return direct_access[r + a * (1 + mfb[0])];
}

std::vector<facul_method_side> const * const & facul_strategies::direct_access_get(unsigned int r, unsigned int a) const
{
    ASSERT_ALWAYS(r <= mfb[0]);
    ASSERT_ALWAYS(a <= mfb[1]);
    return direct_access[r + a * (1 + mfb[0])];
}

/*
  Make a simple strategy for factoring. We start with
  P-1 and P+1 (with x0=2/7), then an ECM curve with low bounds, then
  a bunch of ECM curves with larger bounds. How many methods to do in
  total is controlled by the n parameter: P-1, P+1 and the first ECM
  curve (with small bounds) are always done, then n ECM curves (with
  larger bounds).
  This function is used when you don't give a strategy file.
*/
std::vector<facul_method::parameters>
facul_strategy_oneside::default_strategy (int n)
{
    ASSERT_ALWAYS(n >= 0);

    std::vector<facul_method::parameters> chain;

#if 0
    /* This is relevant only for very weird experiments where we have small
     * factors that stick together.  */
    chain.emplace_back({ PM1_METHOD, 30, 100 });
#endif

    chain.emplace_back(PM1_METHOD, 315, 2205);
    chain.emplace_back(PP1_27_METHOD, 525, 3255);
#ifdef USE_LEGACY_DEFAULT_STRATEGY
    chain.emplace_back(EC_METHOD, 105, 3255, MONTY12, 2, 1);
#else
    chain.emplace_back(EC_METHOD, 105, 3255, MONTYTWED12, 1, 1);
#endif

    if (n--) {
#ifdef USE_LEGACY_DEFAULT_STRATEGY
        chain.emplace_back(EC_METHOD, 315, 5355, BRENT12, 11, 1);
#else
        chain.emplace_back(EC_METHOD, 315, 5355, MONTYTWED12, 2, 1);
#endif
    }

    /* heuristic strategy where B1 is increased by c*sqrt(B1) at each curve
     *
     * Note that only some addition chains were computed in
     * sieve/ecm/bytecode_mishmash_B1_data.h ; those have fewer
     * arithmetic operations than the ones we compute automatically.
     *
     * For this reason, if the sequence of attained values for B1 is
     * changed, the performance mught degrade a little bit.
     */
    double B1 = 105;
    for (int i = 3 ; n-- > 0 ; i++) {
        /* If the sequence of B1 values is modified, it may be a good thing to
         * regenerate bytecode_mishmash_B1_data.h to add precomputed chains for
         * the new B1 values.
         */
        B1 += sqrt (B1);

        /* The factor 50 was determined experimentally with testbench, to find
           factors of 40 bits:
           testbench -p -cof 1208925819614629174706189 -strat 549755813888 549755913888
           This finds 1908 factors (out of 3671 input numbers) with n=24 curves
           and 3.66s.
           With B2=17*B1, and 29 curves, we find 1898 factors in 4.07s.
           With B2=100*B1, and 21 curves we find 1856 factors in 3.76s.
           Thus 50 seems close to optimal.
Warning: changing the value of B2 might break the sieving tests.
*/
        /* we round B2 to (2k+1)*105, thus k is the integer nearest to
           B2/210-0.5 */
        unsigned long B2 = (2 * (unsigned int) ((50.0 * B1) / 210.0) + 1) * 105;

#ifdef USE_LEGACY_DEFAULT_STRATEGY
        /* FIXME: are MONTY12 and MONTYTWED12 swapped, perhaps ? */
        chain.emplace_back(EC_METHOD, B1, B2, MONTYTWED12, i, 1);
#else
        chain.emplace_back(EC_METHOD, B1, B2, MONTY12, i, 1);
#endif
    }

#ifdef USE_MPQS
    chain.emplace_back({ MPQS_METHOD });
#endif

    return chain;
}

facul_strategy_oneside::facul_strategy_oneside (unsigned long fbb, unsigned int lpb, unsigned int mfb, std::vector<facul_method::parameters> const & mps, int verbose)
    : B(fbb),
      lpb(lpb),
      BB((double) B * (double) B),
      BBB((double) B * (double) B * (double) B),
      mfb(mfb)
{
    for(auto const & mp : mps)
        methods.emplace_back(mp, verbose);
}

facul_strategy_oneside::facul_strategy_oneside (unsigned long fbb, unsigned int lpb, unsigned int mfb, int n, int verbose)
    : facul_strategy_oneside(fbb, lpb, mfb, default_strategy(n < 0 ? nb_curves_with_fbb (fbb, lpb, mfb) : n), verbose)
{}


/*
 * Create our strategy book with our default strategy.
 */
facul_strategies::facul_strategies (
        std::array<unsigned long, 2> const & lim,
        std::array<unsigned int, 2> const & lpb,
        std::array<unsigned int, 2> const & mfb,
        std::array<int, 2> ncurves,
        bool perfectly_sieved,
	const int verbose)
    : facul_strategies_base(lim, lpb, mfb, perfectly_sieved)
{
    int max_ncurves = -1;
    for(int side = 0 ; side < 2 ; side++) {
        if (ncurves[side] < 0)
            ncurves[side] = nb_curves_with_fbb (B[side], lpb[side], mfb[side]);
        if (ncurves[side] > max_ncurves)
            max_ncurves = ncurves[side];
    }

    verbose_output_print(0, 2, "# Using default strategy for the cofactorization: ncurves0=%d ncurves1=%d\n", ncurves[0], ncurves[1]);

    /* prepare the chain of methods that we want to use in order to
     * factor a number, irrespective of its side.
     */
    auto chain_parameters = facul_strategy_oneside::default_strategy(max_ncurves);

    /* Add all methods to the cache */
    for(facul_method::parameters const & mp: chain_parameters)
        precompute_method(mp, verbose);

    /* We now need to truncate the list of strategies according to
     * ncurves[0] and ncurves[1]
     *
     * NOTE: This is incompatible with USE_MPQS
     */
    for(int first = 0 ; first < 2 ; first++) {
        /* first == 0 means that r >= a: the rational side is largest.
         * Try to factor it first.
         *
         * Note that facul_strategies::operator() returns
         * uniform_strategy[r < a]
         */
        std::vector<facul_method_side> & u(uniform_strategy[first]);
        for (int z = 0; z < 2; z++) {
            int side = first ^ z;
            int n = ncurves[side] + (chain_parameters.size() - max_ncurves);
            for(facul_method::parameters const & mp: chain_parameters) {
                if (!n--)
                    break;
                u.emplace_back(&precomputed_methods[mp], side);
            }
        }
    }
}
