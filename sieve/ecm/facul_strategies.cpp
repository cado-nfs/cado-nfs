#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <regex.h>

#include "fmt/format.h"
#include "fmt/ranges.h"

#include "facul_strategies.hpp"
#include "facul_ecm.h"
#include "facul_method.hpp"
#include "pm1.h"
#include "pp1.h"
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
  int const T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
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

/* TODO: move elsewhere. */
static const char * parameterization_name(ec_parameterization_t p)
{
    switch(p) {
        case BRENT12:     return "ECM-B12";
        case MONTY12:     return "ECM-M12";
        case MONTY16:     return "ECM-M16";
        case MONTYTWED12: return "ECM-TM12";
        case MONTYTWED16: return "ECM-TM16";
    }
    ASSERT_ALWAYS(0);
    return nullptr;
}

struct strategy_file_parser {/*{{{*/
    using key_type = std::vector<unsigned int>;
    using value_type = std::vector<facul_method::parameters_with_side>;
private:
    class regexp_define_t {/*{{{*/
        regex_t re;
        public:
        using T = std::string;
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
        using T = std::string;
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
        std::regex re;

        static std::regex build_regex(size_t n) {
            std::ostringstream s;
            for (unsigned int i = 0; i < n; ++i) {
                s << (i == 0 ? "^\\[?" : ",")
                  << "[[:space:]]*r" << i << "=([[:digit:]]+)[[:space:]]*";
            }
            s << "]?[[:space:]]*";
            return std::regex{s.str(), std::regex::extended|std::regex::icase};
        }

        public:
        using T = std::vector<unsigned int>;
        regexp_index_t(size_t n) : re{build_regex(n)} {
        }
        bool operator()(T & index, const char * & str) const {
            std::cmatch m;
            if (!regex_search(str, m, re)) {
                return false;
            } else {
                for (size_t i = 1; i < m.size(); ++i) {
                    index[i-1] = std::stoi(m[i].str());
                }
                str += m.length();
                return true;
            }
        }

    };/*}}}*/

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
        using T = facul_method::parameters_with_side;
        regexp_fm_t() {
            // regular expression for the strategy
            const char *str_preg_fm =
                "^\\[?[[:space:]]*"
                "S([[:alnum:]]+):[[:space:]]*"  /* side, like "S0: " */
                "([[:alnum:]_-]+)"               /* method, like "PP1-65" or "ECM-M12" */
                "(\\(([[:digit:]]+)\\))?" /* optional parameter in parentheses */
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
            constexpr const int nmatch = 7;
            regmatch_t p[nmatch];
            if (regexec (&preg_fm, str, nmatch, p, 0) == REG_NOMATCH)
                return false;
            if (p[0].rm_so == p[0].rm_eo)
                return false;
            res.side = std::stoi(std::string(str + p[1].rm_so, str + p[1].rm_eo));
            res.parameter = ULONG_MAX;
            if (p[4].rm_so != p[4].rm_eo)
                res.parameter = std::stoi(std::string(str + p[4].rm_so, str + p[4].rm_eo));
            res.B1 = std::stoi(std::string(str + p[5].rm_so, str + p[5].rm_eo));
            res.B2 = std::stoi(std::string(str + p[6].rm_so, str + p[6].rm_eo));
            std::string const mtoken(str + p[2].rm_so, str + p[2].rm_eo);
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

    struct error : public std::runtime_error {
        explicit error(std::string const & s) : std::runtime_error(s) {}
    };

    /* This only returns the vector of descriptions. The methods are not
     * instantiated yet
     */
    facul_strategies::strategy_file
        operator()(std::vector<unsigned int> const & mfb, FILE * file);
};/*}}}*/

class parameter_sequence_tracker {/*{{{*/
    std::map<unsigned int, std::vector<std::pair<bool, unsigned long>>> seq;
    public:
    parameter_sequence_tracker() {
        decltype(seq)::mapped_type::value_type const z { true, 0 };
        decltype(seq)::mapped_type const Z { z , z };
        seq[BRENT12] = Z;
        seq[MONTY12] = Z;
        seq[MONTY16] = Z;
        seq[MONTYTWED12] = Z;
        /* we have no parameter_from_sequence for MONTY16
         * anyway. FIXME */
        seq[MONTYTWED16] = Z;
    }

    /* This returns true if p is actually the next "normal" parameter in
     * our sequence, given our default choices, for this parameterization
     * and on this side.
     * If true, then our current record of which parameters have been
     * used is updated.
     * If false, we record the fact that we have broken the sequence. All
     * further calls to this method, on this side and for this
     * parameterization, will return false.
     */
    bool follows_sequence(int side, ec_parameterization_t para, unsigned long p)
    {
        if (!seq[para][side].first) return false;
        bool const t = (p == ec_valid_parameter_from_sequence(para, seq[para][side].second));
        if (t) seq[para][side].second++;
        return t;
    }
    unsigned long next_in_sequence(int side, ec_parameterization_t para)
    {
        /* If the sequence has been broken, we can no longer give a
         * "next" value. Of course we could improve this by recording the
         * out-of-sequence parameters that have been chosen, and exclude
         * them from the default choices.
         */
        if (!seq[para][side].first)
            throw strategy_file_parser::error(fmt::format("cannot use a parameter from the default sequence for parameterization {} after a non-default one was set", parameterization_name(para)));
        return ec_valid_parameter_from_sequence(para, seq[para][side].second++);
    }
    void break_sequence(int side, ec_parameterization_t para)
    {
        seq[para][side].first = false;
    }

    static void fill_default_parameters(std::vector<facul_method::parameters_with_side> & v) {
        parameter_sequence_tracker tracker;

        for(auto & fm : v) {
            if (fm.method == EC_METHOD) {
                if (fm.parameter == ULONG_MAX) {
                    fm.parameter = tracker.next_in_sequence(fm.side, fm.parameterization);
                } else {
                    tracker.break_sequence(fm.side, fm.parameterization);
                }
                /* it seems that we always have this on */
                fm.extra_primes = 1;
            }
        }
    }
};/*}}}*/

facul_strategies::strategy_file
strategy_file_parser::operator()(std::vector<unsigned int> const & mfb, FILE * file)
{
    verbose_fmt_print(0, 2, "# Read the cofactorization strategy file\n");
    // first, read linearly.
    std::vector<std::pair<key_type, value_type>> pre_parse;
    std::map<std::string, value_type> macros;

    regexp_index_t regexp_index{mfb.size()};

    value_type * current = nullptr;

    fseek (file, 0, SEEK_SET);
    int lnum = 1;
    try {
        for(char line[10000]; fgets (line, sizeof(line), file) != nullptr ; lnum++)
        {
            const char * str = line;

            for( ; *str != '\0' ; ) {
                if (isspace(*str)) { str++; continue; }

                // This removes comments.
                if (*str == '#') break;

                key_type index_st(mfb.size());
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
                    if (!it.second)
                        throw error(fmt::format("macro {} redefined", macro));
                    current = &it.first->second;
                } else if (regexp_use(macro, str)) {
                    if (current == nullptr)
                        throw error("dangling macro use");
                    auto it = macros.find(macro);
                    if (it == macros.end())
                        throw error(fmt::format("macro {} unknown", macro));
                    auto const & M = it->second;
                    if (current == &M) {
                        throw error(fmt::format("circular use of macro {} file", macro));
                    }
                    current->insert(current->end(), M.begin(), M.end());
                } else if (regexp_fm(fm, str)) {
                    if (current == nullptr)
                        throw error("dangling methods");
                    current->push_back(fm);
                } else {
                    throw error(fmt::format("cannot parse {}", str));
                }
            }
        }
    } catch(error const & e) {
        throw std::runtime_error(fmt::format(
                        "Parse error on line {} of strategies file: {}",
                        lnum, e.what()));
    }

    std::map<key_type, value_type> parsed_file;

    for(auto & c : pre_parse) {
        key_type index_st = c.first;

        /* std::ranges::equal use the std:less_equal for comparisons so what it is
         * really doing is checking if index_st[i] <= mfb[i] for all i
         */
        if(!std::ranges::equal(index_st, mfb,
                       std::less_equal<unsigned int>{}))
            continue; /* it exists i such that index_st[i] > mfb[i] */
        if (c.second.empty())
            continue;

        /* Note that we're now setting all parameters for the (r0,r1)
         * pairs that are given in the file. This might, in effect, cause
         * repeated checking of the same things over and over again.
         */

        try {
            parameter_sequence_tracker::fill_default_parameters(c.second);
        } catch(error const & e) {
            throw std::runtime_error(fmt::format(
                        "Parse error in strategies file while setting parameters for r={}: {}",
                        index_st, e.what()));
        }
        parsed_file.insert(c);
    }

    return parsed_file;
}

static void fprint_one_chain(FILE * file, std::vector<facul_method_side> const & v)
{
    parameter_sequence_tracker tracker;

    /* TODO: refactor at least some of this closer to facul_method */
    for(facul_method_side const & ms : v) {
        int const side = ms.side;
        facul_method const & fm = * ms.method;
        switch(fm.method) {
            case PM1_METHOD:
                fprintf (file, " S%d: PM1,%d,%d\n", side,
                        ((pm1_plan_t*) fm.plan)->B1,
                        ((pm1_plan_t*) fm.plan)->stage2.B2);
                break;
            case PP1_27_METHOD:
                fprintf (file, " S%d: PP1-27,%d,%d\n", side,
                        ((pp1_plan_t*) fm.plan)->B1,
                        ((pp1_plan_t*) fm.plan)->stage2.B2);
                break;
            case PP1_65_METHOD:
                fprintf (file, " S%d: PP1-65,%d,%d\n", side,
                        ((pp1_plan_t*) fm.plan)->B1,
                        ((pp1_plan_t*) fm.plan)->stage2.B2);
                break;
            case EC_METHOD:
                {
                    auto const * e = (ecm_plan_t const *) fm.plan;
                    if (tracker.follows_sequence(ms.side, e->parameterization, e->parameter)) {
                        fprintf (file, " S%d: %s,%d,%d\n", side,
                                parameterization_name(e->parameterization),
                                e->B1,
                                e->stage2.B2);
                    } else {
                        fprintf (file, " S%d: %s(%lu),%d,%d\n", side,
                                parameterization_name(e->parameterization),
                                e->parameter,
                                e->B1,
                                e->stage2.B2);
                    }
                }
                break;
            default:
                ASSERT_ALWAYS(0);
        }
    }
}

template<class UnaryOp>
void facul_strategies::for_each_sizes_combination(UnaryOp f) const
{
    std::vector<unsigned int> v(mfb.size());
    for_each_sizes_combination_inner(0, v, f);
}

template<class UnaryOp>
void facul_strategies::for_each_sizes_combination_inner(unsigned int i, std::vector<unsigned int> & v, UnaryOp f) const
{
    if (i >= mfb.size()) {
        f(v);
    } else {
        for(size_t j = 0; j < mfb[i]; ++j) {
            v[i] = j;
            for_each_sizes_combination_inner(i+1, v, f);
        }
    }
}

void facul_strategies::print(FILE * file) const/*{{{*/
{
    if (file == nullptr)
        return;
    // print info lpb ...
    fmt::print(file, "# (lpb = {}, as...={}, BBB = {::.0Lf})\n", lpb, BB, BBB);
    fmt::print(file, "# mfb = {}\n", mfb);

    // operator() always returns a reference to things that are either in
    // the structure's precomputed_strategies, uniform_strategies, or
    // even static placeholder things in the code. We can be confident
    // that successive calls will return the same pointers.

    std::map<std::vector<facul_method_side> const *, unsigned int> popularity;

    for_each_sizes_combination(
        [this, &popularity] (std::vector<unsigned int> & sizes) {
            popularity[&(*this)(sizes)]++;
        });

    for_each_sizes_combination(
        [this, &popularity, &file] (std::vector<unsigned int> & sizes) {
            std::vector<facul_method_side> const & v = (*this)(sizes);
            unsigned int const p = popularity[&v];
            ASSERT_ALWAYS(p > 0);       // see above.
            if (p >= 2 && v.size() >= 2) {
                fmt::print(file, "define same_as_{}\n", fmt::join(sizes, "_"));
                fprint_one_chain(file, v);
            }
        });

    for_each_sizes_combination(
        [this, &popularity, &file] (std::vector<unsigned int> & sizes) {
            std::vector<facul_method_side> const & v = (*this)(sizes);
            if (v.empty()) return;
            for(size_t i = 0; i < sizes.size(); ++i) {
                fmt::print(file, "r{}={}{}", i, sizes[i],
                                             i+1 < sizes.size() ? ',' : '\n');
            }
            unsigned int const p = popularity[&v];
            ASSERT_ALWAYS(p > 0);       // see above.
            if (p >= 2 && v.size() >= 2) {
                fmt::print(file, "  use same_as_{}\n", fmt::join(sizes, "_"));
            } else {
                fprint_one_chain(file, v);
            }
        });
}/*}}}*/

facul_strategies_base::facul_strategies_base (
        std::vector<unsigned long> const & lim,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & mfb,
        bool perfectly_sieved)
    : B(lim), lpb(lpb), BB(lim.size()), BBB(lim.size()), mfb(mfb)
{
    ASSERT_ALWAYS(lpb.size() == lim.size());
    ASSERT_ALWAYS(mfb.size() == lim.size());
    auto pBB = BB.begin();
    auto pBBB = BBB.begin();
    for(double const b: B) {
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
        std::vector<unsigned long> const & lim,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & mfb,
        bool perfectly_sieved,
        FILE * file,
        const int verbose)
    : facul_strategies(lim, lpb, mfb, perfectly_sieved, strategy_file_parser()(mfb, file), verbose)
{
}
 
/* Create our strategy book from an in-memory version of a file. */
facul_strategies::facul_strategies(
        std::vector<unsigned long> const & lim,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & mfb,
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

    unsigned int M = 1U;
    for (unsigned int m: mfb) {
        M *= (m+1U);
    }
    direct_access.assign(M, nullptr);

    /* fill the direct_access array */
    for(auto const & v : parsed_file) {
        auto const & i = v.first;
        auto const & chain_parameters = v.second;
        auto const * w = &precomputed_strategies.at(chain_parameters);
        direct_access_get(i) = w;
    }
}

std::vector<facul_method_side> const & facul_strategies::operator()(std::vector<unsigned int> const & v) const
{
    if (direct_access.empty()) {
        if (v.size() == 1) /* one side */
            return uniform_strategy[0];
        else if (v.size() == 2) /* two sides */
            return uniform_strategy[v[0] < v[1]];
        else { /* more than two sides */
            std::vector<unsigned int> index(v.size());
            std::iota(index.begin(), index.end(), 0);
            std::ranges::sort(index,
                      [&](const unsigned int& i, const unsigned int& j) {
                          return (v[i] < v[j]);
                      });
            /* index[0] contains the id of the largest side, index[1] the id of
             * the second largest, ...
             * We now need the lexicographic index of this permutation.
             * Note: cost is quadratic in the number of sides.
             * Found the formula at https://stackoverflow.com/questions/12146910/finding-the-lexicographic-index-of-a-permutation-of-a-given-array
             * Adapted it "Horner-style" to avoid the factorials.
             */
            size_t idx = 0;
            for(size_t i = 0; i < index.size()-1; ++i) {
                for(size_t j = i+1; j < index.size(); ++j)
                    idx += (index[j] < index[i] ? 1 : 0);
                idx *= (index.size()-i-1);
            }
            return uniform_strategy[idx];
        }
    } else {
        static std::vector<facul_method_side> const placeholder;
        auto it = direct_access_get(v);
        if (it != nullptr)
            return *it;
        else
            return placeholder;
    }
}

std::vector<facul_method_side> const * & facul_strategies::direct_access_get(std::vector<unsigned int> const & v)
{
    ASSERT_ALWAYS(v.size() == mfb.size());
    unsigned int b = 1, idx = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        idx += v[i]*b;
        b *= mfb[i]+1;
    }
    return direct_access[idx];
}

std::vector<facul_method_side> const * const & facul_strategies::direct_access_get(std::vector<unsigned int> const &v) const
{
    ASSERT_ALWAYS(v.size() == mfb.size());
    unsigned int b = 1, idx = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        idx += v[i]*b;
        b *= mfb[i]+1;
    }
    return direct_access[idx];
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
    // I'm not sure, but it seems that MONTY12 with param=2 (which has
    // nothing special anyway) is not reached by the MONTYTWED12
    // parameterization. Let's pick an arbitrary, other curve from the
    // MONTYTWED12 family, then.
    chain.emplace_back(EC_METHOD, 105, 3255, MONTYTWED12, 2, 1);
#endif

    /* See #30033 ; B1 is not monotonic here, but this is on purpose. The
     * Brent-Suyama curve with sigma=11 (which is isomorphic to the
     * MONTYTWED12 with parameter=1) is exceptional, and yields more
     * primes. Let's put a little bit more effort into it.
     *
     *     sage: load("cofac_utils.sage")
     *     sage: E0,P0=Twed12_parameterization(1, Rationals())
     *     sage: E1,P1=BrentSuyama_parameterization(11, Rationals())
     *     sage: E0.is_isomorphic(E1)
     *     True
     */
    if (n--) {
#ifdef USE_LEGACY_DEFAULT_STRATEGY
        chain.emplace_back(EC_METHOD, 315, 5355, BRENT12, 11, 1);
#else
        chain.emplace_back(EC_METHOD, 315, 5355, MONTYTWED12, 1, 1);
#endif
    }

    /* heuristic strategy where B1 is increased by c*sqrt(B1) at each curve
     *
     * Note that only some addition chains were precomputed in
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
        unsigned long const B2 = (2 * (unsigned int) ((50.0 * B1) / 210.0) + 1) * 105;

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
        std::vector<unsigned long> const & lim,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & mfb,
        std::vector<int> ncurves,
        bool perfectly_sieved,
	const int verbose)
    : facul_strategies_base(lim, lpb, mfb, perfectly_sieved)
{
    int const nsides = lim.size();

    int max_ncurves = -1;
    for(int side = 0 ; side < nsides ; side++) {
        if (ncurves[side] < 0)
            ncurves[side] = nb_curves_with_fbb (B[side], lpb[side], mfb[side]);
        max_ncurves = std::max(max_ncurves, ncurves[side]);
    }

    verbose_fmt_print(0, 2, "# Using default strategy for the cofactorization:");
    for (unsigned int i = 0; i < ncurves.size(); i++) {
        verbose_fmt_print(0, 2, " ncurves{}={}", i, ncurves[i]);
    }
    verbose_fmt_print(0, 2, "\n");

    /* prepare the chain of methods that we want to use in order to
     * factor a number, irrespective of its side.
     */
    auto chain_parameters = facul_strategy_oneside::default_strategy(max_ncurves);

    /* We now need to truncate the list of strategies according to
     * ncurves[i] for 0 <= i < nsides
     *
     * NOTE: This is incompatible with USE_MPQS
     */

    std::vector<std::vector<facul_method::parameters_with_side>> w;

    std::vector<unsigned int> p(nsides);
    std::iota(p.begin(), p.end(), 0);
    do
    {
        /* We iterate over all permutations of [0..nsides-1] in lexicographic
         * order starting from (0,1,...,nsides-1) and using
         * std::next_permutation.
         * Each permutation p corresponds to the case where side p[0] is the
         * largest, side p[1] the second largest, etc..., i.e.,
         *   side p[0] >= side p[1] >= ... >= side p[nsides-1]
         *
         * See the implementation of facul_strategies::operator() to see how the
         * permutation is computed from the sizes of the sides.
         */
        w.emplace_back();
        for(unsigned int side: p) {
            int n = ncurves[side] + (chain_parameters.size() - max_ncurves);
            for(facul_method::parameters const & mp: chain_parameters) {
                if (!n--)
                    break;
                w.back().emplace_back(side, mp);
            }
        }

        /* facul_strategy_oneside::default_strategy is allowed to use
         * default parameters, of course. And it has to abide by the same
         * rules as the strategy files. */
        parameter_sequence_tracker::fill_default_parameters(w.back());
    } while (std::next_permutation(p.begin(), p.end()));

    /* Add all methods to the cache.
     *
     * Note that at least in the present setting, where the method cache
     * is indexed by something that also contains the parameter and so
     * on, we may end up computing an addition chain twice if it so
     * happens that the same B1 is used with two different parameters, or
     * two different EC parameterizations. That is slightly annoying.
     */
    for(auto const & wi: w) {
        for(auto const & mp: wi)
            precompute_method(mp, verbose);
    }

    for(auto const & wi: w) {
        uniform_strategy.emplace_back();
        std::vector<facul_method_side> & cur = uniform_strategy.back();
        for(auto const & mp: wi)
            cur.emplace_back(&precomputed_methods[mp], mp.side);
    }
}
