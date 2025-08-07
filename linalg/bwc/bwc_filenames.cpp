#include "cado.h" // IWYU pragma: keep

#include <string>
#include <sstream>
#include <istream>
#include <locale>
#include <ios>
#include <vector>
#include <algorithm>

#include <dirent.h>

#include "fmt/format.h"

#include "bwc_filenames.hpp"
#include "utils_cxx.hpp"
#include "macros.h"

/* TODO: all pattern functions are missing!
    V
    std::string pattern_column_range();
    std::string pattern_iteration();
    S
    std::string pattern_solution_range();
    std::string pattern_iteration_range();
    Cv
    std::string pattern_column_range();
    std::string pattern_stretch();
    Cd
    std::string pattern_column_range();
    std::string pattern_stretch();
    A
    std::string pattern_column_range();
    std::string pattern_iteration_range();
    F
    std::string pattern_column_range();
    std::string pattern_iteration_range();
    */

// {{{ This is the generic part
/* This is good to recognize digits in file names */
struct digits_locale: std::ctype<char>
{
    digits_locale(): std::ctype<char>(get_table()) {}

    static std::ctype_base::mask const* get_table()
    {
        static std::vector<std::ctype_base::mask>
            rc(std::ctype<char>::table_size,std::ctype_base::space);
        std::fill(&rc['0'], &rc['9'+1], std::ctype_base::digit);
        return rc.data();
    }
};

template<typename T>
bool
bwc_file_base::helper<T>::match(T & v, std::string const & name)
{
    (decomposed_path&) v = decomposed_path(name);
    std::istringstream is(v.basename());
    is.imbue(std::locale(std::locale(), new digits_locale()));
    is >> std::noskipws;
    v.match_fields(is);
    return is.eof() && !is.fail();
}

template<typename T>
std::vector<T>
bwc_file_base::helper<T>::ls(std::string const & dirname)
{
    std::vector<T> res;
    DIR * dir = opendir(dirname.c_str());
    DIE_ERRNO_DIAG(!dir, "opendir(%s)", dirname.c_str());
    struct dirent * de;
    for( ; (de = readdir(dir)) != nullptr ; )
        bwc_file_base::match(res, de->d_name);
    closedir(dir);
    sort(res.begin(), res.end());
    return res;
}


// }}}

bwc_V_file::operator std::string() const// {{{
{
    return fmt::format("V{}-{}.{}", j0, j1, n);
}

std::istream & bwc_V_file::match_fields(std::istream & is)
{
    return is >> expect("V") >> j0 >> expect("-") >> j1 >> expect(".") >> n;
}

std::string bwc_V_file::pattern(unsigned int n)
{
    return fmt::format("V{}-{}.{}", "{}", "{}", n);
}

template struct bwc_file_base::helper<bwc_V_file>;
// }}}
bwc_K_file::operator std::string() const// {{{
{
    return fmt::format("K.sols{}-{}.{}", s0, s1, n);
}

std::istream & bwc_K_file::match_fields(std::istream & is)
{
    return is
        >> expect("K.sols") >> s0 >> expect("-") >> s1
        >> expect(".") >> n;
}

std::string bwc_K_file::pattern(unsigned int n)
{
    return fmt::format("K.sols{}-{}.{}", "{}", "{}", n);
}

template struct bwc_file_base::helper<bwc_K_file>;
// }}}
bwc_S_file::operator std::string() const// {{{
{
    return fmt::format("S.sols{}-{}.{}-{}", s0, s1, n0, n1);
}

std::istream & bwc_S_file::match_fields(std::istream & is)
{
    return is
        >> expect("S.sols") >> s0 >> expect("-") >> s1
        >> expect(".") >> n0 >> expect("-") >> n1;
}

/*
std::string bwc_S_file::pattern(bwc_solution_range const & S)
{
    return fmt::format("S.sols{}-{}.{}-{}", S[0], S[1], "{}", "{}");
}
*/
std::string bwc_S_file::pattern(bwc_iteration_range const & N)
{
    return fmt::format("S.sols{}-{}.{}-{}", "{}", "{}", N[0], N[1]);
}

template struct bwc_file_base::helper<bwc_S_file>;
// }}}
bwc_Cv_file::operator std::string() const// {{{
{
    return fmt::format("Cv{}-{}.{}", j0, j1, stretch);
}

std::istream & bwc_Cv_file::match_fields(std::istream & is)
{
    return is
        >> expect("Cv") >> j0 >> expect("-") >> j1
        >> expect(".") >> stretch;
}

std::string bwc_Cv_file::pattern(unsigned int stretch)
{
    return fmt::format("Cv{}-{}.{}", "{}", "{}", stretch);
}
template struct bwc_file_base::helper<bwc_Cv_file>;
// }}}
bwc_Cd_file::operator std::string() const// {{{
{
    return fmt::format("Cd{}-{}.{}", j0, j1, stretch);
}

std::istream & bwc_Cd_file::match_fields(std::istream & is)
{
    return is
        >> expect("Cd") >> j0 >> expect("-") >> j1
        >> expect(".") >> stretch;
}

std::string bwc_Cd_file::pattern(unsigned int stretch)
{
    return fmt::format("Cd{}-{}.{}", "{}", "{}", stretch);
}
template struct bwc_file_base::helper<bwc_Cd_file>;
// }}}
bwc_Cr_file::operator std::string() const// {{{
{
    return fmt::format("Cr0-{0}.0-{0}",  nchecks);
}

std::istream & bwc_Cr_file::match_fields(std::istream & is)
{
    unsigned int nc2;
    is >> expect("Cr0-") >> nchecks >> expect(".0-") >> nc2;
    if (!is) return is;
    if (nc2 != nchecks)
        is.setstate(std::ios_base::failbit);
    return is;
}

template struct bwc_file_base::helper<bwc_Cr_file>;
// }}}
bwc_Ct_file::operator std::string() const// {{{
{
    return fmt::format("Ct0-{}.0-{}",  nchecks, m);
}

std::istream & bwc_Ct_file::match_fields(std::istream & is)
{
    return is >> expect("Ct0-") >> nchecks >> expect(".0-") >> m;
}

template struct bwc_file_base::helper<bwc_Ct_file>;
// }}}
bwc_A_file::operator std::string() const// {{{
{
    return fmt::format("A{}-{}.{}-{}", j0, j1, n0, n1);
}

std::istream & bwc_A_file::match_fields(std::istream & is)
{
    return is
        >> expect("A") >> j0 >> expect("-") >> j1
        >> expect(".") >> n0 >> expect("-") >> n1;
}

template struct bwc_file_base::helper<bwc_A_file>;
// }}}
bwc_F_file::operator std::string() const// {{{
{
    return fmt::format("F.sols{}-{}.{}-{}", s0, s1, j0, j1);
}

std::istream & bwc_F_file::match_fields(std::istream & is)
{
    return is
        >> expect("F.sols") >> s0 >> expect("-") >> s1
        >> expect(".") >> j0 >> expect("-") >> j1;
}

template struct bwc_file_base::helper<bwc_F_file>;
// }}}
