#ifndef CADO_LINALG_BWC_BWC_FILENAMES_HPP
#define CADO_LINALG_BWC_BWC_FILENAMES_HPP

#include <cstdio>

#include <string>
#include <array>
#include <stdexcept>
#include <ostream>
#include <type_traits>
#include <vector>

#include "fmt/base.h"
#include "fmt/format.h"
#include "fmt/ostream.h"

#include "utils_cxx.hpp"
#include "macros.h"

/* still missing files:
 *      X
 *      K
 * missing code:
 *      pattern_*
 */
struct bwc_file_base : public decomposed_path {// {{{
    template<typename T>
        struct helper {
            static bool match(T &, std::string const & filename);
            static std::vector<T> ls(std::string const & dirname);
        };

    template<typename T>
        static std::vector<T> ls(std::string const & dirname)
        {
            return helper<T>::ls(dirname);
        }

    template<typename T>
        static bool match(T & v, std::string const & filename)
        requires std::is_base_of_v<bwc_file_base, T>
        {
            return helper<T>::match(v, filename);
        }

    /* This makes it possible to add the parsed filename to a container
     * with the good value type, and then return a boolean if somethin
     * has happened
     */
    template<typename Ctr>
        static bool match(Ctr & c, std::string const & filename)
        requires std::is_base_of_v<bwc_file_base, typename Ctr::value_type>
        {
            typename Ctr::value_type v;
            if (!match<typename Ctr::value_type>(v, filename))
                return false;
            c.emplace_back(std::move(v));
            return true;
        }

    template<typename T>
    static T parse(std::string filename) {
        T res;
        if (!match(res, filename))
            throw std::invalid_argument(
                    fmt::format("filename cannot be parsed: {}", filename));
        return res;
    }
};
// }}}

struct bwc_column_range {// {{{
    unsigned int j0 = 0;
    unsigned int j1 = 0;
    bwc_column_range() = default;
    bwc_column_range(unsigned int j0, unsigned int j1) : j0(j0), j1(j1) {}
    unsigned int operator[](unsigned int i) const {
        ASSERT_ALWAYS(i <= 1);
        return i ? j1 : j0;
    }
    bool operator<(bwc_column_range const & b) const /*{{{*/
    {
        int r;
        r = (j0 > b.j0) - (b.j0 > j0); if (r) return r < 0;
        r = (j1 > b.j1) - (b.j1 > j1); if (r) return r < 0;
        return false;
    }/*}}}*/
};
// }}}
struct bwc_solution_range {// {{{
    unsigned int s0 = 0;
    unsigned int s1 = 0;
    bwc_solution_range() = default;
    bwc_solution_range(unsigned int s0, unsigned int s1) : s0(s0), s1(s1) {}
    unsigned int operator[](unsigned int i) const {
        ASSERT_ALWAYS(i <= 1);
        return i ? s1 : s0;
    }
    bool operator<(bwc_solution_range const & b) const /*{{{*/
    {
        int r;
        r = (s0 > b.s0) - (b.s0 > s0); if (r) return r < 0;
        r = (s1 > b.s1) - (b.s1 > s1); if (r) return r < 0;
        return false;
    }/*}}}*/
};
// }}}
struct bwc_iteration_range {// {{{
    unsigned int n0 = 0;
    unsigned int n1 = 0;
    bwc_iteration_range() = default;
    bwc_iteration_range(unsigned int n0, unsigned int n1) : n0(n0), n1(n1) {}
    unsigned int operator[](unsigned int i) const {
        ASSERT_ALWAYS(i <= 1);
        return i ? n1 : n0;
    }
    bool operator<(bwc_iteration_range const & b) const /*{{{*/
    {
        int r;
        r = (n0 > b.n0) - (b.n0 > n0); if (r) return r < 0;
        r = (n1 > b.n1) - (b.n1 > n1); if (r) return r < 0;
        return false;
    }/*}}}*/
};
// }}}

struct bwc_V_file // {{{
    : public bwc_file_base
    , public bwc_column_range
{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_V_file() = default;
    explicit bwc_V_file(std::string const & filename)
        : bwc_V_file(bwc_file_base::parse<bwc_V_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int n = 0;

    static std::string pattern(unsigned int n);

    bool operator<(bwc_V_file const & b) const /*{{{*/
    {
        bwc_column_range const & aj = *this;
        bwc_column_range const & bj = b;
        int r;
        r = (bj < aj) - (aj < bj); if (r) return r < 0;
        r = (b.n < n) - (n < b.n); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_K_file // {{{
    : public bwc_file_base
    , public bwc_solution_range
{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_K_file() = default;
    explicit bwc_K_file(std::string const & filename)
        : bwc_K_file(bwc_file_base::parse<bwc_K_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int n = 0;

    static std::string pattern(unsigned int n);

    bwc_K_file(bwc_solution_range const & s, unsigned int n)
        : bwc_solution_range(s)
        , n(n)
    {
        /* Do this just for consistency */
        decomposed_path::operator=(decomposed_path(std::string(*this)));
    }

    bool operator<(bwc_K_file const & b) const /*{{{*/
    {
        bwc_solution_range const & as = *this;
        bwc_solution_range const & bs = b;
        int r;
        r = (bs < as) - (as < bs); if (r) return r < 0;
        r = (b.n < n) - (n < b.n); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_S_file // {{{
    : public bwc_file_base
    , public bwc_solution_range
    , public bwc_iteration_range
{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_S_file() = default;
    explicit bwc_S_file(std::string const & filename)
        : bwc_S_file(bwc_file_base::parse<bwc_S_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}

    // static std::string pattern(bwc_solution_range const &);
    static std::string pattern(bwc_iteration_range const &);

    bool operator<(bwc_S_file const & b) const /*{{{*/
    {
        bwc_iteration_range const & bn = b;
        bwc_iteration_range const & an = *this;
        bwc_solution_range const & bs = b;
        bwc_solution_range const & as = *this;
        int r;
        r = (bn < an) - (an < bn); if (r) return r < 0;
        r = (bs < as) - (as < bs); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_Cv_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_Cv_file() = default;
    explicit bwc_Cv_file(std::string const & filename)
        : bwc_Cv_file(bwc_file_base::parse<bwc_Cv_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int j0 = 0;
    unsigned int j1 = 0;
    unsigned int stretch = 0;

    static std::string pattern(unsigned int n);

    std::array<unsigned int, 2> column_range() const {
        return { j0, j1 };
    }

    bool operator<(bwc_Cv_file const & b) const /*{{{*/
    {
        int r;
        r = (j1 > b.j1) - (b.j1 > j1); if (r) return r < 0;
        r = (j0 > b.j0) - (b.j0 > j0); if (r) return r < 0;
        r = (stretch > b.stretch) - (b.stretch > stretch); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_Cd_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_Cd_file() = default;
    explicit bwc_Cd_file(std::string const & filename)
        : bwc_Cd_file(bwc_file_base::parse<bwc_Cd_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int j0 = 0;
    unsigned int j1 = 0;
    unsigned int stretch = 0;

    static std::string pattern(unsigned int n);

    std::array<unsigned int, 2> column_range() const {
        return { j0, j1 };
    }

    bool operator<(bwc_Cd_file const & b) const /*{{{*/
    {
        int r;
        r = (j1 > b.j1) - (b.j1 > j1); if (r) return r < 0;
        r = (j0 > b.j0) - (b.j0 > j0); if (r) return r < 0;
        r = (stretch > b.stretch) - (b.stretch > stretch); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_Cr_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_Cr_file() = default;
    explicit bwc_Cr_file(std::string const & filename)
        : bwc_Cr_file(bwc_file_base::parse<bwc_Cr_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int nchecks = 0;

    bool operator<(bwc_Cr_file const & b) const /*{{{*/
    {
        int r;
        r = (nchecks > b.nchecks) - (b.nchecks > nchecks); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_Ct_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_Ct_file() = default;
    explicit bwc_Ct_file(std::string const & filename)
        : bwc_Ct_file(bwc_file_base::parse<bwc_Ct_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int nchecks = 0;
    unsigned int m = 0;

    bool operator<(bwc_Ct_file const & b) const /*{{{*/
    {
        int r;
        r = (nchecks > b.nchecks) - (b.nchecks > nchecks); if (r) return r < 0;
        r = (m > b.m) - (b.m > m); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_A_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_A_file() = default;
    explicit bwc_A_file(std::string const & filename)
        : bwc_A_file(bwc_file_base::parse<bwc_A_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int j0 = 0;
    unsigned int j1 = 0;
    unsigned int n0 = 0;
    unsigned int n1 = 0;

    std::string pattern_column_range();
    std::string pattern_iteration_range();

    std::array<unsigned int, 2> column_range() const {
        return { j0, j1 };
    }
    std::array<unsigned int, 2> iteration_range() const {
        return { n0, n1 };
    }

    bool operator<(bwc_A_file const & b) const /*{{{*/
    {
        int r;
        r = (j1 > b.j1) - (b.j1 > j1); if (r) return r < 0;
        r = (j0 > b.j0) - (b.j0 > j0); if (r) return r < 0;
        r = (n1 > b.n1) - (b.n1 > n1); if (r) return r < 0;
        r = (n0 > b.n0) - (b.n0 > n0); if (r) return r < 0;
        return false;
    }/*}}}*/
};

// }}}
struct bwc_F_file : public bwc_file_base {// {{{
    using decomposed_path::basename;// {{{
    using decomposed_path::dirname;
    using bwc_file_base::match;
    explicit operator std::string() const;
    bwc_F_file() = default;
    explicit bwc_F_file(std::string const & filename)
        : bwc_F_file(bwc_file_base::parse<bwc_F_file>(filename))
    {}
    private:
    friend struct bwc_file_base;
    std::istream & match_fields(std::istream &);

    public:// }}}
    unsigned int j0 = 0;
    unsigned int j1 = 0;
    unsigned int s0 = 0;
    unsigned int s1 = 0;

    std::string pattern_column_range();
    std::string pattern_iteration_range();

    std::array<unsigned int, 2> column_range() const {
        return { j0, j1 };
    }
    std::array<unsigned int, 2> solution_range() const {
        return { s0, s1 };
    }

    bool operator<(bwc_F_file const & b) const /*{{{*/
    {
        int r;
        r = (s1 > b.s1) - (b.s1 > s1); if (r) return r < 0;
        r = (s0 > b.s0) - (b.s0 > s0); if (r) return r < 0;
        r = (j1 > b.j1) - (b.j1 > j1); if (r) return r < 0;
        r = (j0 > b.j0) - (b.j0 > j0); if (r) return r < 0;
        return false;
    }/*}}}*/
};

/* this proxy is handy */
static inline std::unique_ptr<FILE, delete_FILE> fopen_helper(bwc_file_base const & filename, const char * mode, bool accept_failure = false)
{
    return fopen_helper(std::string(filename), mode, accept_failure);
}

// }}}

static inline std::ostream& operator<<(std::ostream& os, bwc_file_base const & x) {
    return os << std::string(x);
}
namespace fmt {
    template <> struct formatter<bwc_file_base>: ostream_formatter {};
    template <> struct formatter<bwc_S_file>: ostream_formatter {};
    template <> struct formatter<bwc_V_file>: ostream_formatter {};
    template <> struct formatter<bwc_K_file>: ostream_formatter {};
    template <> struct formatter<bwc_Cv_file>: ostream_formatter {};
    template <> struct formatter<bwc_Cd_file>: ostream_formatter {};
    template <> struct formatter<bwc_Cr_file>: ostream_formatter {};
    template <> struct formatter<bwc_Ct_file>: ostream_formatter {};
    template <> struct formatter<bwc_A_file>: ostream_formatter {};
    template <> struct formatter<bwc_F_file>: ostream_formatter {};
}

extern template struct bwc_file_base::helper<bwc_S_file>;
extern template struct bwc_file_base::helper<bwc_V_file>;
extern template struct bwc_file_base::helper<bwc_K_file>;
extern template struct bwc_file_base::helper<bwc_Cv_file>;
extern template struct bwc_file_base::helper<bwc_Cd_file>;
extern template struct bwc_file_base::helper<bwc_Cr_file>;
extern template struct bwc_file_base::helper<bwc_Ct_file>;
extern template struct bwc_file_base::helper<bwc_A_file>;
extern template struct bwc_file_base::helper<bwc_F_file>;
#endif	/* LINALG_BWC_BWC_FILENAMES_HPP_ */
