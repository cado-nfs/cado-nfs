#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <climits>
#include <cerrno>
#include <cstdlib>
#include <cctype>
#include <cstring>

#include <utility>
#include <string>
#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <vector>

#include <hwloc.h>
#include <hwloc/bitmap.h>

#include "cpubinding.hpp"
#include "params.h"
#include "macros.h"
#include "portability.h"

/* This causes all messages to be immediately printed to stderr, instead
 * of being captured to a string */
#define xxxCPUBINDING_DEBUG

/* All output is prefixed by this. */
#define PRE "cpubinding: "

void cpubinding_decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "input-topology-file",
            "simulated topology, only for testing");
    param_list_decl_usage(pl, "input-topology-string",
            "simulated topology, only for testing");
    param_list_decl_usage(pl, "cpubinding", "path to a cpubinding.conf file, or explicit CPU binding string");
}

void cpubinding_lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "cpubinding");
}


/* {{{ sugar */
template<class T>
static std::invalid_argument operator<<(std::invalid_argument const& i, T const& a)
{
    std::ostringstream os;
    os << i.what() << a;
    return std::invalid_argument(os.str());
}
template<typename T>
static std::istream& parse_and_copy_to_list(std::istream& is, std::list<T>& L)
{
    L.clear();
    /* The following is neat, but is does nothing to distinguish eof from
     * parse failure */
#if 0
    std::istream_iterator<T> w0(is);
    std::istream_iterator<T> w1;
    std::ranges::copy(w0, w1, back_inserter(L));
#endif
    /* Here we promise to set failbit on failure. */
    std::string t;
    while (getline(is, t, ' ')) {
        if (t.empty()) continue;
        T v;
        if (std::istringstream(t) >> v) {
            L.push_back(v);
        } else {
            is.clear();
            is.setstate(is.rdstate() | std::ios_base::failbit);
            return is;
        }
    }
    /* EOF, but no fail */
    is.clear();
    is.setstate(is.rdstate() | std::ios_base::eofbit);
    // is.setstate(is.rdstate() | std::ios_base::failbit);
    // is.setstate(is.rdstate() & ~std::ios_base::goodbit);
    return is;
}
/* }}} */

/* {{{ class thread_split understands 2-dimension pairs */
class thread_split {
    std::array<int, 2> t;
    public:
    int const& operator[](int i) const { return t[i]; }
    int& operator[](int i) { return t[i]; }
    explicit operator int() const { return t[0] * t[1]; }
    thread_split() : t { 1, 1 } {}
    thread_split(int t0, int t1) : t { t0, t1 } {}
    explicit thread_split(const int tt[2]) : t { tt[0], tt[1] } {}
    auto operator<=>(const thread_split& o) const = default;
    bool operator==(const thread_split& o) const = default;
#ifdef HAVE_LIBSTDCXX_BUG_114153
    bool operator<(const thread_split& o) const { return operator<=>(o) < 0; }
#endif
    thread_split operator*(thread_split const& o) const {
        return { t[0]*o[0], t[1]*o[1] };
    }
    thread_split& operator*=(thread_split const& o) {
        t[0]*=o[0]; t[1]*=o[1]; return *this;
    }
    thread_split operator/(thread_split const& o) const {
        return { t[0]/o[0], t[1]/o[1] };
    }
    thread_split& operator/=(thread_split const& o) {
        t[0]/=o[0]; t[1]/=o[1]; return *this;
    }
    bool is_divisible_by(thread_split const& o) const {
        return t[0]%o[0]==0 && t[1]%o[1]==0;
    }
};
static std::ostream& operator<<(std::ostream& os, thread_split const& t) {
    os << t[0] << 'x' << t[1];
    return os;
}
static std::istream& operator>>(std::istream& is, thread_split& t) {
    /* a "thr=" may be prepended */
    std::string s;
    if (!(is>>s)) return is;

    const char * digits = "0123456789";
    std::string::size_type x00;
    x00 = s.find_first_of(digits);
    if (x00) {
        if (s[x00-1] == '=') {
            /* we might want to do something with s.substr(0, x00-1) */
            // t.key = s.substr(0, x00-1);
        } else {
            is.setstate(std::ios_base::failbit);
            return is;
        }
    }
    std::string::size_type const x01 = s.find_first_not_of(digits, x00);
    if (x01 == std::string::npos
            || x01 >= s.size()
            || !(std::istringstream(s.substr(x00, x01 - x00)) >> t[0]))
    {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    std::string::size_type const x10 = x01 + 1;
    if (!(std::istringstream(s.substr(x10)) >> t[1])) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    return is;
}
/* }}} */

/* {{{ topology_level is the type describing one level in the processor
 * topology tree */
class topology_level {
public:
    std::string object;
    int n = -1;
    bool has_memory = false;
    topology_level() = default;
    topology_level(std::string const& s, int n, bool has_memory = false)
        : object(s)
        , n(n)
        , has_memory(has_memory)
    {
        if (s == "Socket") object="Package";
    }
    friend std::istream& operator>>(std::istream& is, topology_level& t);
    bool operator<(topology_level const& o) const {
        if (object < o.object) return true;
        return (object == o.object && n < o.n);
    }
    bool operator==(topology_level const& o) const { return object == o.object && n == o.n; }
    bool operator!=(topology_level const& o) const { return !operator==(o); }
};
std::istream& operator>>(std::istream& is, topology_level& t)
{
    std::string s;
    if (!(is>>s)) return is;
    /* This makes sense only with hwloc-2.x, but we might want to avoid
     * having code compiled for hwloc-1.x choke on files that were
     * intended for hwloc-2.x
     */
    if (s == "[NUMANode]") {
        t.has_memory = true;
        return is;
    }
    std::string::size_type const colon = s.find(':');
    if (colon != std::string::npos) {
        t.object = s.substr(0, colon);
        if (t.object == "Socket") t.object="Package";
        if (!(std::istringstream(s.substr(colon+1)) >> t.n)) {
            is.setstate(std::ios_base::failbit);
        }
    }
    return is;
}
static std::ostream& operator<<(std::ostream& os, const topology_level& t)
{
    os << t.object << ":" << t.n;
    if (t.has_memory)
        os << " [NUMANode]";
    return os;
}
/* }}} */

using synthetic_topology = std::list<topology_level>;

/* {{{ dealing with hwloc synthetic topology strings */
static synthetic_topology hwloc_synthetic_topology(hwloc_topology_t topology)
{
    /* too bad hwloc itself has no function for this ! */
    hwloc_obj_t obj = hwloc_get_root_obj(topology);

    if (!obj->symmetric_subtree) {
        throw std::invalid_argument("Cannot output asymetric topology in synthetic format.");
    }

    synthetic_topology result;

#if HWLOC_API_VERSION >= 0x020000
    bool const root_has_memory = obj->memory_arity > 0;
#endif

    for(unsigned int arity = obj->arity ; arity ; arity = obj->arity) {
        obj = obj->first_child;
        char t[64];
        int const d = hwloc_obj_type_snprintf(t, sizeof(t), obj, 1);
        if (d >= (int) sizeof(t))
            throw std::overflow_error("Too long hwloc type name.");
        topology_level T(t, arity);
#if HWLOC_API_VERSION >= 0x020000
        T.has_memory = obj->memory_arity > 0;
#endif
        result.push_back(T);
    }

#if HWLOC_API_VERSION >= 0x020000
    if (root_has_memory)
        result.front().has_memory = true;
#endif

    return result;
}



static std::ostream& operator<<(std::ostream& os, synthetic_topology const& t)
{
    std::ranges::copy(t, std::ostream_iterator<synthetic_topology::value_type>(os, " "));
    return os;
}
static std::istream& operator>>(std::istream& is, synthetic_topology& L)
{
    return parse_and_copy_to_list(is, L);
}

/* }}} */

/* {{{ mapping strings are elementary objects of mapping indications */
struct mapping_string {
    std::string object;
    int group;
    thread_split t;
    mapping_string() : group(0) {}
    mapping_string(std::string const& type, int g, thread_split const& t) : object(type), group(g), t(t) {
        if (type == "Socket") { object = "Package"; }
    }
    mapping_string(std::string const& type, int g, int t[2]) : object(type), group(g), t(t) {
        if (type == "Socket") { object = "Package"; }
    }
    mapping_string(std::string const& type, int g): object(type), group(g) {
        if (type == "Socket") { object = "Package"; }
    }
    bool operator<(const mapping_string& o) const {
        if (object != o.object) return object < o.object;
        if (group != o.group)   return group < o.group;
        return t < o.t;
    }
};
static std::ostream& operator<<(std::ostream& os, const mapping_string& ms)
{
    os << ms.object;
    if (ms.group) os << "*" << ms.group;
    return os << "=>" << ms.t;
}
static std::istream& operator>>(std::istream& is, mapping_string& ms)
{
    ms = mapping_string();
    std::string s;
    if (!(is>>s)) {
        return is;
    }
    std::string::size_type const rel = s.find("=>");
    if (rel == std::string::npos || !(std::istringstream(s.substr(rel+2)) >> ms.t)) {
        is.setstate(std::ios_base::failbit);
    } else {
        std::string::size_type const star = s.find('*');
        if (star == std::string::npos) {
            ms.object = s.substr(0, rel);
            if (ms.object == "Socket") ms.object="Package";
        } else if (star < rel) {
            ms.object = s.substr(0, star);
            if (ms.object == "Socket") ms.object="Package";
            if (!(std::istringstream(s.substr(star+1, rel-(star+1))) >> ms.group)) {
                is.setstate(std::ios_base::failbit);
            }
        } else {
            is.setstate(std::ios_base::failbit);
        }
    }
    return is;
}

static std::ostream& operator<<(std::ostream& os, std::list<mapping_string> const& L)
{
    if (L.empty()) return os;
    std::copy(L.begin(), --L.end(), std::ostream_iterator<mapping_string>(os, " "));
    os << L.back();
    return os;

}

static std::istream& operator>>(std::istream& is, std::list<mapping_string>& L)
{
    return parse_and_copy_to_list(is, L);
}


/* }}} */

/* {{{ matching strings are just the same, with jokers */
struct matching_string : public topology_level {
    private:
        using super = topology_level;
    public:
    std::string joker;
    bool operator<(const matching_string& o) const {
        if (!joker.empty() && !o.joker.empty()) { return joker < o.joker; }
        if (!joker.empty()) { return true; }
        if (!o.joker.empty()) { return true; }
        return (const topology_level&)*this < (const topology_level&)o;
    }
};

static std::istream& operator>>(std::istream& is, matching_string& ms)
{
    char c;
    ms = matching_string();
    if (!(is>>c)) return is;
    if (c == '@') {
        return is >> ms.joker;
    } else {
        is.putback(c);
        return is >> (topology_level&) ms;
    }
}

static std::ostream& operator<<(std::ostream& os, const matching_string& ms)
{
    if (!ms.joker.empty()) return os << '@' << ms.joker;
    return os << (const topology_level&) ms;
}

static std::ostream& operator<<(std::ostream& os, std::list<matching_string> const& L)
{
    if (L.empty()) return os;
    std::copy(L.begin(), --L.end(), std::ostream_iterator<matching_string>(os, " "));
    os << L.back();
    return os;
}

static std::istream& operator>>(std::istream& is, std::list<matching_string>& L)
{
    parse_and_copy_to_list(is, L);
    /* Our parser may either encounter matching strings for hwloc-1.x or hwloc-2.x. In the former case, we have NUMANode objects that we must place differently. In the latter case our parser inserted hwloc-2.x memory objects at temporary places, and we must adjust the std::list accordingly. */
    
#if HWLOC_API_VERSION >= 0x020000
    std::list<matching_string> L2 { };
    /* our section matching must be smart enough to do the right thing
     * with NUMANode matchers that we have here and there. Most often, we
     * want our binding logic to remain sensible to the position of the
     * NUMANode information. So we can't just ditch it.
     *
     * See:
     *
     * https://www.open-mpi.org/projects/hwloc/doc/v2.0.3/a00327.php
     *
     */
    int is_version = 0;
    for(; !L.empty() ;) {
        auto &x = L.front();
        if (x.has_memory) {
            if (is_version == 1) {
                throw std::invalid_argument("")
                    << "Invalid matching std::string:"
                    << "can't have both hwloc-1.x and hwloc-2.x syntax";
            } else if (is_version == 2) {
                throw std::invalid_argument("")
                    << "Invalid matching std::string:"
                    << "found memory objects at two different levels";
            }
            is_version = 2;
            if (L2.empty()) {
                throw std::invalid_argument("")
                    << "Invalid matching std::string:"
                    << "memory objects cannot be on top";
            }
            L2.back().has_memory = true;
            L.pop_front();
            continue;
        }
        if (x.object != "NUMANode") {
            L2.splice(L2.end(), L, L.begin());
            continue;
        }
        if (is_version == 2) {
            throw std::invalid_argument("")
                << "Invalid matching std::string:"
                << "can't have both hwloc-1.x and hwloc-2.x syntax";
        } else if (is_version == 1) {
            throw std::invalid_argument("")
                << "Invalid matching std::string:"
                << "found memory objects at two different levels";
        }
        is_version = 1;
        int const n = x.n;
        
        if (n == 1) {
            /* I think that hwloc-1.x *NEVER* outputs this:
             *
             * - If there's one single level in the machine, we have
             *   no NUMANode in the hierarchy, e.g.:
             *     Socket:2 L2Cache:2 L1Cache:2 Core:1 PU:1
             *
             * - If there are several, they're counted as different
             *   objects, and come earlier in the hierarchy than the item
             *   they contain:
             *     NUMANode:2 Package:1
             *
             * The former goes unchanged, of course. The latter is
             * changed to Package:2 [NUMANode]
             *
             * In a situation where we have two items at the level
             * directly under NUMANode, e.g.:
             *     NUMANode:2 Package:2
             *
             * Then we must insert a Group item as follows:
             *     Group:2 [NUMANode] Package:2
             */

            /* The temptation is to just refuse to parse this. Now in
             * fact, it's also quite easy to form something that is
             * synctactically correct anyway, so we might as well do
             * it.
             */
            if (!L2.empty())
                L2.back().has_memory = true;
            continue;
        }
        L.pop_front();
        if (L.empty()) {
            throw std::invalid_argument("")
                << "Invalid matching std::string:"
                <<" NUMANode must be followed by something";
        }
        auto & next = L.front();
        if (!next.joker.empty()) {
            throw std::invalid_argument("")
                << "Invalid matching std::string:"
                <<" NUMANode followed by a joker would have ambiguous meaning with hwloc-2.x";
        }
        next.n *= n;
        next.has_memory = true;
        L2.splice(L2.end(), L, L.begin());
    }
    if (!is_version) {
        ASSERT_ALWAYS(!L2.empty());
        /* then we have an hwloc-1.x version std::string, let's mark its
         * topmost item as having memory */
        L2.front().has_memory = true;
    }
    L=std::move(L2);
#endif

    return is;
}

/* }}} */

using conf_file = std::map<std::list<matching_string>,
            std::map<thread_split,
                std::list<mapping_string>>>;

/* {{{ Dealing with the configuration file */
static std::istream& operator>>(std::istream& f, conf_file& result)
{
    std::string s;
    result.clear();

    std::pair<conf_file::key_type, conf_file::mapped_type> current;

    for(int lnum = 0; getline(f, s) ; lnum++) {
        std::istringstream is(s);
        thread_split t;
        char c;

        if (!(is>>c)) continue;
        if (c == '#') continue;

        if (c == '[') {
            std::string::size_type i0, i1;
            for(i0 = 0 ; i0 < s.size() && isspace(s[i0]) ; i0++) ;
            for(i1 = s.size() ; --i1 < s.size() && isspace(s[i1]) ; ) ;
            if (i0 == std::string::npos) goto conf_file_parse_error;
            if (i1 == std::string::npos) goto conf_file_parse_error;
            i0++;
            if (i0 >= i1) goto conf_file_parse_error;

            if (!current.first.empty()) {
                if (result.find(current.first) != result.end()) {
                    throw std::invalid_argument("")
                        << "Found two sections with header " << current.first;
                }
                result.insert(std::move(current));
            }

            conf_file::key_type key_in_conf;

            if (!(std::istringstream(s.substr(i0, i1 - i0)) >> key_in_conf)) {
                goto conf_file_parse_error;
            }

            current = { key_in_conf, {} };
            continue;
        }
        is.putback(c);
        if (is >> std::skipws >> t) {
            std::string rhs;
            for( ; (is>>c) && isspace(c) ; ) ;
            is.putback(c);
            getline(is, rhs);
            conf_file::mapped_type::mapped_type v;
            if (rhs == "remove") {
                /* then an empty v is fine */
            } else if (!(std::istringstream(rhs) >> v)) {
                /* try special cases */
                goto conf_file_parse_error;
            }
            if (!current.first.empty()) {
                current.second.insert(std::make_pair(t, v));
            }
            continue;
        }

conf_file_parse_error:
        throw std::invalid_argument("") << "parse error on line " << lnum << ": " << s;
    }
    if (!current.first.empty()) {
        if (result.find(current.first) != result.end()) {
            throw std::invalid_argument("")
                << "Found two sections with header " << current.first;
        }
        result.insert(std::move(current));
    }

    /* reaching EOF is *normal* here ! */
    bool const fail = f.rdstate() & std::ios_base::failbit;
    bool const eof = f.rdstate() & std::ios_base::eofbit;
    if (fail && !eof) {
        return f;
    } else if (eof) {
        f.clear();
    } else {
        f.setstate(f.rdstate() | std::ios_base::failbit);
        f.setstate(f.rdstate() & ~std::ios_base::goodbit);
    }
    return f;
}

#if 0
static std::ostream& operator<<(std::ostream& os, conf_file const& conf)
{
    for(auto it : conf) {
        os << "[" << it.first << "]\n";
        for(auto jt : it.second) {
            os << jt.first << " " << jt.second << "\n";
        }
    }
    return os;
}
#endif
/* }}} */

/* {{{ the matching code within the conf file */
static bool compare_to_section_title(std::ostream& os, synthetic_topology & topology, std::list<matching_string> const& title, int& njokers, std::list<mapping_string>& extra)
{
    njokers = 0;
    extra.clear();

    auto t = topology.begin();
    for(auto const & s : title) {
        if (t == topology.end()) {
            /* no exact match possible. */
            return false;
        }
        if (s.joker == "merge_caches") {
            while (t->object.find("Cache") != std::string::npos) {
                int const g = t->n;
                t = topology.erase(t);
                if (t == topology.end()) {
                    throw std::invalid_argument("@merge_caches failure");
                }
                t->n *= g;
            }
            if (t->object.find("Core") == std::string::npos) {
                os << PRE << "Warning: @merge_caches should encounter one or several Caches above a \"Core\"\n";
                os << PRE << "Warning: got " << *t << " instead, weird (but harmless).\n";
            }
            njokers++;
            continue;
        } else if (s.joker == "group_PU") {
            if (t->object != "PU")
                return false;
            /* This sets the group argument to t->n */
            extra.emplace_back("PU", t->n);
            t++;
            njokers++;
            continue;
        } else if (!s.joker.empty()) {
            throw std::invalid_argument("") << "Bad joker " << s;
        }
        /* now compare the topology level *t with the section title token s */
        if (s != *t) return false;
#if HWLOC_API_VERSION >= 0x020000
        if (s.has_memory != t->has_memory) return false;
#endif
        t++;
    }

    return t == topology.end();
}

/* }}} */

/* {{{ pinning_group */
class pinning_group {
    /* owned pointers. This whole class is only about ownership, really.
     * The only thing we care about is to not leave stray pointers in the
     * wild. */
    hwloc_bitmap_t cpu;
    hwloc_bitmap_t mem;
public:
    pinning_group() {
        cpu = hwloc_bitmap_alloc();
        mem = hwloc_bitmap_alloc();
    }
    ~pinning_group() {
        hwloc_bitmap_free(cpu);
        hwloc_bitmap_free(mem);
    }
    pinning_group(pinning_group const& o) {
        cpu = hwloc_bitmap_dup(o.cpu);
        mem = hwloc_bitmap_dup(o.mem);
    }
    pinning_group(hwloc_obj_t p) {
        cpu = hwloc_bitmap_dup(p->cpuset);
        mem = hwloc_bitmap_dup(p->nodeset);
    }
    pinning_group& operator=(pinning_group const& o) {
        if (this == &o) return *this;
        hwloc_bitmap_copy(cpu, o.cpu);
        hwloc_bitmap_copy(mem, o.mem);
        return *this;
    }
    /* define mergeing operations */
    pinning_group operator+(pinning_group const& o) const {
        pinning_group const r;
        hwloc_bitmap_or(r.cpu, cpu, o.cpu);
        hwloc_bitmap_or(r.mem, mem, o.mem);
        return r;
    }
    pinning_group& operator+=(pinning_group const& o) {
        hwloc_bitmap_or(cpu, cpu, o.cpu);
        hwloc_bitmap_or(mem, mem, o.mem);
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& o, pinning_group const & p);
    int pin(hwloc_topology_t topology, int flags) const {
        return hwloc_set_cpubind(topology, cpu, flags);
    }
};
std::ostream& operator<<(std::ostream& o, pinning_group const & p) {
    char * c, * m;
    hwloc_bitmap_list_asprintf(&c, p.cpu);
    hwloc_bitmap_list_asprintf(&m, p.mem);
    o << "cpu:" << c << " ; mem:" << m;
    free(c);
    free(m);
    return o;
}
/* }}} */


class cpubinder {
    /* {{{ pinning_group_matrices */
    /* This is an internal type for stage_mapping, really. We have a std::list of
     * matrices of pinning groups. All pinning groups are distinct. At the
     * beginning of the stage_mapping processing, they collectively represent
     * the whole system. Later, as chop_off gets called (which might be
     * never), we have a restricted view.
     * The two inner dimensions are fixed by t, the outer one
     * is implied by the length of m.
     */
    class pinning_group_matrices {
        public:
        thread_split t;
        std::vector<pinning_group> m;
        int outer_dimension() const { return m.size() / (int) t; }
    //private:
        pinning_group_matrices(thread_split const& t) : t(t) {}
    public:
        /* flat constructors */
        pinning_group_matrices() = default;
        pinning_group_matrices(pinning_group const& p) : m(1,p) {}
        pinning_group_matrices(std::vector<pinning_group> const& p) : m(p) {}
        /* i <= t[0], j <= t[1] */
        pinning_group& xs(int k, int i, int j) {
            ASSERT_ALWAYS(k >= 0 && k < outer_dimension());
            ASSERT_ALWAYS(i >= 0 && i < t[0]);
            ASSERT_ALWAYS(j >= 0 && j < t[1]);
            return m[(k*t[0]+i)*t[1]+j];
        }
        const pinning_group& xs(int k, int i, int j) const {
            ASSERT_ALWAYS(k >= 0 && k < outer_dimension());
            ASSERT_ALWAYS(i >= 0 && i < t[0]);
            ASSERT_ALWAYS(j >= 0 && j < t[1]);
            return m[(k*t[0]+i)*t[1]+j];
        }
        pinning_group& xs(int i, int j) { return this->xs(0, i, j); }
        const pinning_group& xs(int i, int j) const { return this->xs(0, i, j); }

        pinning_group_matrices& coarsen(thread_split n) {
            if (outer_dimension() % (int) n) {
                throw std::invalid_argument("")
                        << "Cannot coarsen a std::list of " << outer_dimension()
                        << " matrices of pinning groups in blocks of size " << n;
            }
            pinning_group_matrices result(t*n);
            result.m.assign(m.size(), pinning_group());
            for(int k = 0 ; k < outer_dimension() / (int) n ; k++) {
                /* place (int) n blocks */
                for(int n0 = 0 ; n0 < n[0] ; n0++) {
                    for(int n1 = 0 ; n1 < n[1] ; n1++) {
                        int const nn = n0 * n[1] + n1;
                        /* place the elementary blocks at the right places */
                        for(int t0 = 0 ; t0 < t[0] ; t0++) {
                            for(int t1 = 0 ; t1 < t[1] ; t1++) {
                                result.xs(k, n0 * t[0] + t0, n1 * t[1] + t1) =
                                    (pinning_group const&) xs(k*(int)n+nn, t0, t1);
                            }
                        }
                    }
                }
            }
            std::swap(*this, result);
            return *this;
        }
        pinning_group_matrices& chop_off(int n) {
            pinning_group_matrices result(t);
            ASSERT_ALWAYS(outer_dimension() > 0 && outer_dimension() % n == 0);
            for(int k = 0 ; k < outer_dimension() ; k+=n) {
                std::copy(&xs(k,0,0), &xs(k+1,0,0), back_inserter(result.m));
            }
            std::swap(*this, result);
            return *this;
        }
    };
    /* }}} */
    std::ostream& os;
    /* initialized by ctor, freed by dtor */
    hwloc_topology_t topology;

    /* filled by ::read_param_list */
    conf_file cf;
    
    /* This is set by read_param_list */
    bool fake = false;

    /* filled by ::find and ::force */
    synthetic_topology stopo;
    thread_split thr;
    std::list<mapping_string> mapping;

    /* filled by ::stage */
    thread_split coarse;
    pinning_group_matrices coarse_slots;

    public:
    /* we don't want the hwloc private thing be copied around without
     * notice. */
    cpubinder(cpubinder const&) = delete;
    cpubinder&operator=(cpubinder const&) = delete;
    cpubinder(std::ostream& os) : os(os) { hwloc_topology_init(&topology); }
    ~cpubinder() { hwloc_topology_destroy(topology); }
    void read_param_list(cxx_param_list & pl, int want_conf_file);
    bool find(thread_split const& thr);
    void force(thread_split const& ,const char * desc);
    void set_permissive_binding();
    void stage();
    /* this one may be called in MT context */
    void apply(int i, int j) const;
};


void cpubinder::read_param_list(cxx_param_list & pl, int want_conf_file)
{
    /* the first two arguments here are not parsed in the cado-nfs
     * context. It's only used by the helper binary I have for testing
     * the topology matching code */
    const char * topology_file = param_list_lookup_string(pl, "input-topology-file");
    const char * topology_string = param_list_lookup_string(pl, "input-topology-string");
    const char * cpubinding_conf = param_list_lookup_string(pl, "cpubinding");

    /* If we arrive here, then cpubinding_conf is not something which
     * looks like a forced binding, so this should be a config file
     */
    if (cpubinding_conf && want_conf_file) {
        if (std::ifstream(cpubinding_conf) >> cf) {
            os << PRE << "Read configuration from " << cpubinding_conf << "\n";
        } else {
            throw std::invalid_argument("")
                << "Could not read cpubinding conf file "
                << cpubinding_conf;
        }
        // cerr << "Configuration:\n" << cf << "\n";
    }

#if HWLOC_API_VERSION < 0x020000
    unsigned long flags = 0;
#if HWLOC_API_VERSION >= 0x010700
    flags = hwloc_topology_get_flags(topology);
#endif  /* HWLOC_API_VERSION >= 0x010700 */
    /* we must make sure to remove these flags, but it's likely that
     * they're off by default anyway */
    flags &= ~(HWLOC_TOPOLOGY_FLAG_IO_DEVICES | HWLOC_TOPOLOGY_FLAG_IO_BRIDGES);
    hwloc_topology_set_flags(topology, flags);
#endif

    /* {{{ retrieve the topology */
    if (topology_file) {
        int const rc = hwloc_topology_set_xml(topology, topology_file);
        ASSERT_ALWAYS_OR_THROW(rc >= 0, std::invalid_argument);
        fake = true;
    } else if (topology_string) {
        /* hwloc-1.4.1 does not seem to understand "NUMANode" when
         * parsing synthetic strings. With 1.9.1 it works fine. I don't
         * know when exactly that changed.
         * 1.4.1 wants "node" in the synthetic std::string instead. 1.9.1
         * still groks that. Let's do a transformation for pre-1.9, that
         * should keep us safe.
         */
        int rc = -1;
        {
            std::istringstream is(topology_string);
            synthetic_topology stopo;
            if (is >> stopo) {
#if HWLOC_API_VERSION < 0x010b00
                for(auto& x: stopo) {
#if HWLOC_API_VERSION < 0x010900
                    if (x.object == "NUMANode") { x.object="node"; }
#endif
                    if (x.object == "Package") { x.object="Socket"; }
                }
#endif
                std::ostringstream os;
                os << stopo;
                std::string const ss(os.str());
                const char * v = ss.c_str();
                rc = hwloc_topology_set_synthetic(topology, v);
                if (rc < 0)
                    std::cerr << "hwloc_topology_set_synthetic("<< v <<") [mangled for hwloc<1.9] failed\n";
            } else {
                std::cerr << "pre-hwloc-1.9 mangling for hwloc_topology_set_synthetic("<< topology_string <<") failed\n";
                rc = -1;
            }
        }
        ASSERT_ALWAYS(rc >= 0);
        fake = true;
    } else {
        fake = false;
    }
    hwloc_topology_load(topology);
    /* }}} */
}

void cpubinder::set_permissive_binding()
{
    /* stopo and thr must have been set */
    os << PRE << "Selecting permissive mapping:"
        << " " << mapping << "\n";
    mapping.clear();
    mapping_string const top(stopo.front().object, stopo.front().n, thr);
    mapping.push_back(top);
}

/* {{{ finds a mapping for thr. Return true if successful. */
/* This sets the fields [stopo] [thr] [mapping].
 */
bool cpubinder::find(thread_split const& thr)
{
    this->thr = thr;
    stopo = hwloc_synthetic_topology(topology);
    os << PRE << "Hardware: " << stopo << "\n";
    os << PRE << "Target split: " << thr << "\n";

    struct {
        int n;
        std::pair<conf_file::key_type, conf_file::mapped_type> it;
        std::list<mapping_string> e;
        synthetic_topology s;
    } best;
    int nm=0;
    best.n = INT_MAX;
    for(auto const & it : cf) {
        int n;
        std::list<mapping_string> e;
        synthetic_topology s = stopo;
        if (compare_to_section_title(os, s, it.first, n, e)) {
            nm++;
            if (n < best.n) {
                best.s = s;
                best.n = n;
                best.it = it;
                best.e = e;
            } else if (n == best.n) {
                os << PRE << "Found two matches with same accuracy level\n";
                os << PRE << "First match:\n" << best.it.first << "\n";
                os << PRE << "Second match:\n" << it.first << "\n";
                os << PRE << "First match wins.\n";
            }

        }
    }
    if (best.n == INT_MAX)
        return false;
    os << PRE << "config match: [" << best.it.first << "]\n";
    if (nm > 1) {
        os << PRE << "(note: " << (nm-1) << "other (possibly looser) matches)\n";
    }
    auto jt = best.it.second.find(thr);
    if (jt == best.it.second.end())
        return false;
    mapping = jt->second;
    if (!mapping.empty()) {
        for( ; !best.e.empty() ; best.e.pop_front()) {
            if (mapping.back().object != best.e.front().object)
                break;
            os << PRE << "not appending " << best.e.front() << " since " << mapping.back().object << " is present\n";
        }
        mapping.splice(mapping.end(), best.e);
    } else {
        os << PRE << "Found \"remove\" mapping in config file.\n";
        set_permissive_binding();
    }
    stopo = best.s;
    return true;
}
/* }}} */

void cpubinder::force(thread_split const& t, const char * desc)/*{{{*/
{
    thr = t;
    stopo = hwloc_synthetic_topology(topology);
    if (strcmp(desc, "remove") == 0 || strlen(desc) == 0) {
        os << PRE << "As per the provided std::string \""<<desc<<"\","
            << " applying permissive mapping:"
            << " " << mapping << "\n";
        set_permissive_binding();
    } else {
        std::istringstream(desc) >> mapping;
    }
}
/*}}}*/

/* {{{ stage_mapping */
/* This fills the coarse and coarse_slots fields */
void cpubinder::stage()
{
    /* First gather all PUs in an ordered sequence which matches the
     * topology tree. Note that the next_cousin member function is really
     * perfect for that */
    int const depth = hwloc_topology_get_depth(topology);
    /*
    int npu = hwloc_get_nbobjs_by_depth(topology, depth-1);
    int g=1;
    for(auto it = stopo.begin() ; it != stopo.end() ; g*=it++->n) ;
    ASSERT_ALWAYS(g == npu);
    */
    std::vector<pinning_group> slots;
    for(hwloc_obj_t pu = hwloc_get_obj_by_depth(topology, depth-1, 0) ;
            pu != nullptr;
            pu = pu->next_cousin) slots.emplace_back(pu);


    auto rt = stopo.rbegin()  ;
    auto jt = mapping.rbegin();

    /* rtn is the number of not-yet-assigned (groups of) objects at the
     * current level */
    int rtn = rt->n;

    bool stars = true;

    os << PRE << "Reduced topology: " << stopo << "\n";
    os << PRE << "Applying mapping: " << mapping << "\n";

#if HWLOC_API_VERSION >= 0x020000
    int const numa_depth = hwloc_get_memory_parents_depth(topology);
    std::string numa_replace;
    {
        hwloc_obj_t pu = hwloc_get_obj_by_depth(topology, numa_depth, 0);
        char t[64];
        int const d = hwloc_obj_type_snprintf(t, sizeof(t), pu, 1);
        if (d >= (int) sizeof(t))
            throw std::overflow_error("Too long hwloc type name.");
        numa_replace = t;
    }
    for(auto & x : mapping) {
        if (x.object == "NUMANode") {
            x.object = numa_replace;
            os << PRE << "Found legacy hwloc-1.x mapping std::string,"
                << "replacing NUMANode by " << numa_replace << "\n";
        }
    }
#endif

    for( ;jt != mapping.rend() ; jt++) {
        if (!jt->group) {
            if (stars)
                coarse_slots = pinning_group_matrices(slots);
            stars = false;
        }
        if (!stars) {
            if (jt->group) {
                throw std::invalid_argument("") << "Wrongly placed * in " << mapping;
            }
            if ((int) jt->t == 1) continue;
        }
        /* If we have a match, that's good, let's keep it, even if
         * the remaining item of this name in the topology has count
         * one. If we don't have a match, then we're allowed to move
         * up until we find one. If we jump over levels of non-trivial
         * arity in that process, we'll act accordingly.
         */
        int hidden_grouping = 1;
        for( ; rt != stopo.rend() && rt->object != jt->object ; ) {
            if (rtn > 1) {
                /* we have no way to silence this warning at the moment.
                 * Maybe add an explicit syntax like Core/2, or something ?
                 */
                if (stars) {
                    os << PRE << "Warning: while applying " << *jt
                        << ": this implicitly merges " << rtn
                        << " " << rt->object << "-level (groups of) objects\n";
                    hidden_grouping *= rtn;
                } else {
                    os << PRE << "Warning: while applying " << *jt
                        << ": we have only " << rtn
                        << " " << rt->object << " scheduled"
                        << " out of " << rt->n << "\n";
                    /* well, do it, then ! */
                    coarse_slots.chop_off(rtn);
                }
            }
            if (++rt == stopo.rend())
                break;

            rtn = rt->n;
        }

        if (rt == stopo.rend())
            throw std::invalid_argument("")
                << "Hit end of hardware description"
                << " while applying " << *jt;
        bool check;

        if (stars) {
            /* we coarsen the thread group */
            coarse *= jt->t;
            check = rtn % jt->group == 0 && thr.is_divisible_by(coarse);
        } else {
            check = rtn % (int) jt->t == 0;
        }

        if (!check) {
            throw std::invalid_argument("")
                    << *rt << " ("<<rtn<<" left)"
                    << " is not correct while applying " << *jt;
        }

        if (stars) {
            /* we need to coarsen the PU groups by a factor of [group] */
            for(unsigned int i = 0, j = 0 ; j != slots.size() ; i++) {
                slots[i] = slots[j++];
                for(int k = 1 ; k < jt->group * hidden_grouping ; k++) {
                    slots[i] += slots[j++];
                }
            }
            slots.erase(slots.begin() + slots.size() / (jt->group * hidden_grouping), slots.end());
            /* and we reduce the number of items in the topology */
            rtn /= jt->group;
        } else {
            coarse_slots.coarsen(jt->t);
            rtn /= (int) jt->t;
        }
    }
    if (stars)
        coarse_slots = pinning_group_matrices(slots);
    stars = false;
    for( ; ; ) {
        if (rtn > 1) {
            /* we have no way to silence this warning at the moment.
             * Maybe add an explicit syntax like Core/2, or something ?
             */
            os << PRE << "Warning: completed mapping uses only " << (rt->n/rtn) << " " << rt->object << " out of " << rt->n << "\n";
            coarse_slots.chop_off(rtn);
        }
        if (++rt == stopo.rend())
            break;
        rtn = rt->n;
    }

    ASSERT_ALWAYS(coarse_slots.outer_dimension() == 1);

    // os << PRE << "Coarse slots now organized in " << coarse_slots.outer_dimension() << " matrices of size " << coarse_slots.t << "\n";

    /* We should have a matrix of thread groups whose dimension is
     * thr/coarse, each grouping coarse threads */

    os << PRE << "Threads organized as " << thr/coarse << " blocks of dimensions " << coarse << "\n";

    if (coarse_slots.t != thr/coarse) {
        throw std::invalid_argument("") << "mapping does achieve desired split"
            << " (" << coarse_slots.t*coarse << ", wanted " << thr << ")\n";
    }
    for(int i = 0 ; i < thr[0] ; i+=coarse[0]) {
        for(int j = 0 ; j < thr[1] ; j+=coarse[1]) {
            std::ostringstream oos;
            if ((int) coarse > 1) {
                oos << "threads ("<<i<<","<<j<<")"
                    << " to ("<<i+coarse[0]-1<<","<<j+coarse[1]-1<<")";
            } else {
                oos << "thread ("<<i<<","<<j<<")";
            }
            pinning_group const p = coarse_slots.xs(i/coarse[0], j/coarse[1]);
            os << PRE << "" << oos.str() << " -> " << p << "\n";
        }
    }
    if (fake)
        os << PRE << "NOTE: since this is a fictitious hardware description, pinning will not be done for real\n";
    os << PRE << "Done.\n";
}
/* }}} */



/* This returns an opaque pointer to data which will be used to perform
 * the actual cpu binding. This function must be called in
 * single-threaded context.
 * 
 * This returns nullptr if cpubinding failed.
 */

void * cpubinding_get_info(char ** messages, cxx_param_list & pl, unsigned int tt0, unsigned int tt1)
{
    const char * conf = param_list_lookup_string(pl, "cpubinding");
    thread_split const thr(tt0, tt1);

    if (conf && (strcmp(conf, "no") == 0 || strcmp(conf, "none")==0)) {
        if (messages) {
            *messages = strdup("cpubinding disabled by cmdline switch\n");
        }
        return nullptr;
    }

#ifdef  CPUBINDING_DEBUG
    std::ostream& os(std::cerr);
#else   /* CPUBINDING_DEBUG */
    std::ostringstream os;
#endif  /* CPUBINDING_DEBUG */

    auto * cb = new cpubinder(os);

    int const force = conf && (strstr(conf, "=>") || strcmp(conf, "remove") == 0 || strlen(conf) == 0);

    try {
        cb->read_param_list(pl, !force);
        if (force) {
            cb->force(thr, conf);
            cb->stage();
        } else if (cb->find(thr)) {
            cb->stage();
        } else {
            os << PRE << "no mapping found\n";
            delete cb;
            cb = nullptr;
        }
    } catch (std::invalid_argument const& e) {
        os << PRE << "Failed on error:\n"
            << PRE << "  " << e.what() << "\n";
        delete cb;
        cb = nullptr;
    }

    if (messages) {
#ifndef CPUBINDING_DEBUG
        *messages = strdup(os.str().c_str());
#else
        *messages = nullptr;
#endif  /* CPUBINDING_DEBUG */
    }
    return static_cast<void*>(cb);
}

void cpubinder::apply(int i, int j) const
{
    ASSERT_ALWAYS(i < thr[0]);
    ASSERT_ALWAYS(j < thr[1]);
    int const ti = i / coarse[0];
    int const tj = j / coarse[1];
    pinning_group const p = coarse_slots.xs(ti, tj);
    // cout << "Pinning thread ("<<i<<","<<j<<") to " << p << "\n";
    int const rc = p.pin(topology, HWLOC_CPUBIND_THREAD);
    if (rc < 0) {
        std::cerr << "Pinning thread ("<<i<<","<<j<<") to " << p << ": " << strerror(errno) << "\n";
    }
}

/* perform the actual pinning. This must be called for each thread */
void cpubinding_do_pinning(void * pinning_info_pre, int i, int j)
{
    if (pinning_info_pre == nullptr) return;
    auto* cb = static_cast<cpubinder*>(pinning_info_pre);

    cb->apply(i, j);
}

/* free the opaque pointer */
void cpubinding_free_info(void * pinning_info_pre, unsigned int, unsigned int)
{
    if (pinning_info_pre == nullptr) return;
    auto* cb = static_cast<cpubinder*>(pinning_info_pre);
    delete cb;
}
