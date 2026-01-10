#ifndef CADO_HWLOC_AUX_HPP
#define CADO_HWLOC_AUX_HPP

#include <cstdlib>
#include <cstdio>

#include <type_traits>
#include <utility>

#include <hwloc.h>

#include "fmt/base.h"
#include "macros.h"

extern "C" {
int hwloc_aux_get_depth_from_string(hwloc_topology_t topology, const char * desc);
}

/* Note that most hwloc structures are pimpls (hwloc_topology_t is a
 * pointer to the incomplete type struct hwloc_topology -- same goes
 * for hwloc_obj_t for instance). It has various annoying
 * consequences, not the least of which is the fact that we can't
 * type-pun. (e.g. make cxx_hwloc_blah contain a struct hwloc_blah[1],
 * and given an interface that returns an hwloc_blah_t ==
 * struct hwloc_blah * with name foo, return the reference formed by
 * *(cxx_hwlo_blah *)foo. We need temporary objects to do that, which is
 * a pain
 */

struct cxx_hwloc_topology {
    hwloc_topology_t x = nullptr;
    static_assert(std::is_same_v<hwloc_topology *, hwloc_topology_t>);
    /* it's ugly. Unfortunately the hwloc interface is C only and has
     * undesirable non-const prototypes. We want to work around this and
     * propose const alternatives */
    private:
    hwloc_topology_t x_nonconst() const {
        return const_cast<hwloc_topology_t>(x);
    }
    public:
    cxx_hwloc_topology() { hwloc_topology_init(&x); }
    ~cxx_hwloc_topology() { hwloc_topology_destroy(x); }
    explicit cxx_hwloc_topology(hwloc_topology const * o) {
        hwloc_topology_dup(&x, const_cast<hwloc_topology *>(o));
    }
    cxx_hwloc_topology(cxx_hwloc_topology && o) noexcept
        : x(o.x)
    {
        o.x = nullptr;
    }
    cxx_hwloc_topology& operator=(cxx_hwloc_topology && o) noexcept {
        std::swap(x, o.x);
        return *this;
    }
    cxx_hwloc_topology(cxx_hwloc_topology const & o) {
        hwloc_topology_dup(&x, o.x);
    }
    cxx_hwloc_topology & operator=(cxx_hwloc_topology const & o) {
        if (&o != this) {
            hwloc_topology_destroy(x);
            hwloc_topology_dup(&x, o.x);
        }
        return *this;
    }
    cxx_hwloc_topology & operator=(hwloc_topology_t ox) {
        if (ox != x) {
            hwloc_topology_destroy(x);
            hwloc_topology_dup(&x, ox);
        }
        return *this;
    }
    /* alright. we "should" try to propagate constness as we do below.
     * But the sad thing is that nothing in the hwloc interface is really
     * const-aware as of now. So it's really a major pain to do so, and
     * we're making our life more difficult than it is.
     */
    /*
    operator hwloc_topology_t() { return x; }
    operator hwloc_topology const *() const { return x; }
    hwloc_topology_t operator->() { return x; }
    hwloc_topology const * operator->() const { return x; }
    */
    /* the following are totally cheating, but it works well, and it is
     * no more dumb than the hwloc interface itself.
     */
    operator hwloc_topology_t() const { return x; }
    hwloc_topology_t operator->() const { return x; }

    int get_nbobjs_by_type(hwloc_obj_type_t type) const {
        return hwloc_get_nbobjs_by_type(x_nonconst(), type);
    }
    unsigned int get_nbobjs_by_depth(int depth) const {
        return hwloc_get_nbobjs_by_depth(x_nonconst(), depth);
    }
    int get_type_depth(hwloc_obj_type_t type) const {
        return hwloc_get_type_depth(x_nonconst(), type);
    }
    hwloc_obj_t get_obj_by_depth(int depth, unsigned int idx) const {
        return hwloc_get_obj_by_depth(x_nonconst(), depth, idx);
    }
};

struct cxx_hwloc_bitmap {
    hwloc_bitmap_t x;
    cxx_hwloc_bitmap()
        : x(hwloc_bitmap_alloc())
    {}
    ~cxx_hwloc_bitmap() { hwloc_bitmap_free(x); }
    cxx_hwloc_bitmap(cxx_hwloc_bitmap && o)
        : x(o.x)
    {
        o.x = nullptr;
    }
    cxx_hwloc_bitmap& operator=(cxx_hwloc_bitmap && o) noexcept {
        hwloc_bitmap_t y = x;
        x = o.x;
        o.x = y;
        return *this;
    }
    cxx_hwloc_bitmap(cxx_hwloc_bitmap const & o)
        : x(hwloc_bitmap_dup(o.x))
    { }
    cxx_hwloc_bitmap & operator=(cxx_hwloc_bitmap const & o) noexcept {
        hwloc_bitmap_copy(x, o.x);
        return *this;
    }
    explicit cxx_hwloc_bitmap (hwloc_const_bitmap_t ox)
        : x(hwloc_bitmap_dup(ox))
    { }
    cxx_hwloc_bitmap & operator=(hwloc_const_bitmap_t ox) {
        hwloc_bitmap_copy(x, ox);
        return *this;
    }
    operator hwloc_bitmap_t() { return x; }
    operator hwloc_const_bitmap_t() const { return x; }
    hwloc_bitmap_t operator->() { return x; }
    hwloc_const_bitmap_t operator->() const { return x; }
    /* We include these for convenience and illustration, but really, as
     * we always do with these wrapper types, we advocate the use of the
     * C bindings instead. We don't want to replicate the full api here.
     */
    cxx_hwloc_bitmap operator|(cxx_hwloc_bitmap const & o) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_or(res.x, x, o.x);
        return res;
    }
    cxx_hwloc_bitmap operator|(hwloc_bitmap_t ox) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_or(res.x, x, ox);
        return res;
    }
    cxx_hwloc_bitmap operator&(cxx_hwloc_bitmap const & o) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_and(res.x, x, o.x);
        return res;
    }
    cxx_hwloc_bitmap operator&(hwloc_bitmap_t ox) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_and(res.x, x, ox);
        return res;
    }
    cxx_hwloc_bitmap operator^(cxx_hwloc_bitmap const & o) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_xor(res.x, x, o.x);
        return res;
    }
    cxx_hwloc_bitmap operator^(hwloc_bitmap_t ox) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_xor(res.x, x, ox);
        return res;
    }
    cxx_hwloc_bitmap operator!() const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_not(res.x, x);
        return res;
    }
    cxx_hwloc_bitmap& operator|=(cxx_hwloc_bitmap const & o) {
        hwloc_bitmap_or(x, x, o.x);
        return *this;
    }
    cxx_hwloc_bitmap& operator|=(hwloc_bitmap_t ox) {
        hwloc_bitmap_or(x, x, ox);
        return *this;
    }
    cxx_hwloc_bitmap& operator&=(cxx_hwloc_bitmap const & o) {
        hwloc_bitmap_and(x, x, o.x);
        return *this;
    }
    cxx_hwloc_bitmap& operator&=(hwloc_bitmap_t ox) {
        hwloc_bitmap_and(x, x, ox);
        return *this;
    }
    cxx_hwloc_bitmap& operator^=(cxx_hwloc_bitmap const & o) {
        hwloc_bitmap_xor(x, x, o.x);
        return *this;
    }
    cxx_hwloc_bitmap& operator^=(hwloc_bitmap_t ox) {
        hwloc_bitmap_xor(x, x, ox);
        return *this;
    }
    auto operator<=>(cxx_hwloc_bitmap const & o) const {
        return hwloc_bitmap_compare(x, o.x) <=> 0;
    }
    auto operator==(cxx_hwloc_bitmap const & o) const {
        return operator<=>(o) == 0;
    }
};
using cxx_hwloc_cpuset = cxx_hwloc_bitmap;
using cxx_hwloc_nodeset = cxx_hwloc_bitmap;

namespace fmt {
    template <> struct formatter<cxx_hwloc_bitmap> : formatter<string_view> {
        auto format(cxx_hwloc_bitmap const & e, format_context& ctx) const -> format_context::iterator
        {
            char * s = nullptr;
            int const rc = hwloc_bitmap_asprintf(&s, e);
            ASSERT_ALWAYS(rc >= 0);
            format_to(ctx.out(), "{}", s);
            free(s);
            return ctx.out();
        }
    };
} /* namespace fmt */

#endif	/* CADO_HWLOC_AUX_HPP */
