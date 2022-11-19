#ifndef HWLOC_AUX_H_
#define HWLOC_AUX_H_

#include <hwloc.h>

#ifdef __cplusplus
extern "C" {
#endif


int hwloc_aux_get_depth_from_string(hwloc_topology_t topology, const char * desc);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
struct cxx_hwloc_bitmap {
    hwloc_bitmap_t x;
    cxx_hwloc_bitmap() { x = hwloc_bitmap_alloc(); }
    ~cxx_hwloc_bitmap() { hwloc_bitmap_free(x); }
    cxx_hwloc_bitmap(hwloc_const_bitmap_t o) {
        x = hwloc_bitmap_dup(o);
    }
    cxx_hwloc_bitmap(cxx_hwloc_bitmap && o) {
        x = o.x;
        o.x = NULL;
    }
    cxx_hwloc_bitmap& operator=(cxx_hwloc_bitmap && o) {
        hwloc_bitmap_t y = x;
        x = o.x;
        o.x = y;
        return *this;
    }
    cxx_hwloc_bitmap(cxx_hwloc_bitmap const & o) {
        x = hwloc_bitmap_dup(o.x);
    }
    cxx_hwloc_bitmap & operator=(cxx_hwloc_bitmap const & o) {
        hwloc_bitmap_copy(x, o.x);
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
    cxx_hwloc_bitmap operator&(cxx_hwloc_bitmap const & o) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_and(res.x, x, o.x);
        return res;
    }
    cxx_hwloc_bitmap operator^(cxx_hwloc_bitmap const & o) const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_xor(res.x, x, o.x);
        return res;
    }
    cxx_hwloc_bitmap operator!() const {
        cxx_hwloc_bitmap res;
        hwloc_bitmap_not(res.x, x);
        return res;
    }
    bool operator<(cxx_hwloc_bitmap const & o) const {
        return hwloc_bitmap_compare(x, o.x) < 0;
    }
};
typedef cxx_hwloc_bitmap cxx_hwloc_cpuset;
typedef cxx_hwloc_bitmap cxx_hwloc_nodeset;
#endif

#endif	/* HWLOC_AUX_H_ */
