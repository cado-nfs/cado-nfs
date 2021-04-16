#ifndef FB_TYPES_H
#define FB_TYPES_H

// pragma no prototypes

/* Elementary data types for the factor base */

#include <stdint.h>
#include <inttypes.h>
#include "las-config.h"

typedef uint32_t fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT PRIu32
#define FBPRIME_MAX UINT32_MAX
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT PRIu32
typedef fbprime_t redc_invp_t;



/* Within one factor base, there is exactly one (index, offset) tuple per
   factor base entry. */
/* Each slice in a factor base has a unique index */
typedef uint32_t slice_index_t;
/* Each factor base entry withing a slice has a unique offset */
typedef uint16_t slice_offset_t;

// FIXME: could probably go somewhere else...
// Small struct for sublattice info:
// One sieves only positions congrent to (i0,j0) mod m
struct sublat_s {
    uint32_t m=0; // 0 means no sublattices.
    uint32_t i0=0;
    uint32_t j0=0;
};
typedef struct sublat_s sublat_t;

#ifdef __cplusplus
/* This structure encodes an element of P1(Z/p^k) for some p^k (which is
 * not embedded in the structure).
 *
 * If proj == false, r encodes the element (r:1)
 *      in this case, we have 0<=r<p^k
 *
 * If proj == true, r encodes the element (1:r)
 *      in this case, we have p|r and 0<=r/p<p^{k-1}
 *
 */

template<typename T>
struct fb_root_p1_t {
    T r;
    bool proj;
    inline bool is_affine() const { return !proj; }
    inline bool is_projective() const { return proj; }
    fb_root_p1_t(T const & r, bool proj = false) : r(r), proj(proj) {}
    fb_root_p1_t(fb_root_p1_t const &) = default;
    fb_root_p1_t(fb_root_p1_t &&) = default;
    fb_root_p1_t() = default;
    fb_root_p1_t& operator=(fb_root_p1_t const &) = default;
    fb_root_p1_t& operator=(fb_root_p1_t &&) = default;
    static fb_root_p1_t affine_root(fbprime_t r) { return r; }
    static fb_root_p1_t projective_root(fbprime_t r) { return fb_root_p1_t { r, true }; }
    inline fbprime_t to_old_format(fbprime_t p) const { return r + (proj ? p : 0); }
    static fb_root_p1_t<T> from_old_format(fbprime_t old_r, fbprime_t q) {
        return fb_root_p1_t<T>((old_r >= q) ? (old_r - q) : old_r, old_r >= q);;
    }
};
template<typename T, typename U>
inline bool operator==(fb_root_p1_t<T> const & a, fb_root_p1_t<U> const & b)
{
    return a.r == b.r && a.proj == b.proj;
}
template<typename T, typename U>
inline bool operator!=(fb_root_p1_t<T> const & a, fb_root_p1_t<U> const & b)
{
    return !(a == b);
}

/* This is the type that we'll be using most of the time */
typedef fb_root_p1_t<fbprime_t> fb_root_p1;

#endif


#endif
