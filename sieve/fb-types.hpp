#ifndef CADO_FB_TYPES_HPP
#define CADO_FB_TYPES_HPP

// pragma no prototypes

/* Elementary data types for the factor base */

#include <cstdint>
#include <cinttypes>

using fbprime_t = uint32_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT PRIu32
#define FBPRIME_MAX UINT32_MAX
#define FBPRIME_BITS 32
using fbroot_t = fbprime_t;
#define FBROOT_FORMAT PRIu32
using redc_invp_t = fbprime_t;



/* Within one factor base, there is exactly one (index, offset) tuple per
   factor base entry. */
/* Each slice in a factor base has a unique index */
using slice_index_t = uint32_t;
/* Each factor base entry within a slice has a unique offset */
using slice_offset_t = uint16_t;

// FIXME: could probably go somewhere else...
// Small struct for sublattice info:
// One sieves only positions congruent to (i0,j0) mod m
struct sublat_t {
    uint32_t m=0; // 0 means no sublattices.
    uint32_t i0=0;
    uint32_t j0=0;

    void adjustIJ(int & i, unsigned int & j) const
    {
        if (m != 0) {
            i = i*m + i0;
            j = j*m + j0;
        }
    }
};

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
    bool proj = false;
    bool is_affine() const { return !proj; }
    bool is_projective() const { return proj; }
    fb_root_p1_t(T const & r, bool proj = false) : r(r), proj(proj) {}
    fb_root_p1_t() = default;
    /*
    fb_root_p1_t(fb_root_p1_t const &) = default;
    fb_root_p1_t(fb_root_p1_t &&) = default;
    ~fb_root_p1_t() = default;
    fb_root_p1_t& operator=(fb_root_p1_t const &) = default;
    fb_root_p1_t& operator=(fb_root_p1_t &&) = default;
    */
    static fb_root_p1_t affine_root(fbprime_t r) { return r; }
    static fb_root_p1_t projective_root(fbprime_t r) { return fb_root_p1_t { r, true }; }
    fbprime_t to_old_format(fbprime_t p) const { return r + (proj ? p : 0); }
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
using fb_root_p1 = fb_root_p1_t<fbprime_t>;

#endif


#endif
