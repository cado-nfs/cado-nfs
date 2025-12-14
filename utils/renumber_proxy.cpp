#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstddef>

#include <utility>       // for pair
#include <vector>        // for vector

#include "cado_poly.h"
#include "mpz_poly.h"
#include "renumber.hpp"
#include "renumber_proxy.h"
#include "typedefs.h"

static inline const renumber_t * deref(renumber_proxy_srcptr R)
{
    return (const renumber_t *) R->x;
}

static inline renumber_t * deref(renumber_proxy_ptr R)
{
    return (renumber_t *) R->x;
}

void renumber_table_init(renumber_proxy_ptr R, cado_poly_ptr P)
{
    cxx_cado_poly const ff(P);
    R->x = (void*) new renumber_t(ff);
}

void renumber_table_clear(renumber_proxy_ptr R)
{
    delete deref(R);
}

void renumber_table_set_lpb(renumber_proxy_ptr R, const unsigned int * lpb, size_t n)
{
    std::vector<unsigned int> const L(lpb, lpb+n);
    deref(R)->set_lpb(L);
}

void renumber_table_read_from_file(renumber_proxy_ptr R, const char * filename, int for_dl)
{
    deref(R)->read_from_file(filename, for_dl);
}

int renumber_table_get_format(renumber_proxy_srcptr R)
{
    return deref(R)->get_format();
}

unsigned int renumber_table_get_lpb(renumber_proxy_srcptr R, int side)
{
    return deref(R)->get_lpb(side);
}
unsigned int renumber_table_get_max_lpb(renumber_proxy_srcptr R)
{
    return deref(R)->get_max_lpb();
}
unsigned int renumber_table_get_min_lpb(renumber_proxy_srcptr R)
{
    return deref(R)->get_min_lpb();
}
uint64_t renumber_table_get_size(renumber_proxy_srcptr R)
{
    return deref(R)->size();
}
unsigned int renumber_table_get_nb_polys(renumber_proxy_srcptr R)
{
    return deref(R)->get_nb_polys();
}
mpz_poly_srcptr renumber_table_get_poly(renumber_proxy_srcptr R, int side)
{
    return deref(R)->get_poly(side);
}
int renumber_table_get_poly_deg(renumber_proxy_srcptr R, int side)
{
    return deref(R)->get_poly_deg(side);
}
int renumber_table_get_rational_side(renumber_proxy_srcptr R)
{
    return deref(R)->get_rational_side();
}
index_t renumber_table_get_max_index(renumber_proxy_srcptr R)
{
    return deref(R)->get_max_index();
}
index_t renumber_table_get_max_cached_index(renumber_proxy_srcptr R)
{
    return deref(R)->get_max_cached_index();
}
index_t number_of_additional_columns(renumber_proxy_srcptr R)
{
    return deref(R)->number_of_additional_columns();
}

size_t renumber_table_get_sides_of_additional_columns(renumber_proxy_srcptr R, int * sides, size_t * n)
{
    auto v = deref(R)->get_sides_of_additional_columns();
    size_t i = 0;
    for( ; i < v.size() && i < *n ; ++i) {
        sides[i] = v[i];
    }
    *n = i;
    return v.size();
}
index_t number_of_bad_ideals(renumber_proxy_srcptr R)
{
    return deref(R)->number_of_bad_ideals();
}
size_t renumber_table_get_memory_size(renumber_proxy_srcptr R)
{
    return deref(R)->get_memory_size();
}

int renumber_table_index_is_bad(renumber_proxy_srcptr R, index_t * first_index, index_t h)
{
    index_t f;
    int const n = deref(R)->is_bad(f, h);
    if (!n) return n;
    if (first_index) *first_index = f;
    return n;
}

int renumber_table_p_r_side_is_bad(renumber_proxy_srcptr R, index_t * first_index, p_r_values_t p, p_r_values_t r, int side)
{
    index_t f;
    int const n = deref(R)->is_bad(f, p, r, side);
    if (!n) return n;
    if (first_index) *first_index = f;
    return n;
}

bool renumber_table_index_is_additional_column(renumber_proxy_srcptr R, index_t h)
{
    return deref(R)->is_additional_column(h);
}

index_t renumber_table_index_from_p_r (renumber_proxy_srcptr R, p_r_values_t p, p_r_values_t r, int side)
{
    return deref(R)->index_from_p_r(p, r, side);
}

int renumber_table_p_r_side_get_inertia (renumber_proxy_srcptr R, p_r_values_t p, p_r_values_t r, int side)
{
    return deref(R)->inertia_from_p_r(p, r, side);
}

bool renumber_table_p_r_from_index(renumber_proxy_srcptr R, p_r_values_t * pp, p_r_values_t * pr, int * pside, index_t h)
{
    try {
        renumber_t::p_r_side const x = deref(R)->p_r_from_index(h);
        if (pp) *pp = x.p;
        if (pr) *pr = x.r;
        if (pside) *pside = x.side;
        return true;
    } catch (renumber_t::corrupted_table const &) {
        return false;
    }
}

bool renumber_table_indices_from_p_a_b(renumber_proxy_srcptr R, index_t * first, int * exps, size_t * nexps, p_r_values_t p, p_r_values_t r, int side, int e, int64_t a, uint64_t b)
{
    try {
        renumber_t::p_r_side const x { p, r, side };
        auto v = deref(R)->indices_from_p_a_b(x, e, a, b);
        if (first) *first = v.first;
        if (exps && nexps) {
            size_t i = 0;
            for( ; i < *nexps && i < v.second.size() ; i++) {
                exps[i] = v.second[i];
            }
            *nexps = i;
        }
        return true;
    } catch (renumber_t::corrupted_table const &) {
        return false;
    }
}
