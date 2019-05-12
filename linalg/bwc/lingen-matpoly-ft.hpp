#ifndef LINGEN_MATPOLY_FT_H_
#define LINGEN_MATPOLY_FT_H_

#include "lingen-matpoly.hpp"
#include "flint-fft/fft.h"
#include "lingen-substep-schedule.h"

struct matpoly_ft {
    abdst_field ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
    const struct fft_transform_info * fti = NULL;
    void * data = NULL;
    int check_pre_init() const { return data == NULL; }
    matpoly_ft(abdst_field ab, unsigned int m, unsigned int n, const struct fft_transform_info * fti);
    matpoly_ft() = default;
    ~matpoly_ft();
    void zero();
    void dft(matpoly const & p);
    void add(matpoly_ft const & t0, matpoly_ft const & t1);
    void sub(matpoly_ft const & t0, matpoly_ft const & t1);
    void mul(matpoly_ft const & t0, matpoly_ft const & t1);
    void addmul(matpoly_ft const & t0, matpoly_ft const & t1);
    void ift(matpoly & p);
    void ift_mp(matpoly & p, unsigned int shift);
    void to_export();
    void to_import();
};

/* In a way, this is the only real API exported by this module */
void matpoly_mul_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void matpoly_mul_caching(matpoly & c, matpoly const & a, matpoly const & b, const struct lingen_substep_schedule * S) { return matpoly_mul_caching_adj(c, a, b, UINT_MAX, S); }

void matpoly_mp_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void matpoly_mp_caching(matpoly & c, matpoly const & a, matpoly const & b, const struct lingen_substep_schedule * S) { return matpoly_mp_caching_adj(c, a, b, UINT_MAX, S); }

#endif	/* LINGEN_MATPOLY_FT_H_ */
