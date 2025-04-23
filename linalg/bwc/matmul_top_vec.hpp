#ifndef CADO_MATMUL_TOP_VEC_HPP
#define CADO_MATMUL_TOP_VEC_HPP

#include <cstdint>

#include <string>
#include <memory>

#include "gmp_aux.h"
#include "parallelizing_info.hpp"
#include "arith-generic.hpp"
#include "macros.h"

struct matmul_top_data;

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_vec {
    arith_generic * abase = nullptr;
    parallelizing_info_ptr pi = nullptr;
    int d = -1;
    pi_datatype_ptr pitype = nullptr;

    unsigned int n = 0;          // total size in items
    unsigned int i0 = 0;
    unsigned int i1 = 0;

    arith_generic::owned_vector owned_v;

    // v == owned_v except for shared vector, where this holds only at
    // pi->wr[d]->trank == 0
    arith_generic::elt * v = nullptr;

    /* pointers to other vectors are held in shared memory areas (shared
     * among cores, of course)
     */
    typedef std::unique_ptr<
                mmt_vec *,
                shared_free_deleter<mmt_vec *>> pointer_to_others;

    /* pi->wr[d]->ncores siblings ; only in case all cores in the
     * communicator have their own data area v */
    pointer_to_others siblings;

    /* pi->m->ncores siblings, always */
    pointer_to_others mpals;

    /* pi->wr[0]->ncores and */
    /* pi->wr[1]->ncores siblings, always. */
    pointer_to_others wrpals[2];

    private:
    void set_pointer_to_others(bool with_siblings);

    public:

    mmt_vec() = default;
    mmt_vec(matmul_top_data & mmt, arith_generic * abase, pi_datatype_ptr pitype, int d, int flags, unsigned int n);
    mmt_vec(mmt_vec const &) = delete;
    mmt_vec& operator=(mmt_vec const &) = delete;
    ~mmt_vec() = default;

    mmt_vec& operator=(mmt_vec &&) = delete;

    /* We need to do some stuff here! */
    mmt_vec(mmt_vec && o) noexcept;

    size_t rsbuf_items = 0;      // auto-expanded on demand.

    // ATTENTION: elements are laid out flat in here! rsbuf[0][i] is not
    // a thing.
    arith_generic::owned_vector rsbuf[2];            // only for RS_CHOICE == RS_CHOICE_MINE
    int consistency = 0;    /* 0 == inconsistent ; 1 == partial ; 2 == full */

    // XXX private ?
    mmt_vec & sibling(unsigned int i) { return siblings ? *(siblings.get()[i]) : *this; }
    mmt_vec const & sibling(unsigned int i) const { return siblings ? *(siblings.get()[i]) : *this; }
};

/* some handy macros to access portions of mmt_vec's */
#define SUBVEC(v,w,offset) (v)->abase->vec_subvec((v)->abase, (v)->w, offset)
#define SUBVEC_const(v,w,offset) (v)->abase->vec_subvec_const((v)->abase, (v)->w, offset)

/* These are flags for the distributed vectors. For the moment we have
 * only one flag */
#define THREAD_SHARED_VECTOR    1

extern size_t mmt_my_own_size_in_items(mmt_vec const & v);
extern size_t mmt_my_own_size_in_bytes(mmt_vec const & v);
extern size_t mmt_my_own_offset_in_items(mmt_vec const & v);
extern size_t mmt_my_own_offset_in_items(mmt_vec const & v, unsigned int);
extern size_t mmt_my_own_offset_in_bytes(mmt_vec const & v);
extern arith_generic::elt * mmt_my_own_subvec(mmt_vec & v);
extern arith_generic::elt * mmt_my_own_subvec(mmt_vec & v, unsigned int);
extern arith_generic::elt const * mmt_my_own_subvec(mmt_vec const & v);

/* do not use this function if you want consistency when the splitting
 * changes ! */
extern void mmt_vec_set_random_inconsistent(mmt_vec & v, cxx_gmp_randstate & rstate);
extern unsigned long mmt_vec_hamming_weight(mmt_vec const & y);
extern void mmt_vec_set_x_indices(mmt_vec & y, std::vector<uint32_t> const & gxvecs, unsigned int j0, unsigned int j1, unsigned int nx);
extern void mmt_vec_set_expanded_copy_of_local_data(mmt_vec & y, const void * v, unsigned int n);

extern void mmt_own_vec_set(mmt_vec & w, mmt_vec const & v);
extern void mmt_own_vec_set2(mmt_vec const & z, mmt_vec & w, mmt_vec const & v);
// extern void mmt_vec_swap(mmt_vec & w, mmt_vec & v);
extern void mmt_full_vec_set(mmt_vec & w, mmt_vec const & v);
extern void mmt_full_vec_set_zero(mmt_vec & v);
extern void mmt_vec_set_basis_vector_at(mmt_vec & v, int k, unsigned int j);
extern void mmt_vec_set_basis_vector(mmt_vec & v, unsigned int j);
extern void mmt_vec_add_basis_vector_at(mmt_vec & v, int k, unsigned int j);
extern void mmt_vec_add_basis_vector(mmt_vec & v, unsigned int j);
#if 0
extern void matmul_top_fill_random_source_generic(matmul_top_data & mmt, size_t stride, mmt_vec & v, int d);
#endif
extern int mmt_vec_load(mmt_vec & v, std::string const & name, unsigned int itemsondisk, unsigned int block_position) ATTRIBUTE_WARN_UNUSED_RESULT;
extern int mmt_vec_save(mmt_vec & v, std::string const & name, unsigned int itemsondisk, unsigned int block_position);
extern void mmt_vec_clear_padding(mmt_vec & v, size_t unpadded, size_t padded);

static inline int mmt_vec_is_shared(mmt_vec const & v) {
    return v.siblings.get() == nullptr;
}

extern void mmt_vec_share_across_threads(mmt_vec & v);

extern void mmt_vec_set_random_through_file(mmt_vec & v, std::string const & name, unsigned int itemsondisk, cxx_gmp_randstate & rstate, unsigned int block_position);

// private ?
extern void mmt_vec_downgrade_consistency(mmt_vec & v);

#endif	/* MATMUL_TOP_VEC_HPP_ */
