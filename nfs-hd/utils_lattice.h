#ifndef UTILS_SIEVE_H
#define UTILS_SIEVE_H 

#include <stdint.h>
#include "utils.h"
#include "int64_vector.h"
#include "double_vector.h"
#include "mat_int64.h"
#include "array.h"
#include "ideal.h"
#include "mat_double.h"
#include "list_int64_vector.h"
#include "list_double_vector.h"
#include "list_int64_vector_index.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Compute the acute Gauss reduction of two vector (v0_root and v1_root) if they
 *  are in the same plane. Return the determinant of U if U is not set as NULL,
 *  0 otherwise.
 *
 * v0: the first output vector.
 * v1: the second output vector.
 * U: unimodular matrix to go from the two first coordinate of v0_root and
 *  v1_root to the two first coordinate of v0 and v1. (v0_root, v1_root) * U =
 *  (v0, v1).
 * v0_root: the first input vector.
 * v1_root: the second input vector.
 */
int gauss_reduction(int64_vector_ptr v0, int64_vector_ptr v1, mat_int64_ptr U,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

/*
 * Do the same thing as gauss_reduction, but set all of the other two first
 * coordinates of v0 and v1 to zero.
 */
int gauss_reduction_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

/*
 * Skew LLL (vector in columns).
 * Work for square matrix.
 *
 * MSLLL: the ouput matrix.
 * M_root: the initial matrix.
 * skewness: skewness (for the rows).
 */
void skew_LLL(mat_int64_ptr MSLLL, mat_int64_srcptr M_root,
    int64_vector_srcptr skewness, FILE * errstd);

/*
 * Use the Franke-Kleinjung to modify the two first coordinates of v0_root and
 *  v1_root, if they respect the conditions. See sieve/las-plattice.h
 *  (03/31/2015). Return 1 if the reduction is done, 0 otherwise.
 *
 * v0: the output first vector.
 * v1: the output second vector.
 * v0_root: the initial first vector.
 * v1_root: the initial second vector
 * I: minimal gap between the two x coordinate of v0 and v1. Usually, I = 2*H0.
 */
int reduce_qlattice(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I);

/*
 * Same thing as reduce_qlattice above, but set the other two first coordinates
 *  to zero for v0 and v1.
 */
int reduce_qlattice_zero(int64_vector_ptr v0, int64_vector_ptr v1,
                         int64_vector_srcptr v0_root,
                         int64_vector_srcptr v1_root, int64_t I);

/*
 * The SV4 (for instance, is just SV3 but I keep the name for instance) to find
 *  some short vector in the lattice defined by (v0_root, v1_root, v2)  with a
 *  z coordinate equal to 1. Very dependent on the form of the vector: v0 = (a,
 *  b, 0, …, 0), v1 = (c, d, 0, …, 0) and v2 = (e, 0, 1, 0, …, 0).
 *
 * SV: the list of short vector close to v2.
 * v0_root, v1_root, v2: three vectors that gives a basis of a lattice.
 */
void SV4(list_int64_vector_ptr SV, int64_vector_srcptr v0_root,
         int64_vector_srcptr v1_root, int64_vector_srcptr v2);

/*
 * Same thing as SV4, but the basis vector is in a Mqr matrix (vector in
 * columns).
 */
void SV4_Mqr(list_int64_vector_ptr SV, mat_int64_srcptr Mqr);

/*
 * Add an FK vector (e0 or e1) if v is outside of the sieving region defined by
 * H to have the x coordinate of v in [-H0, H0[.
 *
 * v: current vector.
 * e0: a vector given by the Franke-Kleinjung algorithm.
 * e1: a vector given by the Franke-Kleinjung algorithm.
 * H: sieving bound.
 */
void add_FK_vector(int64_vector_ptr v, list_int64_vector_srcptr list,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H);

/*
 * Enumerate with the Franke-Kleinjung algorithm for x increasing. v->c[0]
 *  is in [A, A + I[. Return 0 if v0 is added, 1 if v1 is added, 2 if v0 + v1 is added.
 *
 * v: v_old + l v0 + m v1, (l, m) in [0, 1]^2. v->c[0] > v_old->c[0].
 * v_old: the starting point vector.
 * v0: a vector given by reduce_qlattice, v0->c[0] < 0.
 * v1: a vector given by reduce_qlattice, v1->c[0] > 0.
 * A: how the interval is centered.
 * I: the width of the interval.
 */
unsigned int enum_pos_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I);

/*
 * Enumerate with the same technique as in Franke-Kleinjung algorithm but for x
 *  decreasing. v->c[0] is in ]-A - I, -A]. Return 0 if v0 is added, 1 if v1 is added,
 *  2 if v0 + v1 is added.
 *
 * v: v_old + l v0 + m v1, (l, m) in [0, 1]^2. v->c[0] < v_old->c[0].
 * v_old: the starting point vector.
 * v0: a vector given by reduce_qlattice, v0->c[0] < 0.
 * v1: a vector given by reduce_qlattice, v1->c[0] > 0.
 * A: how the interval is centered.
 * I: the width of the interval.
 */
unsigned int enum_neg_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I);

/*
 * Compute the contribution of adding v0 or v1 to the index of the array in
 *  which we store the norm to go from v to v + l * v0 + m * v1.
 *
 * coord_v0: the contribution of v0.
 * coord_v1: the contribution of v1.
 * v0: the vector v0 given by the Franke-Kleinjung algorithm.
 * v1: the vector v0 given by the Franke-Kleinjung algorithm.
 * H: the sieving bound.
 * number_element: number of element in the array.
 */
void coordinate_FK_vector(uint64_t * coord_v0, uint64_t * coord_v1,
    int64_vector_srcptr v0, int64_vector_srcptr v1, sieving_bound_srcptr H,
    uint64_t number_element);

/*
 * Compute the contribution of adding a vector v in array_ptr.
 */
uint64_t index_vector(int64_vector_srcptr v, sieving_bound_srcptr H,
    uint64_t number_element);

/*
 * Go from the plane z=d to the plane z=d+1 and keep the x coordinate of vs in
 *  [-H0, H0[.
 *
 * vs: the starting point vector in the plane z=d (at beginning) and z=d+1 (at
 *  the end).
 * SV: the possible short vector to go from z=d to z=d+1.
 * e0: a Franke-Kleinjung vector.
 * e1: a Franke-Kleinjung vector.
 * H: the sieving bound.
 */
void plane_sieve_next_plane(int64_vector_ptr vs, list_int64_vector_srcptr SV,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H,
    int boolean);

/*
 * Compute the Gram-Schmidt orthogonalisation of the vector in list_old to
 *  produce list_new and the m matrix, the matrix with the result of the
 *  orthogonal projection (see page 164 of "Algorithms for the Shortest and
 *  Closest Lattice Vector Problems" by Guillaume Hanrot, Xavier Pujol, and
 *  Damien Stehlé, IWCC 2011).
 *
 * list_new: the Gram-Schmidt basis.
 * m: the matrix with the coefficient \mu_{i, j}.
 * list_old: the original basis.
 */
void double_vector_gram_schmidt(list_double_vector_ptr list_new,
    mat_double_ptr m, list_double_vector_srcptr list_old);

/*
 * TODO: is s_out really needed? Can be s?
 * When the space sieve is not able to go from a plane to the next plane, call
 *  plane_sieve_1_incomplete to find the next reached plane.
 *
 * s_out: new starting point vector.
 * s: old starting point vector.
 * Mqr: the Mqr matrix.
 * H: sieving bound.
 * list_FK: the 2 Franke-Kleinjung vectors.
 * list_SV: the list of short vectors with a z coordinate equal to 1.
 */
void plane_sieve_1_incomplete(int64_vector_ptr s_out, int64_vector_srcptr s,
    MAYBE_UNUSED mat_int64_srcptr Mqr, sieving_bound_srcptr H,
    list_int64_vector_srcptr list_FK, list_int64_vector_srcptr list_SV);

/*
 * Beginning of the space sieve. ie, store in list_vec the 3 vectors given by
 *  skew LLL on Mqr and return 1 if a vector with a z coordinate equal to 1 is
 *  in list.
 */
unsigned int space_sieve_1_init(list_int64_vector_index_ptr list_vec,
    list_int64_vector_index_ptr list_vec_zero, ideal_1_srcptr r,
    mat_int64_srcptr Mqr, sieving_bound_srcptr H, uint64_t number_element,
    unsigned int * skew_lll_fail, FILE * file_skew_lll_stat, FILE * errstd);
/*
 * During the plane sieve, compute reduce q lattice and SV4 and store vectors in
 *  the appropriate lists.
 */
int space_sieve_1_plane_sieve_init(list_int64_vector_ptr list_SV,
    list_int64_vector_ptr list_FK, list_int64_vector_index_ptr list_vec,
    list_int64_vector_index_ptr list_vec_zero, MAYBE_UNUSED ideal_1_srcptr r,
    sieving_bound_srcptr H, mat_int64_srcptr Mqr,
    unsigned int vector_1, uint64_t number_element, unsigned int * new_vec);

/*
 * Store in s_tmp the starting point with the smallest z coordinate that are in
 *  the sieving region.
 */
unsigned int space_sieve_1_next_plane_seek(int64_vector_ptr s_tmp,
    unsigned int * index_vec, unsigned int * s_change,
    list_int64_vector_srcptr list_s, list_int64_vector_index_srcptr list_vec,
    sieving_bound_srcptr H, int64_vector_srcptr s);

void lll_Mqr(mat_int64_ptr C, mat_int64_srcptr A, FILE * errstd);

void slll_Mqr(mat_int64_ptr C, mat_int64_srcptr A, int64_vector_srcptr,
    FILE * errstd);

void lll_Mqr_unimodular(mat_int64_ptr U, mat_int64_srcptr A, FILE * errstd);

#ifdef __cplusplus
}
#endif
#endif // UTILS_SIEVE_H
