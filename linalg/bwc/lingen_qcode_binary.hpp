#ifndef LINGEN_QCODE_BINARY_HPP_
#define LINGEN_QCODE_BINARY_HPP_

#include <cstddef>
#include <gmp.h>
#include "lingen_bmstatus.hpp"
#include "mpfq_fake.hpp"
#include "lingen_matpoly_binary.hpp"

/* We have two interfaces here. The first one is the one that is common
 * with qcode_prime, and that we want the future lingen program to use.
 *
 * The second one, which is activated if
 * LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE is #defined prior to
 * #including this file, is really legacy code, and it's only exposed for
 * the legacy lingen_binary code. (as a matter of fact, the
 * implementation of the future interface does build upon the legacy
 * interface presently, but that is not a reason to have it exposed).
 */
extern int
bw_lingen_basecase(bmstatus & bm, matpoly & pi, matpoly & E);
extern void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);

#ifdef LINGEN_QCODE_BINARY_TRAMPOLINE_INTERFACE
/* This trampoline structure is no longer useful, really. At some point
 * we had both C and C++ 
 */
struct lingen_qcode_data_s {
    /* we have a matrix of size m*b, and one of size b*b. The second
     * dimension is not called n in order to avoid confusion with the n
     * in BW ; in fact we have b = m + n */
    unsigned int m, b;
    unsigned int t;
    unsigned long length, outlength;
    unsigned int luck_mini;

    unsigned int * local_delta;

    /* where we grab our input and store our output */
    /* Note that we don't own the data corresponding to the innermost
     * level (while we do own the outermost level for iptrs and optrs).
     */
    unsigned int * delta;
    int * ch;

    const unsigned long ** iptrs;
    unsigned long ** optrs;
};


typedef struct lingen_qcode_data_s lingen_qcode_data[1];
typedef struct lingen_qcode_data_s * lingen_qcode_data_ptr;
typedef const struct lingen_qcode_data_s * lingen_qcode_data_srcptr;

void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength);
void lingen_qcode_clear(lingen_qcode_data_ptr qq);

unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j);

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq);

static inline void lingen_qcode_hook_delta(lingen_qcode_data_ptr qq, unsigned int * delta)
{
    qq->delta = delta;
}

static inline void lingen_qcode_hook_chance_list(lingen_qcode_data_ptr qq, int * ch)
{
    qq->ch = ch;
}

static inline void lingen_qcode_hook_input(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->iptrs[i * qq->b + j] = poly;
}

static inline void lingen_qcode_hook_output(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->optrs[i * qq->b + j] = poly;
}
#endif

#endif	/* LINGEN_QCODE_BINARY_HPP_ */
