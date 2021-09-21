#ifndef PLATTICE_PROXY_HPP_
#define PLATTICE_PROXY_HPP_

#include <map>
#include "las-plattice.hpp"

struct plattice_proxy : public plattice_info {
    /* use this default ctor only for the comparison with the reference
     * routines.
     */
    plattice_proxy () = default;

    using plattice_info::mi0;
    using plattice_info::j0;
    using plattice_info::i1;
    using plattice_info::j1;
    using plattice_info::check_post_conditions;

#ifndef NDEBUG
#define ASSERT_THROW(e, c) do { if (!(c)) throw (e)(#c); } while (0)
#else
#define ASSERT_THROW(e, c)
#endif
    struct error : public std::runtime_error {
        error(const char * s) : std::runtime_error(s) {}
    };
#define ASSERT_PLATTICE(c) ASSERT_THROW(error, c)

#include "las-reduce-plattice-using-64bitmul.hpp"

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    /* The simplistic code is the fallback that gets used when we have no
     * inline assembly, so when we do have assembly, we want to include
     * it here */
#include "las-reduce-plattice-simplistic.hpp"
#endif

#include "las-reduce-plattice-swapping-loop.hpp"

#include "las-reduce-plattice-swapping-loop2.hpp"

#include "las-reduce-plattice-two-legs.hpp"

#include "las-reduce-plattice-mimick-production-noasm.hpp"


    void instrumented_two_legs(uint32_t I, std::map<int, unsigned long> & T) {
        /* This is the main reduce_plattice loop */
        for( ;; ) {
            if (i1 < I) {
                /* an "UNLIKELY" macro here actually has an adverse
                 * effect...  */
                if (i1 == 0) {
                    // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                    j0 = j1 - j0;
                    reduce_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / i1;
                mi0 -= a * i1;
                j0  += a * j1;
                return;
            }
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
                T[k]++;
            }
            if (mi0 < I) {
                /* an "UNLIKELY" macro here actually has an adverse
                 * effect...  */
                if (mi0 == 0) {
                    mi0 = i1;
                    i1 = j0 ; j0 = j1 ; j1 = i1;
                    i1 = 0;
                    reduce_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / mi0;
                i1 -= a * mi0;
                j1 += a * j0;
                return;
            }
            {
                int k = i1 / mi0; i1 -= k * mi0; j1 += k * j0;
                T[k]++;
            }
        }
    }

    /* XXX
     * beware: this constructor takes I, but it shadows a constructor in
     * the production code which takes only logI !!!
     */
    plattice_proxy(const unsigned long q, const unsigned long r, bool proj, uint32_t I) : plattice_info()
    {
        initial_basis(q, r, proj);
        /* At this point, (mi0,j0) represents itself, i.e. a vector with
         * two positive coordinates.
         * Note that j0==0
         */
        // ASSERT_ALWAYS(check_pre_conditions(I));
        bool needs_special_treatment = (i1 == 0 || (j1 > 1 && mi0 < I));
        if (needs_special_treatment) {
            reduce_with_vertical_vector(I);
            return;
        }
        reduce(I);
        // simplistic(I);
        // using_64bit_mul(I);
        // swapping_loop(I);
    }

    public:
    using plattice_info::initial_basis;
    using plattice_info::reduce;
#ifndef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using plattice_info::reduce_plattice_simplistic;
#endif
    using plattice_info::reduce_with_vertical_vector;
    using plattice_info::needs_special_treatment;

    // friend void instrumented_two_legs(plattice *pli, const unsigned long q, const unsigned long r, bool proj, uint32_t I, std::map<int, unsigned long> & T);
};


#endif	/* PLATTICE_PROXY_HPP_ */
