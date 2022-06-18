#include "cado.h"
#include "arith-generic.hpp"
#include "arith-mod2.hpp"
#include "arith-modp.hpp"
#include "bwc_config.h"
#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
#include <dlfcn.h>
#include "solib-naming.h"
#endif

template<typename T>
struct arith_wrapper: public arith_generic, public T {
    private:
        /* These functions are all no-ops */
        static typename T::elt * cast(arith_generic::elt * x) {
            return static_cast<typename T::elt *>(static_cast<arith_concrete_base::elt *>(x));
        }
        static typename T::elt const * cast(arith_generic::elt const * x) {
            return static_cast<typename T::elt const *>(static_cast<arith_concrete_base::elt const *>(x));
        }
        static typename T::elt & cast(arith_generic::elt & x) {
            return static_cast<typename T::elt &>(static_cast<arith_concrete_base::elt &>(x));
        }
        static typename T::elt const & cast(arith_generic::elt const & x) {
            return static_cast<typename T::elt const &>(static_cast<arith_concrete_base::elt const &>(x));
        }
        static arith_generic::elt * uncast(typename T::elt * x) {
            return static_cast<arith_generic::elt *>(static_cast<arith_concrete_base::elt *>(x));
        }
        static arith_generic::elt const * uncast(typename T::elt const * x) {
            return static_cast<arith_generic::elt const *>(static_cast<arith_concrete_base::elt const *>(x));
        }
        static arith_generic::elt & uncast(typename T::elt & x) {
            return static_cast<arith_generic::elt &>(static_cast<arith_concrete_base::elt &>(x));
        }
        static arith_generic::elt const & uncast(typename T::elt const & x) {
            return static_cast<arith_generic::elt const &>(static_cast<arith_concrete_base::elt const &>(x));
        }
    public:

    template<typename... Args> arith_wrapper(Args&&... args) : T(std::forward<Args>(args)...) {}
    virtual T * concrete() override { return dynamic_cast<T *>(this); }
    virtual T const * concrete() const override { return dynamic_cast<T const *>(this); }

    virtual void vec_add_and_reduce(elt * dst, elt const * b, size_t n) const override {
        concrete()->vec_add_and_reduce(cast(dst), cast(b), n);
    }
    virtual void add_and_reduce(elt & dst, elt const & b) const override {
        concrete()->add_and_reduce(cast(dst), cast(b));
    }
    virtual void sub_and_reduce(elt & dst, elt const & b) const override {
        concrete()->sub_and_reduce(cast(dst), cast(b));
    }

    virtual void vec_set(elt * x, elt const * a, size_t n) const override {
        concrete()->vec_set(cast(x), cast(a), n);
    }
    virtual void vec_neg(elt * x, elt const * a, size_t n) const override {
        concrete()->vec_neg(cast(x), cast(a), n);
    }
    virtual void vec_set_zero(elt * dst, size_t n) const override {
        concrete()->vec_set_zero(cast(dst), n);
    }
    virtual bool vec_is_zero(elt const * a, size_t n) const override {
        return concrete()->vec_is_zero(cast(a), n);
    }
    virtual void vec_set_random(elt * dst, size_t n, gmp_randstate_ptr rstate) const override {
        concrete()->vec_set_random(cast(dst), n, rstate);
    }
    virtual elt * vec_subvec(elt * a, size_t n) const override {
        return uncast(concrete()->vec_subvec(cast(a), n));
    }
    virtual elt const * vec_subvec(elt const * a, size_t n) const override {
        return uncast(concrete()->vec_subvec(cast(a), n));
    }
    virtual size_t vec_elt_stride(size_t n) const override {
        return concrete()->vec_elt_stride(n);
    }
    virtual void vec_dotprod(elt & w, elt const * u, elt const * v, size_t n) const override
    {
        concrete()->vec_dotprod(cast(w), cast(u), cast(v), n);
    }
    virtual void vec_addmul_and_reduce(elt * w, elt const * u, elt const & v, size_t n) const override
    {
        concrete()->vec_addmul_and_reduce(cast(w), cast(u), cast(v), n);
    }
    virtual int vec_cmp(elt const * a, elt const * b, size_t k) const override {
        return concrete()->vec_cmp(cast(a), cast(b), k);
    }
    virtual int cmp(elt const & a, elt const & b) const override {
        return concrete()->cmp(cast(a), cast(b));
    }
    virtual int cmp(elt const & a, unsigned long b) const override {
        return concrete()->cmp(cast(a), b);
    }

    virtual bool is_zero(elt const & x) const override { return concrete()->is_zero(cast(x)); }
    virtual void simd_set_ui_at(elt & a, size_t k, int v) const override {
        concrete()->simd_set_ui_at(cast(a), k, v);
    }
    virtual void simd_add_ui_at(elt & a, size_t k, int v) const override {
        concrete()->simd_set_ui_at(cast(a), k, v);
    }
    virtual int simd_hamming_weight(elt const & a) const override {
        return concrete()->simd_hamming_weight(cast(a));
    }
    virtual void vec_simd_set_ui_at(elt * a, size_t k, int v) const override {
        concrete()->vec_simd_set_ui_at(cast(a), k, v);
    }
    virtual void vec_simd_add_ui_at(elt * a, size_t k, int v) const override {
        concrete()->vec_simd_set_ui_at(cast(a), k, v);
    }
    virtual int vec_simd_hamming_weight(elt const * a, size_t n) const override {
        return concrete()->vec_simd_hamming_weight(cast(a), n);
    }
    virtual int vec_simd_find_first_set(elt & a, elt const * p, size_t n) const override {
        return concrete()->vec_simd_find_first_set(cast(a), cast(p), n);
    }
    virtual std::ostream& cxx_out(std::ostream& os, elt const & x) const override {
        return concrete()->cxx_out(os, cast(x));
    }
    virtual elt * alloc(size_t n = 1) const override {
        return uncast(concrete()->alloc(n));
    }
    virtual void free(elt * dst) const override {
        return concrete()->free(cast(dst));
    }
    virtual void set(elt & dst, elt const & src) const override {
        concrete()->set(cast(dst), cast(src));
    }
    virtual void set(elt & dst, cxx_mpz const & src) const override {
        concrete()->set(cast(dst), src);
    }
    virtual void neg(elt & dst, elt const & src) const override {
        concrete()->neg(cast(dst), cast(src));
    }
    virtual std::string impl_name() const override {
        return T::impl_name();
    }
    virtual size_t simd_groupsize() const override {
        return T::simd_groupsize;
    }
    virtual bool is_characteristic_two() const override {
        return T::is_characteristic_two;
    }
    virtual mpz_srcptr characteristic() const override {
        return concrete()->characteristic();
    }
    virtual void reduce(elt & a) const override {
        concrete()->reduce(cast(a));
    }
};

#ifdef ARITH_LAYER
extern "C" void * arith_layer(mpz_srcptr p, int simd_groupsize) {
#ifdef ARITH_MOD2
    return new arith_wrapper<arith_mod2::gf2<ARITH_SIMD_GROUPSIZE>>(p, simd_groupsize);
#elif defined(ARITH_MODP)
    return new arith_wrapper<arith_modp::gfp<ARITH_PRIME_WIDTH>>(p, simd_groupsize);
#else
#error "ARITH_LAYER requires ARITH_MOD2 or ARITH_MODP be defined"
#endif
}
#endif


arith_generic * arith_generic::instance(mpz_srcptr p, int simd_groupsize)
{
#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    std::string libname;
    if (mpz_cmp_ui(p, 2) == 0) {
        ASSERT_ALWAYS(simd_groupsize % 64 == 0);
        libname = fmt::format(FMT_STRING(SOLIB_PREFIX "arithmetic_b{}" SOLIB_SUFFIX), simd_groupsize);
    } else {
        ASSERT_ALWAYS(simd_groupsize == 1);
        ASSERT_ALWAYS(mpz_cmp_ui(p, 0) > 0);
        size_t nwords = mpz_size(p);
        libname = fmt::format(FMT_STRING(SOLIB_PREFIX "arithmetic_p{}" SOLIB_SUFFIX), nwords);
    }

    /* It might still be the case that we can have variable-width
     * interfaces after all. And these could conceivably be used as
     * fallbacks if the proper solib is not found.
     */
    void * handle = dlopen(libname.c_str(), RTLD_NOW);
    if (handle == NULL) {
        fprintf(stderr, "loading %s: %s\n", libname.c_str(), dlerror());
        abort();
    }
    typedef void * (*f_t)(mpz_srcptr, unsigned int);
    f_t f = (f_t) dlsym(handle, "arith_layer");
    if (f == NULL) {
        fprintf(stderr, "loading %s: %s\n", libname.c_str(), dlerror());
        abort();
    }
    return reinterpret_cast<arith_generic *>((*f)(p, simd_groupsize));
#else
    /* Then we need code that pulls in all implementation (from a
     * restrictive list), and selects which one should be returned.
     */
#define DO_b(g) do {							\
    if (mpz_cmp_ui(p, 2) == 0 && simd_groupsize == g)			\
        return new arith_wrapper<arith_mod2::gf2<g>>(p, simd_groupsize);\
    } while (0)

#define DO_p(s) do {							\
    if (mpz_size(p) == (s) && simd_groupsize == 1)			\
        return new arith_wrapper<arith_modp::gfp<s>>(p, simd_groupsize);\
    } while (0)

    COOKED_ARITHMETIC_BACKENDS;
    /*
    DO_b(64);
    DO_b(128);
    DO_b(192);
    DO_b(256);
    DO_p(1);
    DO_p(2);
    DO_p(3);
    DO_p(4);
    DO_p(5);
    DO_p(6);
    */

    abort();
#endif
}
