#ifndef CADO_MATMUL_HPP
#define CADO_MATMUL_HPP

#include <cstddef>
#include <cstdarg>
#include <cstdint>

/* This header is common to the different matrix product implementations
 */
#include "params.h"
#include "arith-generic.hpp"
#include "matrix_u32.hpp"
#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep

struct matmul_public {
    /* The fields here must be exposed by all implementations */
    std::array<unsigned int, 2> dim { 0, 0 };
                                /* dim[0] is nrows, dim[1] is ncols. Really. */
                                /* However, data needs not be stored in
                                 * that order */
    uint64_t ncoeffs { 0 };

    bool store_transposed { false }; /* For matrix building purposes, this
                                   indicates whether the implementation
                                   layer expects the flat matrix data in
                                   row major order (0) or transposed, in
                                   column major order (1). */
    std::array<int, 2> iteration { 0, 0 };
                                /* [0]: number of vector times matrix products
                                 * [1]: number of matrix times vector products
                                 */

    std::string locfile;
    bool no_save_cache { false };    /* if true, cache file is not saved */

    std::string cachefile_name;
    std::string local_cache_copy;

    std::string report_string;
};

struct matmul_interface : public matmul_public {
    protected:
    explicit matmul_interface(matmul_public && P)
        : matmul_public(std::move(P))
    {}
    public:

    int reload_cache();
    void save_cache();
    void save_to_local_copy();

    virtual void build_cache(matrix_u32 &&) = 0;
    virtual int reload_cache_private() = 0;
    virtual void save_cache_private() = 0;

    /* d == 1: this is a matrix-times-vector product.
     * d == 0: this is a vector-times-matrix product.
     *
     * In terms of number of rows and columns of vectors, this has
     * implications of course.
     *
     * A matrix times vector product takes a src vector of mm->dim[1]
     * coordinates, and outputs a dst vector of mm->dim[0] coordinates.
     *
     * This external interface is _not_ concerned with the value of the
     * mm->store_transposed flag. That one only relates to how the
     * implementation layer expects the data to be presented when the
     * cache is being built.
     */
    virtual void mul(void * dst, const void * src, int d) = 0;

    virtual void report(double) {}

    /* matmul_aux is used on some occasions for obtaining private
     * informations on the matmul structure. It's really a virtual method
     * hiding its name.
     *
     * _aux() is only the top-level. auxv is the one that does something.
     */
#define MATMUL_AUX_GET_READAHEAD        10
#define MATMUL_AUX_ZERO_STATS        11

    virtual void aux(int, ...) {} // NOLINT(cert-dcl50-cpp)

    virtual ~matmul_interface() = default;

    matmul_interface(matmul_interface const &) = delete;
    matmul_interface& operator=(matmul_interface const &) = delete;

    matmul_interface(matmul_interface &&) = default;
    matmul_interface& operator=(matmul_interface &&) = default;

    /* This function is the main code multiplexing point. It gets
     * called with a pointer to an abstract virtual base, and the actual
     * instantiation is a derived class that also inherits from a concrete
     * class. Internally, the low-level implementation will only see the
     * concrete class (with some wild type-punning).
     *
     * arith_generic (the vbase) has a virtual void * concrete() const member
     * function that is defined in the derived class, and returns the
     * concrete instance.
     *
     * Note that this _must_ be a shared_ptr, not a unique_ptr, for a few
     * reasons.
     *  - interfaces _can_ be shared across multiple threads in the
     *  context of interleaving. (but obviously not so for threads that
     *  deal with different matrix bits). Therefore, down the line, we're
     *  going to need a shared_ptr anyway.
     *  - we may or may not need a custom deleter. We need one if we're
     *  using shared libraries. In which case, it's slightly simpler to
     *  mask everything under a generous shared_ptr
     *
     * Pay attention to the caveat that the removal of the shared library
     * must appear _after_ the instance has been deleted and not, for
     * example, as a dtor action using a handle that would be stored in
     * matmul_public. In the latter case, stack unwinding _after_ the
     * dlclose would go to Babylon.
     */
    static std::shared_ptr<matmul_interface> create( // formerly matmul_init
            arith_generic * x,
            unsigned int nr,
            unsigned int nc,
            std::string const & locfile,
            std::string const & impl,
            cxx_param_list & pl,
            int optimized_direction);

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    struct remove_shared_lib_after_object {
        /* in the case of shared lib loading, we need to remove the
         * shared library as well.
         * Note that we may have a gazillion handles to the shared
         * library, in fact.
         */
        void * solib_handle;
        void operator()(matmul_interface * mm);
    };
#endif  /* BUILD_DYNAMICALLY_LINKABLE_BWC */
};


extern void matmul_decl_usage(cxx_param_list & pl);

extern void matmul_lookup_parameters(cxx_param_list & pl);

typedef matmul_interface * matmul_interface_ctor_t(matmul_public &&, arith_generic *, cxx_param_list &, int);


#endif	/* MATMUL_HPP_ */
