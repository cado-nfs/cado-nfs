#ifndef LAS_WHERE_AM_I_PROXY_HPP_
#define LAS_WHERE_AM_I_PROXY_HPP_

#include <memory>
#include <stdlib.h>
struct cxx_param_list;
class nfs_work;

/* Yes, where_am_I is an _incomplete_ type, on purpose. We only
 * manipulate opaque pointers. The compilation units that actually _do_
 * something with that are:
 *  - separate functions that do allocation/free ; there's a trace and a
 *    notrace version
 *  - template instantiations that also call more specific code in the
 *    las-debug.cpp module.
 *
 */

struct where_am_I {
    where_am_I();
    where_am_I(where_am_I const &);
    where_am_I & operator=(where_am_I const &);
    ~where_am_I();
    static void decl_usage(cxx_param_list &);
    static void interpret_parameters(cxx_param_list &);
    static void begin_special_q(nfs_work const &);
    private:
    struct impl;  // forward declaration of the implementation class
    // unique-ownership opaque pointer to the forward-declared
    // implementation class:
    typedef std::unique_ptr<impl> pimpl_t;
    pimpl_t pimpl;
    public:
    impl * operator->() { return pimpl.get(); };
    impl const * operator->() const { return pimpl.get(); };
};

/* Please, don't use this in tight loops. This just to avoid code bloat.
 * In cases where the GOT linking saves us a megabyte of duplicated code,
 * why not go for it ? */
int extern_trace_on_spot_ab(int64_t a, uint64_t b);

#endif	/* LAS_WHERE_AM_I_PROXY_HPP_ */
