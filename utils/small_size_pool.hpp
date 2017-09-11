#ifndef SMALL_SIZE_POOL_HPP_
#define SMALL_SIZE_POOL_HPP_

#include "macros.h"
#include <string.h>
#include <vector>
#include <map>
#include <limits>
#include <memory>
#ifdef DEBUG_SMALL_SIZE_POOL
#include <string>
#include <ostream>
#endif
#ifdef HAVE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace std { using boost::shared_ptr; using boost::make_shared; }
#endif

/* in a way, the code here can be regarded as a simple-minded
 * implementation of malloc(), which comes with some extra benefits and
 * some extra requirements.
 *
 * For each "pointer", we always want to know both the number of pointed
 * items (that is, the size of the area), as well as the internal pointer
 * value, which is an integer type (morally a size_t, althoughthe second
 * argument to the template class can be used to make it tighter). This
 * internal pointer value compares equal to 0 whenever there's no data.
 *
 * The data allocated on the heap is only global per range size (the
 * third optional parameter can be used to make this coarser. By default
 * it's set to 1. Change at your own risk).
 *
 * There is no extra allocation cost per pointed range.
 *
 * From a structure such as:
 *      small_size_pool<int, size_t, 4> G;
 * we will get int* ranges which hold strings of integer values whose
 * lengths are runtime-decided (good practice recommends that the runtime
 * application is able to guarantee that these lengths take relatively
 * few different values).
 *
 * The 4 indicates that sizes which round up similarly to the same next
 * multiple of 4 are considered at the same time, and allocated
 * similarly.
 *
 * The size_t argument indicates that we expect that a size_t will be
 * necessary to enumerate all the strings of a given length (or
 * congruence class of length). For most uses though, size_t is way
 * overkill, you're probably better off with uint32_t.
 *
 * Typical calls are:
 *      size_t value = G.get(17);
 *      int * v = G(17, value);
 *      G.free(17, value);
 *
 *      size_t othervalue = G.get(15);
 *      int * w = G(15, othervalue);
 *      G.free(15, othervalue);
 *
 *      size_t width = 22;
 *      size_t bigvalue = G.get(22);
 *      = G.realloc(width, bigvalue, 11);
 *      // now width == 11, and bigvalue has changed
 *      G.free(width, bigvalue);
 */
/* Make sure you do not take pointers to this structure, as there may be
 * reallocs */
template<typename T, typename S = typename std::vector<T>::size_type> struct single_size_pool {
    typedef S size_type;
    S width;
    typedef std::vector<T> container_type;
    static size_type min_width() { return (sizeof(S) + sizeof(T) - 1) / sizeof(T); }
    private:
    container_type data;
    static constexpr size_type size_max = std::numeric_limits<S>::max();
    size_type holes;
    single_size_pool(single_size_pool const&);
    single_size_pool& operator=(single_size_pool const&);
    public:
    single_size_pool(S width) : width(width), data(width), holes(0) {
        memset(&data[0], 0, width * sizeof(T));
        ASSERT_ALWAYS(width >= min_width());
    }
    T * operator[](size_type i) { return &data[i*width]; }
    T const * operator[](size_type i) const { return &data[i*width]; }
    size_type get_fresh() {
        size_type & next_free = *(size_type*)(&data[0]);
        if (next_free) {
            size_type ret = next_free;
            next_free = *(size_type*)(&data[ret*width]);
            holes--;
            return ret;
        } else {
            data.insert(data.end(), width, T());
            return data.size()/width-1;
        }
    }
    size_type release(size_type p) {
        ASSERT_ALWAYS(p > 0 && p*width < data.size());
        size_type & next_free = *(size_type*)(&data[0]);
        *(size_type*)(&data[p*width]) = next_free;
        holes++;
        next_free = p;
        return 0;
    }
    size_type size() const { return data.size() / width; }
    size_t allocated_bytes() const { return data.capacity()*sizeof(T) + sizeof(*this); }

#ifdef DEBUG_SMALL_SIZE_POOL
    std::ostream& print(std::ostream& o) const {
        o << "width " << width
            << ", allocated " << size()
            << ", capacity " << data.capacity() / width
            << " (" << (allocated_bytes() >> 20) << " MB)"
            << ", holes " << holes;
        /*
        if (holes) {
            o << ":";
            for(size_type next = *(size_type*)(&data[0]) ; next ; ) {
                o << " " << next;
                next= *(size_type*)(&data[next*width]);
            }
        }
        */
        o << "\n";
        return o;
    }
#endif
};
#ifdef DEBUG_SMALL_SIZE_POOL
template<typename T, typename S>
std::ostream& operator<<(std::ostream& o, single_size_pool<T,S> const& p) {
    return p.print(o);
}
template<typename T> struct small_size_pool_printer;
#endif


template<typename T, typename S = typename std::vector<T>::size_type, int coarse=1>
struct small_size_pool {
    typedef T value_type;
    typedef S size_type;
protected:
    typedef small_size_pool<T, S, coarse> self;
    static_assert(coarse > 0, "\"coarse\" template parameter must be >0");
    typedef single_size_pool<T, S> spool;
    std::map<S, std::shared_ptr<spool>> pools;
    static S get_coarse(S size) {
        size = std::max(size, spool::min_width());
        return size-1 + coarse - ((size-1) % coarse);
    }
public:
    void clear() { pools.clear(); }
    S alloc(S const& size) {
        if (!size) return 0;
        int csize = get_coarse(size);
        if (!pools[csize])
            pools[csize] = std::make_shared<spool>(csize);
        return pools[csize]->get_fresh();
    }
    T * operator()(S size, S value) {
        if (!size || !value) return NULL;
        int csize = get_coarse(size);
        ASSERT_ALWAYS(pools[csize]);
        return (*pools[csize])[value];
    }
    T * operator()(std::pair<S, S> const& p) {
        return (*this)(p.first, p.second);
    }
    const T * operator()(S size, S value) const {
        if (!size || !value) return NULL;
        int csize = get_coarse(size);
        return (*pools.at(csize))[value];
    }
    const T * operator()(std::pair<S, S> const& p) const {
        return (*this)(p.first, p.second);
    }
    void free(S & size, S & value) {
        ASSERT_ALWAYS(size || !value);
        if (!size || !value) return;
        int csize = get_coarse(size);
        ASSERT_ALWAYS(pools[csize]);
        value = pools[csize]->release(value);
        size = 0;
    }
    void free(std::pair<S, S> & p) {
        free(p.first, p.second);
    }
    void realloc(S & size, S & value, S newsize) {
        if (!value || !size) {
            value = alloc(size = newsize);
            return;
        }
        if (!newsize) {
            free(size, value);
            return;
        }
        ASSERT_ALWAYS(size);
        ASSERT_ALWAYS(newsize);
        int csize = get_coarse(size);
        int cnewsize = get_coarse(newsize);
        if (csize == cnewsize) {
            size = newsize;
            return;
        }
        ASSERT_ALWAYS(pools[csize]);
        ASSERT_ALWAYS(value);
        if (!pools[cnewsize])
            pools[cnewsize] = std::make_shared<spool>(cnewsize);
        const T * oldptr = (*pools[csize])[value];
        S newvalue = pools[cnewsize]->get_fresh();
        T * newptr = (*pools[cnewsize])[newvalue];
        std::copy(oldptr, oldptr + std::min(size, newsize), newptr);
        pools[csize]->release(value);
        value = newvalue;
        size = newsize;
    }
    void realloc(std::pair<S, S> & p, S newsize) {
        realloc(p.first, p.second, newsize);
    }
    size_t allocated_bytes() const {
        size_t res = sizeof(*this);
        for(auto const & x : pools)
            res += x.second->allocated_bytes();
        return res;
    }
#ifdef DEBUG_SMALL_SIZE_POOL
    friend class small_size_pool_printer<self>;
    small_size_pool_printer<self> printer(std::string const& s) const { return small_size_pool_printer<self>(*this, s); }
    std::ostream& print(std::ostream& o) const { return o << printer(""); }
#endif
};
#ifdef DEBUG_SMALL_SIZE_POOL
template<typename T>
struct small_size_pool_printer {
    T const & P;
    std::string prefix;
    small_size_pool_printer(T const& P, std::string const& s) : P(P), prefix(s) {}
    std::ostream& print(std::ostream& o) const {
        typename T::size_type total = 0;
        for(auto const & x : P.pools) {
            o << prefix << *x.second;
            total += x.second->allocated_bytes();
        }
        o << prefix << "Total allocated size for all widths: " << (total >> 20) << " MB\n";
        return o;
    }
};
template<typename T, typename S, int c>
std::ostream& operator<<(std::ostream& o, small_size_pool<T,S,c> const& p) {
    return p.print(o);
}
template<typename P>
std::ostream& operator<<(std::ostream& o, small_size_pool_printer<P> const& p) {
    return p.print(o);
}
#endif

#endif	/* SMALL_SIZE_POOL_HPP_ */
