#ifndef FB_SLICE_WEIGHT_HPP_
#define FB_SLICE_WEIGHT_HPP_

#include "fb.hpp"

template<typename FB_ENTRY_TYPE>
class fb_slice_weight_estimator {
    typedef typename std::vector<FB_ENTRY_TYPE>::const_iterator it_t;
    typedef fb_slice<FB_ENTRY_TYPE> sl_t;

    /* Compute an upper bound on this slice's weight by assuming all *
     * primes in the slice are equal to the first (and thus smallest)
     * one. This relies on the vector being sorted. */
    inline double max(sl_t const &, it_t a, it_t b) const {
        return a->weight() * (b - a);
    }
    inline double max(sl_t const & sl) const { return max(sl, sl._begin, sl._end); }

    /* Estimate weight by the average of the weights of the two
     * endpoints, i.e., by trapezoidal rule. */
    inline double avg(sl_t const &, it_t a, it_t b) const {
        return (a->weight() + (--b)->weight()) / 2. * (b - a);
    }
    inline double avg(sl_t const & sl) const { return avg(sl, sl._begin, sl._end); }

    /* Estimate weight by Simpson's rule on the weight of the two
     * endpoints and of the midpoint */
    inline double simpson(sl_t const &, it_t a, it_t b) const {
      it_t mid = a + (b-a) / 2;
      return (a->weight() + 4.*mid->weight() + (--b)->weight()) / 6. * (b-a);
    }
    inline double simpson(sl_t const & sl) const { return simpson(sl, sl._begin, sl._end); }

    /* Estimate weight by using Merten's rule on the primes at the two
     * endpoints */
    inline double mertens(sl_t const &, it_t a, it_t b) const {
        return log(log((--b)->get_q())) - log(log(a->get_q()));
    }
    inline double mertens(sl_t const & sl) const { return mertens(sl, sl._begin, sl._end); }

    /* Compute weight exactly with a sum over all entries */
    inline double exact(sl_t const &, it_t a, it_t b) const {
        double s = 0;
        for(it_t x = a ; x != b ; ++x) s += x->weight();
        return s;
    }
    inline double exact(sl_t const & sl) const { return exact(sl, sl._begin, sl._end); }

    /* this is in the C file */
    double compare(sl_t const &, it_t a, it_t b) const;
    double compare(sl_t const & sl) const;

    public:
    double operator()(sl_t const &, it_t a, it_t b) const;

    inline double operator()(sl_t const & sl) const {
        return operator()(sl, sl._begin, sl._end);
    }
};

/* For general vectors, we compute the weight the hard way, via a sum over
   all entries. General vectors are quite small so this should not take long */
template <>
class
fb_slice_weight_estimator<fb_entry_general> {
    typedef typename fb_slice<fb_entry_general>::fb_entry_vector::const_iterator it_t;
    typedef fb_slice<fb_entry_general> sl_t;
    public:
    double operator()(sl_t const &, it_t a, it_t b) const {
        double s = 0;
        for(it_t x = a ; x != b ; ++x) s += x->weight();
        return s;
    }
    inline double operator()(sl_t const & sl) const {
        return operator()(sl, sl._begin, sl._end);
    }
};

#endif	/* FB_SLICE_WEIGHT_HPP_ */
