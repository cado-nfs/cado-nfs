#ifndef CXX_MISC_HPP_
#define CXX_MISC_HPP_

/* miscellaneous stuff for C++ */

/* This is handy for range-based for loops that also need a counter.
 */
template<typename T>
struct increment_counter_on_dtor {
    T & a;
    increment_counter_on_dtor(T & a) : a(a) {}
    ~increment_counter_on_dtor() { ++a; }
};

/* usage example:
 *
    size_t counter = 0;
    for(auto x : foo) {
        increment_counter_on_dtor<size_t> _dummy(counter);
        if (x % 2 == 0) continue;
        std::cout << "[" << counter << "]: " << x << "\n";
    }

 * of course, &x-begin(foo) is also acceptable, but that works only for
 * random-access iterators.
 */

#endif	/* CXX_MISC_HPP_ */
