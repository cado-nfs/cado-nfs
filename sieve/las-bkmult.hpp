#ifndef LAS_BKMULT_HPP_
#define LAS_BKMULT_HPP_

#include <map>                     // for map<>::key_type, operator!=, _Rb_t...
#include <string>                  // for string, allocator
#include <utility>                 // for pair
#include "clonable-exception.hpp"  // for clonable_exception

class bkmult_specifier {
    double base = 1.0;
    typedef std::map<std::pair<int, char>, double> dict_t;
    dict_t dict;
    public:
    typedef dict_t::key_type key_type;
    static std::string printkey(dict_t::key_type const& key) {
        char c[3] = { (char) ('0' + key.first), key.second, '\0' };
        return std::string(c);
    }
    template<typename T> static dict_t::key_type getkey() {
        return dict_t::key_type(T::level(), T::rtti[0]);
    }
    template<typename T> double get() const { return get(getkey<T>()); }
    double const & get(dict_t::key_type const& key) const {
        auto xx = dict.find(key);
        if (xx != dict.end()) return xx->second;
        return base;
    }
    double grow(dict_t::key_type const& key, double d) {
        double v = get(key) * d;
        return dict[key] = v;
    }
    template<typename T> double get(T const &) const { return get<T>(); }
    template<typename T> double operator()(T const &) const { return get<T>(); }
    template<typename T> double operator()() const { return get<T>(); }
    bkmult_specifier() = default;
    bkmult_specifier(double x) : base(x) {}
    bkmult_specifier(const char * specifier);
    std::string print_all() const;
};

struct buckets_are_full : public clonable_exception {
    bkmult_specifier::key_type key;
    int bucket_number;
    int reached_size;
    int theoretical_max_size;
    std::string message;
    ~buckets_are_full();
    buckets_are_full(buckets_are_full const &);
    buckets_are_full(bkmult_specifier::key_type const&, int b, int r, int t);
    virtual const char * what() const noexcept { return message.c_str(); }
    bool operator<(buckets_are_full const& o) const {
        return (double) reached_size / theoretical_max_size < (double) o.reached_size / o.theoretical_max_size;
    }
    virtual clonable_exception * clone() const { return new buckets_are_full(*this); }
};


#endif	/* LAS_BKMULT_HPP_ */
