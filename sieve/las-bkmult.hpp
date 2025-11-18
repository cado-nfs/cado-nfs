#ifndef CADO_LAS_BKMULT_HPP
#define CADO_LAS_BKMULT_HPP

#include <map>                     // for map<>::key_type, operator!=, _Rb_t...
#include <string>                  // for string, allocator
#include <utility>                 // for pair
#include "clonable-exception.hpp"  // for clonable_exception

class bkmult_specifier {
    double base = 1.0;
    using dict_t = std::map<std::pair<int, char>, double>;
    dict_t dict;
    public:
    using key_type = dict_t::key_type;
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
    ~buckets_are_full() override = default;
    buckets_are_full(buckets_are_full const &) noexcept = default;
    buckets_are_full(bkmult_specifier::key_type const&, int b, int r, int t);
    const char * what() const noexcept override { return message.c_str(); }
    bool operator<(buckets_are_full const& o) const {
        return (double) reached_size / theoretical_max_size < (double) o.reached_size / o.theoretical_max_size;
    }
    clonable_exception * clone() const override { return new buckets_are_full(*this); }
};


#endif	/* CADO_LAS_BKMULT_HPP */
