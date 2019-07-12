#ifndef PLINGEN_TUNING_CACHE_HPP_
#define PLINGEN_TUNING_CACHE_HPP_

#include <tuple>
#include <map>

#include "plingen-tuning.hpp"  // lingen_round_operand_size

struct plingen_tuning_cache {
    struct basecase_key : public std::tuple<size_t, unsigned int, unsigned int, size_t, int> {
        typedef std::tuple<size_t, unsigned int, unsigned int, size_t, int> super;
        basecase_key() = default;
        template<typename... Args>
            basecase_key(Args&&... args) : super(args...) {}
    };
    struct mul_key : public std::tuple<size_t, size_t, size_t> {
        typedef std::tuple<size_t, size_t, size_t> super;
        template<typename... Args>
            mul_key(Args&&... args) : super(args...) {}
    };
    struct mp_key : public std::tuple<size_t, size_t, size_t> {
        typedef std::tuple<size_t, size_t, size_t> super;
        template<typename... Args>
            mp_key(Args&&... args) : super(args...) {}
    };
    typedef double basecase_value;
    typedef std::array<double, 4> mul_value;
    typedef std::array<double, 4> mp_value;

    struct coarse_compare {
        bool operator()(size_t const &a, size_t const& b) const {
            size_t ca = lingen_round_operand_size(a);
            size_t cb = lingen_round_operand_size(b);
            return ca < cb;
        }
        bool operator()(basecase_key const &a, basecase_key const& b) const {
            basecase_key ca = a;
            basecase_key cb = b;
            std::get<3>(ca) = lingen_round_operand_size(std::get<3>(ca));
            std::get<3>(cb) = lingen_round_operand_size(std::get<3>(cb));
            return ca < cb;
        }
        bool operator()(mul_key const &a, mul_key const& b) const {
            mul_key ca = a;
            mul_key cb = b;
            std::get<1>(ca) = lingen_round_operand_size(std::get<1>(ca));
            std::get<2>(ca) = lingen_round_operand_size(std::get<2>(ca));
            std::get<1>(cb) = lingen_round_operand_size(std::get<1>(cb));
            std::get<2>(cb) = lingen_round_operand_size(std::get<2>(cb));
            return ca < cb;
        }
        bool operator()(mp_key const &a, mp_key const& b) const {
            mp_key ca = a;
            mp_key cb = b;
            std::get<1>(ca) = lingen_round_operand_size(std::get<1>(ca));
            std::get<2>(ca) = lingen_round_operand_size(std::get<2>(ca));
            std::get<1>(cb) = lingen_round_operand_size(std::get<1>(cb));
            std::get<2>(cb) = lingen_round_operand_size(std::get<2>(cb));
            return ca < cb;
        }
    };

    std::map<basecase_key, double, coarse_compare> basecase_cache;
    std::map<mul_key, std::array<double, 4>, coarse_compare> mul_cache;
    std::map<mp_key, std::array<double, 4>, coarse_compare> mp_cache;

    void load(const char * timing_cache_filename);
    void save(const char * timing_cache_filename);

    bool has(basecase_key const & K) const { return basecase_cache.find(K) != basecase_cache.end(); };
    bool has(mul_key const & K) const { return mul_cache.find(K) != mul_cache.end(); };
    bool has(mp_key const & K) const { return mp_cache.find(K) != mp_cache.end(); };
    basecase_value & operator[](basecase_key const & K) { return basecase_cache[K]; }
    mul_value & operator[](mul_key const & K) { return mul_cache[K]; }
    mp_value & operator[](mp_key const & K) { return mp_cache[K]; }
};

std::istream& operator>>(std::istream& is, plingen_tuning_cache::basecase_key &);
std::istream& operator>>(std::istream& is, plingen_tuning_cache::mul_key &);
std::istream& operator>>(std::istream& is, plingen_tuning_cache::basecase_key &);
std::ostream& operator<<(std::ostream& is, plingen_tuning_cache::basecase_key const &);
std::ostream& operator<<(std::ostream& is, plingen_tuning_cache::mul_key const &);
std::ostream& operator<<(std::ostream& is, plingen_tuning_cache::basecase_key const &);

#endif	/* PLINGEN_TUNING_CACHE_HPP_ */
