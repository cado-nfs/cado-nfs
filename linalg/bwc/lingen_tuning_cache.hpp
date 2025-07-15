#ifndef CADO_LINGEN_TUNING_CACHE_HPP
#define CADO_LINGEN_TUNING_CACHE_HPP

#include <cstddef>

#include <iosfwd>
#include <map>
#include <array>
#include <list>
#include <utility>

#include "lingen_mul_substeps_base.hpp"
#include "lingen_substep_schedule.hpp"
#include "lingen_round_operand_size.hpp"

struct lingen_tuning_cache {
    struct basecase_key {
        size_t np;
        unsigned int m, n;
        size_t length;
        int nthreads;
    };
    struct mul_or_mp_key {
        op_mul_or_mp_base::op_type_t op_type;
        lingen_substep_schedule::fft_type_t fft_type;
        size_t np, na, nb;
    };
    typedef double basecase_value;
    typedef std::array<std::list<std::pair<unsigned int, double>>, 4> mul_or_mp_value;
    // typedef mul_or_mp_value mul_value;
    // typedef mul_or_mp_value mp_value;

    struct coarse_compare {
        bool operator()(size_t const &a, size_t const& b) const {
            size_t ca = lingen_round_operand_size(a);
            size_t cb = lingen_round_operand_size(b);
            return ca < cb;
        }
        bool operator()(basecase_key const &a, basecase_key const& b) const {
            size_t cla = lingen_round_operand_size(a.length);
            size_t clb = lingen_round_operand_size(b.length);
            return cla < clb;
        }
        bool operator()(mul_or_mp_key const &a, mul_or_mp_key const& b) const {
            if (a.op_type < b.op_type) return true;
            if (a.op_type > b.op_type) return false;
            if (a.fft_type < b.fft_type) return true;
            if (a.fft_type > b.fft_type) return false;
            size_t cana = lingen_round_operand_size(a.na);
            size_t cbna = lingen_round_operand_size(b.na);
            if (cana < cbna) return true;
            if (cana > cbna) return false;
            size_t canb = lingen_round_operand_size(a.nb);
            size_t cbnb = lingen_round_operand_size(b.nb);
            if (canb < cbnb) return true;
            if (canb > cbnb) return false;
            return false;
        }
    };

    std::map<basecase_key, basecase_value, coarse_compare> basecase_cache;
    std::map<mul_or_mp_key, mul_or_mp_value, coarse_compare> mul_or_mp_cache;

    void load(const char * timing_cache_filename);
    void save(const char * timing_cache_filename);

    bool has(basecase_key const & K) const { return basecase_cache.find(K) != basecase_cache.end(); };
    bool has(mul_or_mp_key const & K) const {
        return mul_or_mp_cache.find(K) != mul_or_mp_cache.end();
    };
    basecase_value & operator[](basecase_key const & K) { return basecase_cache[K]; }
    mul_or_mp_value & operator[](mul_or_mp_key const & K) { return mul_or_mp_cache[K]; }
};

std::istream& operator>>(std::istream& is, lingen_tuning_cache::basecase_key &);
std::istream& operator>>(std::istream& is, lingen_tuning_cache::mul_or_mp_key &);
std::ostream& operator<<(std::ostream& os, lingen_tuning_cache::basecase_key const &);
std::ostream& operator<<(std::ostream& os, lingen_tuning_cache::mul_or_mp_key const &);

#endif	/* LINGEN_TUNING_CACHE_HPP_ */
