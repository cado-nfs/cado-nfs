#ifndef UTILS_BUCKET_HASH_HPP_
#define UTILS_BUCKET_HASH_HPP_

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "macros.h"
#include "misc.h"
#include "utils_cxx.hpp"

namespace cado
{

namespace bucket_hash_details
{
/*{{{ silly utility */
template <int n> struct template_log2 {
    static_assert((n & (n - 1)) == 0, "n must be a power of two");
    static_assert(n > 0, "n must be positive");
    static constexpr int const value = 1 + template_log2<n / 2>::value;
};
template <> struct template_log2<1> {
    static constexpr int const value = 0;
};
/*}}}*/

template <typename T, bool owns> struct data_field_type {
    using type = std::unique_ptr<T[]>;
};
template <typename T> struct data_field_type<T, false> {
    using type = T *;
};
template <typename T, bool owns>
using data_field_type_t = data_field_type<T, owns>::type;

struct bucket_full : std::runtime_error {
    bucket_full()
        : std::runtime_error("bucket full")
    {
    }
};

} /* namespace bucket_hash_details */

/*{{{ passed as a parameter to the template below */
struct polyselect_shash_config {
    static constexpr int const Nbuckets = 256;
    using input_type = int64_t;
    using aux_type = uint32_t;
    using tie_breaker_type = uint32_t;
    static constexpr unsigned int const log2_open_hash_extra_size = 2;
    static constexpr unsigned int const open_hash_tail_overrun_protection = 16;
    static tie_breaker_type tie_breaker(input_type x) { return x + (x >> 32); }
};
/*}}}*/

/* The first class citizen is bucket_hash, with owns=true. This class
 * owns and manages its storage.
 *
 * A bucket_hash_array is a generalization of bucket_hash. It's a bunch
 * of them pieced together, where a priori we expect different threads to
 * touch different (sub-) hash tables. The subtletly is that the storage is
 * managed at the bucket_hash_array level instead, hence the bucket_hash
 * elements are instantiated with owns=false. A bucket_hash_array also has
 * the capacity to adjust its internal hash tables in order to adjust to
 * the number of worker threads possibly changing over time.
 */
template<typename config> class bucket_hash_array;

template<typename config, bool owns=true>
class bucket_hash {
    friend class bucket_hash_array<config>;
    static constexpr const int Nbuckets = config::Nbuckets;
    static_assert((Nbuckets & (Nbuckets-1)) == 0, "Nbuckets must be a power of two");
    static_assert(Nbuckets > 0, "Nbuckets must be positive");
    using T = typename config::input_type;
    /* aux_type is used in the finer-grain search */
    using X = typename config::aux_type;
    using Ht = typename config::tie_breaker_type;

    /* this is used for actually finding collisions. There are two
     * options. Either we guarantee exact collisions, at the expense of
     * storing the full x, or we allow some looseness, which saves a
     * little bit of memory and costs us a little computation downstream.
     * The memory savings are marginal anyway.
     */
    struct HX_t {
        T i = 0;
        X k = 0;
    };

    static constexpr const unsigned int log2_nbuckets = bucket_hash_details::template_log2<Nbuckets>::value;
    using U = std::make_unsigned_t<T>;
    size_t bucket_size;
    bucket_hash_details::data_field_type_t<T, owns> data;
    bucket_hash_details::data_field_type_t<X, owns> pdata;
    T * base[Nbuckets + 1];
    T * current[Nbuckets + 1];
    T const * data_start() const
        requires(owns)
    {
        return data.get();
    }
    T const * data_start() const
        requires(!owns)
    {
        return data;
    }

    /* We forbid copies, since own the pointer (data), and we're too lazy
     * (and don't want) to add provisions for copying it
     */

    static size_t get_bucket_size(size_t expected_entries)
    {
        size_t alloc_size = expected_entries + 2 * std::sqrt(expected_entries);
        alloc_size += alloc_size / 4;
        return iceildiv(alloc_size, Nbuckets * 128) * 128;
    }

    static size_t get_secondary_size(size_t bucket_size)
    {
        return next_power_of_2(bucket_size)
               << config::log2_open_hash_extra_size;
    }

  public:
    /* returns the memory cost of storing this table, assuming that it is
     * constructed with this expected_entries parameter */
    static size_t expected_memory_usage(size_t expected_entries)
    {
        size_t bucket_size = get_bucket_size(expected_entries);
        size_t alloc_size = bucket_size * Nbuckets;
        size_t A_size = get_secondary_size(bucket_size) +
                        config::open_hash_tail_overrun_protection;
        ;
        return alloc_size * sizeof(T) + A_size * sizeof(Ht);
    }

    explicit bucket_hash(size_t expected_entries)
        requires(owns)
        : bucket_size(get_bucket_size(expected_entries))
        , data(std::make_unique<T[]>(get_bucket_size(expected_entries) * Nbuckets))
        , pdata(std::make_unique<X[]>(get_bucket_size(expected_entries) * Nbuckets))
    {
        T * p = data.get();
        for (int i = 0; i < Nbuckets + 1; i++) {
            current[i] = base[i] = p;
            p += bucket_size;
        }
    }

    /* this little "private key" idiom allows us to use our private ctor
     * even in emplace_back.
     */
  private:
    struct private_tag {
    };

  public:
    explicit bucket_hash(private_tag &&, size_t expected_entries, T * data,
                         X * pdata = nullptr)
        requires(!owns)
        : bucket_size(get_bucket_size(expected_entries))
        , data(data)
        , pdata(pdata)
    {
        T * p = data;
        for (int i = 0; i < Nbuckets + 1; i++) {
            current[i] = base[i] = p;
            p += bucket_size;
        }
    }

    void reset()
    {
        for (int i = 0; i < Nbuckets; i++)
            current[i] = base[i];
    }
    bool no_overflow() const
    {
        for (int i = 0; i < Nbuckets; i++)
            if (current[i] > base[i + 1])
                return false;
        return true;
    }
    void push(T const & x)
    {
        U const u = x;
        /* It's a matter of taste if we want to shift by log2_nbuckets
         * now or later. We'll compare bucket by bucket, so the equality
         * of all low-order bit is guaranteed.
         */
        U const q = u >> log2_nbuckets;
        U const r = u & (Nbuckets - 1);
#ifndef NDEBUG
        if (UNLIKELY(current[r] >= base[r + 1]))
            throw bucket_hash_details::bucket_full();
#endif
        *current[r]++ = q;
    }

    void push(std::pair<X, T> const & pi)
    {
        /* p is actually any label data of type aux_type. We do not
         * analyze it.
         */
        auto const & [p, x] = pi;
        U const u = x;
        U const q = u >> log2_nbuckets;
        U const r = u & (Nbuckets - 1);
#ifndef NDEBUG
        if (UNLIKELY(current[r] >= base[r + 1]))
            throw bucket_hash_details::bucket_full();
#endif
        pdata[current[r] - data_start()] = p;
        *current[r]++ = q;
    }

    template <typename RangeType>
    /* look for a collision in any of the buckets numbered from i0 to i1
     * within the union of the hash tables in B */
    static bool has_collision(RangeType const & B, unsigned int i0, unsigned int i1)
    {
        size_t sum_bucket_sizes = 0;
        for (auto const & b: B)
            sum_bucket_sizes += b.bucket_size;

        /* Allocate an open hash table that is 4 times as large as what goes in
         * a typical bucket.
         */
        size_t const A_size = get_secondary_size(sum_bucket_sizes);
        size_t const A_mask = A_size - 1;
        std::vector<Ht> A;
        for(auto i : std::views::iota(i0, i1)) {
            A.assign(A_size + config::open_hash_tail_overrun_protection, 0);
            for(auto const & Bj : B) {
                for(auto const & x : std::ranges::subrange(Bj.base[i], Bj.current[i])) {
                    /* where do we insert x in A ?
                     *
                     * Note that we've done the shifting before storing x, so
                     * it's not needed here.
                     */
                    Ht * Th = A.data() + (x & A_mask);
                    Ht const key = config::tie_breaker(x);
                    for (; *Th; Th++)
                        if (*Th == key)
                            return true;
                    /* XXX Note that there's a possibility of overrunning
                     * here. We have an overrun protection guard above, but
                     * it's not necessarily enough */
                    *Th = key;
                }
            }
        }
        return false;
    }
    template<typename RangeType, typename Iterator>
    static size_t find_collision(RangeType const & B, unsigned int i0, unsigned int i1, Iterator q)
    requires requires { *q++ = { X(), X(), T() }; }
    {
        size_t sum_bucket_sizes = 0;
        for (auto const & b: B)
            sum_bucket_sizes += b.bucket_size;

        size_t const A_size = get_secondary_size(sum_bucket_sizes);
        size_t const A_mask = A_size - 1;
        std::vector<HX_t> A;
        size_t ncollisions = 0;
        for (auto i: std::views::iota(i0, i1)) {
            A.assign(A_size + config::open_hash_tail_overrun_protection, {});
            for (auto const & Bj: B) {
                for (auto const & x:
                     std::ranges::subrange(Bj.base[i], Bj.current[i])) {
                    X k2 = Bj.pdata[&x - Bj.data_start()];
                    auto * Th = A.data() + (x & A_mask);
                    for (; Th->i; Th++) {
                        /* We do the comparison on the full (64-bit) x
                         * here, not on the 32-bit snapshot.
                         */
                        if (Th->i == x) {
                            /* record the collision between k1,i and k2,i. We
                             * have to recover the original i first.
                             */
                            auto k1 = Th->k;
                            *q++ = {k1, k2, (x << log2_nbuckets) | i};
                            ncollisions++;
                        }
                    }
                    *Th = {x, k2};
                }
            }
        }
        return ncollisions;
    }
    /* There can be sense in interleaving address computation (of Th,
     * which can include prefetching) and action. Here's an example of
     * how this would go. However all our attempts to make that efficient
     * have been unfruitful. All our measurements indicate that
     * performance is made worse, not better, and that the actual
     * interleaving level does not have much impact.
     */
    template <int interleaving>
    bool has_collision_prefetch()
    {
        size_t const A_size = get_secondary_size(bucket_size);
        size_t const A_mask = A_size - 1;
        std::vector<Ht> A;
        for(int i = 0 ; i < Nbuckets ; i++) {
            A.assign(A_size + config::open_hash_tail_overrun_protection, 0);
            auto read_address = [&](T x) {
                Ht * Th = A.data() + (x & A_mask);
                __builtin_prefetch(Th, 1, 3);
                return Th;
            };

            auto do_test = [&](Ht * Th, T x) {
                Ht const key = config::tie_breaker(x);
                for (; *Th; Th++)
                    if (*Th == key)
                        return true;
                *Th = key;
                return false;
            };

            auto const * b = base[i];
            for (; (current[i] - b) % interleaving; ++b) {
                if (do_test(read_address(*b), *b))
                    return true;
            }
            std::array<Ht *, interleaving> Th {};
            std::array<T, interleaving> x {};
            for (int i = 0; i < interleaving; i++)
                Th[i] = read_address(x[i] = b[i]);
            for (; b != current[i]; b += interleaving) {
                for (int i = 0; i < interleaving; i++) {
                    if (do_test(Th[i], x[i]))
                        return true;
                    Th[i] = read_address(x[i] = b[interleaving + i]);
                }
            }
        }
        return false;
    }

    bool has_collision()
    {
        auto r = std::ranges::subrange(this, this + 1);
        return has_collision(r, 0, Nbuckets);
    }

    /* this mimicks polyselect_shash2_find_collision_multi.
     *
     * the coding scheme in that function is a bit weird. In spirit, it
     * resembles what we now do with config::tie_breaker, with a minor
     * difference.
     *
     * */
    template<typename Iterator>
    size_t find_collisions(Iterator q) const
    // requires { *q++ = { std::declval<X>(), std::declval<X>(), std::declval<T>() }; }
    requires requires { *q++ = { X(), X(), T() }; }
    {
        auto r = std::ranges::subrange(this, this+1);
        return find_collisions(r, 0, Nbuckets, q);
    }

};

template<typename config>
class bucket_hash_array {
    public:
    static constexpr const int Nbuckets = config::Nbuckets;
    private:
    using T = typename config::input_type;
    using X = typename config::aux_type;
    using Ht = typename config::tie_breaker_type;
    using bhash_t = bucket_hash<config, false>;
    public:
    private:

    size_t total_expected_entries;
    size_t total_alloc_size;
    std::unique_ptr<T[]> storage;
    std::unique_ptr<X[]> pstorage;
    std::vector<bhash_t> B;

    static size_t get_alloc_size(size_t total_expected_entries,
                                 unsigned int multi)
    {
        size_t s = 0;
        for (unsigned int i = 1; i <= multi; i++) {
            /* expected number of entries per table when we use i
             * different tables.
             */
            size_t const Ei = iceildiv(total_expected_entries, i);
            s = std::max(s, Nbuckets * i * bhash_t::get_bucket_size(Ei));
        }
        return s;
    }

  public:
    bool has_collision(unsigned int i0, unsigned int i1)
    {
        return bhash_t::has_collision(B, i0, i1);
    }

    std::string stats() const {
        size_t min = SIZE_MAX, max = 0;
        double s1 = 0, s2 = 0;
        size_t ts = 0;
        for(auto const & Bj : B) {
            for(int i = 0 ; i < Nbuckets ; i++) {
                size_t s = Bj.current[i] - Bj.base[i];
                ts += s;
                s2 += (double) s * (double) s;
                min = std::min(min, s);
                max = std::max(max, s);
            }
        }
        s1 = double_ratio(ts, B.size() * Nbuckets);
        s2 = double_ratio(s2, B.size() * Nbuckets);
        s2 = sqrt(s2 - s1 * s1);
        return fmt::format("[{}*{} samples, {} expected] min={} max={} total={} mean={:.2f} sdev={:.2f}",
                B.size(), Nbuckets, total_expected_entries, min, max, ts, s1, s2);

    }

    template<typename Iterator>
    size_t find_collision(unsigned int i0, unsigned int i1, Iterator q)
    requires requires { *q++ = { X(), X(), T() }; }
    {
        return bhash_t::find_collision(B, i0, i1, q);
    }

    bucket_hash_array(size_t total_expected_entries, unsigned int multi)
        : total_expected_entries(total_expected_entries)
        , total_alloc_size(get_alloc_size(total_expected_entries, multi))
        , storage(std::make_unique<T[]>(total_alloc_size))
        , pstorage(std::make_unique<X[]>(total_alloc_size))
    {
        B.reserve(multi);
        reconfigure(multi);
    }

    void reconfigure(unsigned int multi)
    {
        size_t const s = total_alloc_size / multi;
        B.clear();
        ASSERT_ALWAYS(multi <= B.capacity());
        size_t const Ei = iceildiv(total_expected_entries, multi);
        ASSERT_ALWAYS(Nbuckets * bhash_t::get_bucket_size(Ei) <= s);
        for(unsigned int i = 0 ; i < multi ; i++)
            B.emplace_back(
                    typename bhash_t::private_tag(),
                    total_expected_entries / multi,
                    storage.get() + i * s,
                    pstorage.get() + i * s);
    }

    bhash_t & operator[](unsigned int i) { return B[i]; }
    bhash_t const & operator[](unsigned int i) const { return B[i]; }
};

} /* namespace cado */

#endif	/* UTILS_BUCKET_HASH_HPP_ */
