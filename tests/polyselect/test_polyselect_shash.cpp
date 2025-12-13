#include "cado.h" // IWYU pragma: keep

/* This is a c++ implementation of the exact same algorithm as in
 * polyselect_shash.c. (see comments there for the description of the
 * algorithm)
 *
 * Speed is like 5-10% slower.
 */

#define EMIT_ADDRESSABLE_shash_add
#define EXPOSE_DEPRECATED_polyselect_shash_find_collision

#include <ctime>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "polyselect_shash.h"
#include "macros.h"
#include "misc.h"

/*{{{ silly utility */
template<int n>
struct template_log2 {
    static_assert((n & (n-1)) == 0, "n must be a power of two");
    static_assert(n > 0, "n must be positive");
    static constexpr const int value = 1 + template_log2<n/2>::value;
};
template<>
struct template_log2<1> {
    static constexpr const int value = 0;
};
/*}}}*/

/*{{{ passed as a parameter to the template below */
struct polyselect_shash_config {
    using input_type = int64_t;
    using tie_breaker_type = uint32_t;
    static constexpr const unsigned int log2_open_hash_extra_size = 2;
    static constexpr const unsigned int open_hash_tail_overrun_protection = 16;
    static tie_breaker_type tie_breaker(input_type x) {
        return x + (x >> 32);
    }
};
/*}}}*/

template<typename config, int Nbuckets>
class bucket_hash {
    static_assert((Nbuckets & (Nbuckets-1)) == 0, "Nbuckets must be a power of two");
    static_assert(Nbuckets > 0, "Nbuckets must be positive");
    using T = typename config::input_type;
    using Ht = typename config::tie_breaker_type;
    static constexpr const unsigned int log2_nbuckets = template_log2<Nbuckets>::value;
    using U = typename std::make_unsigned<T>::type;
    T * data;
    T * base[Nbuckets + 1];
    T * current[Nbuckets + 1];
    size_t bucket_size;

    /* We forbid copies, since own the pointer (data), and we're too lazy
     * (and don't want) to add provisions for copying it
     */
    bucket_hash(bucket_hash const &) = delete;
    bucket_hash& operator= (bucket_hash const &) = delete;

    static size_t get_alloc_size(size_t expected_entries) {
        size_t alloc_size = expected_entries + 2 * std::sqrt(expected_entries);
        alloc_size *= 1.125;
        size_t const bucket_size = iceildiv(alloc_size, Nbuckets);
        alloc_size = Nbuckets * bucket_size;
        return alloc_size;
    }

    static size_t get_secondary_size(size_t bucket_size) {
        return next_power_of_2(bucket_size) << config::log2_open_hash_extra_size;
    }

    public:

    /* returns the memory cost of storing this table, assuming that it is
     * constructed with this expected_entries parameter */
    static size_t expected_memory_usage(size_t expected_entries) {
        size_t alloc_size = get_alloc_size(expected_entries);
        size_t bucket_size = alloc_size / Nbuckets;
        size_t A_size = get_secondary_size(bucket_size) + config::open_hash_tail_overrun_protection;;
        return alloc_size * sizeof(T) + A_size * sizeof(Ht);
    }

    bucket_hash(size_t expected_entries) {
        size_t const alloc_size = get_alloc_size(expected_entries);
        bucket_size = alloc_size / Nbuckets;
        data = new T[alloc_size];
        T * p = data;
        for(int i = 0 ; i < Nbuckets + 1 ; i++) {
            current[i] = base[i] = p;
            p += bucket_size;
        }
    }
    void reset() {
        for(int i = 0 ; i < Nbuckets ; i++)
            current[i] = base[i];
    }
    ~bucket_hash() {
        delete[] data;
    }
    void push(T const& x) {
        U const u = x;
        /* It's a matter of taste if we want to shift by log2_nbuckets
         * now or later. We'll compare bucket by bucket, so the equality
         * of all low-order bit is guaranteed.
         */
        U const q = u >> log2_nbuckets;
        U const r = u & (Nbuckets - 1);
        // if (UNLIKELY(current[r] >= base[r + 1])) throw std::runtime_error("argh");
        *current[r]++ = q;
    }
    bool has_collision() {
        /* Allocate an open hash table that is 4 times as large as what goes in
         * a typical bucket.
         */
        size_t const A_size = get_secondary_size(bucket_size);
        std::vector<Ht> A;
        for(int i = 0 ; i < Nbuckets ; i++) {
            A.assign(A_size + config::open_hash_tail_overrun_protection, 0);
            for(auto b = base[i] ; b != current[i] ; ++b) {
                T const x = *b;
                /* where do we insert x in A ?
                 *
                 * Note that we've done the shifting before storing x, so
                 * it's not needed here.
                 */
                size_t const where = (x) & (A_size -1);
                Ht const key = config::tie_breaker(x);
                Ht * Th = A.data() + where;
                for( ; *Th ; Th++)
                    if (*Th == key) {
                        return true;
                    }
                /* XXX Note that there's a possibility of overrunning
                 * here. We have an overrun protection guard above, but
                 * it's not necessarily enough */
                *Th = key;
            }
        }
        return false;
    }
};

int main()
{
    polyselect_shash_t H;

    cxx_gmp_randstate rstate;

    constexpr int64_t umax = INT64_C(2000000000000);
    constexpr int pushed_entries = 1000000;
    constexpr int expected_entries = pushed_entries;

    std::string collisions_c_code;
    {
        gmp_randseed_ui(rstate, 1);
        std::ostringstream os;

        polyselect_shash_init(H, expected_entries);

        int found = 0;
        int i;
        clock_t const st = clock();
        for(i = 0 ; i < 100 ; i++) {
            polyselect_shash_reset(H);
            for(size_t i = 0 ; i < pushed_entries ; i++)
                polyselect_shash_add(H, i64_random(rstate) % (2*umax) - umax);
            if (polyselect_shash_find_collision(H)) {
                found++;
                os << " " << i;
            }
        }
        std::string const s = os.str();
        printf("%d %d %.2f%s\n", i, found, (double) (clock() - st) / CLOCKS_PER_SEC, s.c_str());

        polyselect_shash_clear(H);

        collisions_c_code = s;
    }

    std::string collisions_cxx_code;
    {
        gmp_randseed_ui(rstate, 1);
        std::ostringstream os;
        bucket_hash<polyselect_shash_config, 256> H(expected_entries);
        int found = 0;
        int i;
        clock_t const st = clock();
        for(i = 0 ; i < 100 ; i++) {
            H.reset();
            for(size_t i = 0 ; i < pushed_entries ; i++) {
                // try {
                    H.push(i64_random(rstate) % (2*umax) - umax);
                    /*
                } catch(std::runtime_error const& e) {
                    fprintf(stderr, "overflow after %zu push's\n", i);
                }
                */
            }
            if (H.has_collision()) {
                found++;
                os << " " << i;
            }
            // printf("%d %d\n", i, found);
        }
        std::string const s = os.str();
        printf("%d %d %.2f%s\n", i, found, (double) (clock() - st) / CLOCKS_PER_SEC, s.c_str());

        collisions_cxx_code = s;
    }

    /* This code is functionally equivalent, but considerably slower */
    if (0) {
        gmp_randseed_ui(rstate, 1);
        std::ostringstream os;
        std::vector<int64_t> H;
        int found = 0;
        int i;
        clock_t const st = clock();
        for(i = 0 ; i < 100 ; i++) {
            H.clear();
            for(size_t i = 0 ; i < pushed_entries ; i++) {
                    H.push_back(i64_random(rstate) % (2*umax) - umax);
            }
            std::ranges::sort(H);
            for(auto it = ++H.begin() ; it != H.end() ; ++it) {
                if (it[0] == it[-1]) {
                    found++;
                    os << " " << i;
                    break;
                }
            }
            // printf("%d %d\n", i, found);
        }
        std::string const s = os.str();
        printf("%d %d %.2f%s\n", i, found, (double) (clock() - st) / CLOCKS_PER_SEC, s.c_str());
    }

    if (collisions_c_code != collisions_cxx_code) {
        fprintf(stderr, "The two implementations don't give matching results\n");
        exit(EXIT_FAILURE);
    }
    
    return 0;
}

