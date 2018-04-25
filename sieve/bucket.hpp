#ifndef BUCKET_HPP_
#define BUCKET_HPP_

/*
 * Bucket sieving: radix-sort sieve updates as they are created.
 */

#include <stdint.h>
#include "cado-endian.h"
#define xxxSAFE_BUCKETS
#ifdef SAFE_BUCKETS
#include <exception>
#include <stdio.h>
#include <limits>
#include <array>
#include <string>
#include "portability.h"
#endif
#include "misc.h"
#include "fb-types.h"
#include "fb.hpp"
#include "las-debug.hpp"
#include "threadpool.hpp"

#include "electric_alloc.h"

/*
 * This bucket module provides a way to store elements (that are called
 * updates), while partially sorting them, according to some criterion (to
 * be defined externally): the incoming data is stored into several
 * buckets. The user says for each data to which bucket it belongs. This
 * module is supposed to perform this storage in a cache-friendly way and
 * so on...
 */

/*
 * Main available commands are 
 *   push_bucket_update(i, x)  :   add x to bucket number i
 *   get_next_bucket_update(i) :   iterator to read contents of bucket number i
 *
 * See the MAIN FUNCTIONS section below for the complete interface, with
 * prototypes of exported functions.
 */

/********  Data structure for the contents of buckets **************/

/* In principle, the typedef for the bucket_update can be changed without
 * affecting the rest of the code. 
 */

/*
 * For the moment, we store the bucket updates and a 16-bit field
 * that can contain, for instance, the low bits of p.
 */

template <typename TARGET_TYPE, typename SOURCE_TYPE>
TARGET_TYPE limit_cast(const SOURCE_TYPE &b)
{
  ASSERT_EXPENSIVE(b >= std::numeric_limits<TARGET_TYPE>::min());
  ASSERT_EXPENSIVE(b <= std::numeric_limits<TARGET_TYPE>::max());
  return static_cast<TARGET_TYPE>(b);
}

/* 16-bits */
class shorthint_t {
public:
  slice_offset_t hint;
  shorthint_t() {}
  shorthint_t(const slice_offset_t slice_offset)
    : hint(slice_offset) {}
  shorthint_t(const fbprime_t p MAYBE_UNUSED,
              const slice_offset_t slice_offset,
              const slice_index_t slice_index MAYBE_UNUSED)
    : hint(slice_offset) {}
  static constexpr char const * rtti = "shorthint_t";
};

/* sizeof(slice_index_t), that is 4 bytes, + 2 hint bytes == 6 bytes */
/* When purging a bucket, we don't store pointer arrays to indicate where in
   the purged data a new slice begins, as each slice will have only very few
   updates surviving. Instead, we re-write each update to store both slice
   index and offset. */
class longhint_t {
public:
  slice_index_t index;
  slice_offset_t hint;
  longhint_t(){}
  longhint_t(const slice_offset_t slice_offset,
             const slice_index_t slice_index)
    : index(slice_index), hint(slice_offset) {}
  longhint_t(const fbprime_t p MAYBE_UNUSED,
             const slice_offset_t slice_offset,
             const slice_index_t slice_index)
    : index(slice_index), hint(slice_offset) {}
  static constexpr char const * rtti = "longhint_t";
};

/* An update with the complete prime, generated by line re-sieving */
class primehint_t {
public:
  fbprime_t p;
  primehint_t(){}
  primehint_t(const fbprime_t p)
    : p(p) {}
  primehint_t(const fbprime_t p,
              const slice_offset_t slice_offset MAYBE_UNUSED,
              const slice_index_t slice_index MAYBE_UNUSED)
    : p(p) {}
};

/* A bucket update type has two template parameters: the level of the bucket
   sieving where the update was created, and the type of factor base prime
   hint it stores, used to speed up factoring of survivors.

   The level is 1, 2, or 3. The data type for the position where the update
   hits within a bucket is 16, 24, and 32 bits wide, resp.

   When LOG_BUCKET_REGION > 16, however, the position will have a few bits
   more, so the types for levels 1 and 3 must be changed accordingly.
   This creates, of course, a large memory overhead. */
 
template<int LEVEL> struct bucket_update_size_per_level;
#if LOG_BUCKET_REGION <= 16
template<> struct bucket_update_size_per_level<1> { typedef uint16_t type; };
/* TODO: create a fake 24-bit type as uint8_t[3]. */
template<> struct bucket_update_size_per_level<2> { typedef uint32_t type; };
template<> struct bucket_update_size_per_level<3> { typedef uint32_t type; };
#else
template<> struct bucket_update_size_per_level<1> { typedef uint32_t type; };
template<> struct bucket_update_size_per_level<2> { typedef uint32_t type; };
template<> struct bucket_update_size_per_level<3> { typedef uint64_t type; };
#endif

template <int LEVEL, typename HINT> class bucket_update_t;

#define bu_explicit(LEVEL, HINT, ALIGNMENT_ATTRIBUTE)		\
    template <>								\
    class bucket_update_t<LEVEL, HINT> : public HINT {			\
    public:								\
      typedef HINT hint_type;                                           \
      static inline int level() { return LEVEL; }                       \
      typedef typename bucket_update_size_per_level<LEVEL>::type br_index_t;\
      br_index_t x;							\
      bucket_update_t(){};						\
      bucket_update_t(const uint64_t _x, const fbprime_t p,		\
        const slice_offset_t slice_offset, const slice_index_t slice_index)     \
        : HINT(p, slice_offset, slice_index),				\
          x(limit_cast<br_index_t>(_x))					\
        {}								\
    } ALIGNMENT_ATTRIBUTE


/* it's admittedly somewhat unsatisfactory. I wish I could find a better
 * way. Maybe with alignas ? Is it a trick I can play at template scope
 * with CRTP or so ?
 *
 * (see also https://gcc.gnu.org/bugzilla/show_bug.cgi?id=48138)
 */
bu_explicit(1, shorthint_t, ATTR_ALIGNED(4));
static_assert(sizeof(bucket_update_t<1, shorthint_t>) == 4, "wrong size");

/* 4-byte slice index, 2-byte slice offset, and then 2 bytes of position
 * (bucket_update_size_per_level<1> is 16-bits).
 */
bu_explicit(1, longhint_t, ATTR_ALIGNED(8));
static_assert(sizeof(bucket_update_t<1, longhint_t>) == 8, "wrong size");

bu_explicit(1, primehint_t, ATTR_ALIGNED(8));
static_assert(sizeof(bucket_update_t<1, primehint_t>) == 8, "wrong size");

bu_explicit(2, shorthint_t, ATTR_ALIGNED(8));
static_assert(sizeof(bucket_update_t<2, shorthint_t>) == 8, "wrong size");

bu_explicit(2, longhint_t, ATTR_ALIGNED(16));
static_assert(sizeof(bucket_update_t<2, longhint_t>) == 16, "wrong size");

bu_explicit(3, shorthint_t, ATTR_ALIGNED(8));
static_assert(sizeof(bucket_update_t<3, shorthint_t>) == 8, "wrong size");

/******** Bucket array typedef **************/
/******** Bucket array implementation **************/
template <int LEVEL, typename HINT>
class bucket_array_t : private NonCopyable {
    public:
  static const int level = LEVEL;
  typedef bucket_update_t<LEVEL, HINT> update_t;
  typedef HINT hint_type;
    private:
  update_t *big_data = 0;
  size_t big_size = 0;                  // size of bucket update memory

  update_t ** bucket_write = 0;         // Contains pointers to first empty
                                        // location in each bucket
  update_t ** bucket_start = 0;         // Contains pointers to beginning of
                                        // buckets
  update_t ** bucket_read = 0;          // Contains pointers to first unread
                                        // location in each bucket
  slice_index_t * slice_index = 0;      // For each slice that gets sieved,
                                        // new index is added here
  update_t ** slice_start = 0;          // For each slice there are
                                        // n_bucket pointers, each
                                        // pointer tells where in the
                                        // corresponding bucket the
                                        // updates from that slice start
public:
  uint32_t n_bucket = 0;                // Number of buckets
private:
  size_t   size_b_align = 0;            // cacheline-aligned room for a
                                        // set of n_bucket pointers

  size_t   nr_slices = 0;               // Number of different slices
  size_t   alloc_slices = 0;            // number of checkpoints (each of size
                                        // size_b_align) we have allocated

  static const slice_index_t initial_slice_alloc = 256;
  static const slice_index_t increase_slice_alloc = 128;

  /* Get a pointer to the pointer-set for the i_slice-th slice */
  update_t ** get_slice_pointers(const slice_index_t i_slice) const {
    ASSERT_ALWAYS(i_slice < nr_slices);
    ASSERT_ALWAYS(size_b_align % sizeof(update_t *) == 0);
    return (slice_start + i_slice * size_b_align / sizeof(update_t *));
  }
  void realloc_slice_start(size_t);
  void log_this_update (const update_t update, uint64_t offset,
                        uint64_t bucket_number, where_am_I& w) const;
public:
  size_t nb_of_updates(const int i) const {
      ASSERT((uint32_t) i < n_bucket);
      return bucket_write[i] - bucket_start[i];
  }
  size_t room_allocated_for_updates(const int i) const {
      ASSERT((uint32_t) i < n_bucket);
      return bucket_start[i+1] - bucket_start[i];
  }
  /* Constructor sets everything to zero, and does not allocate memory.
     allocate_memory() does all the allocation. */
  bucket_array_t() = default;
  /* Destructor frees memory, if memory was allocated. If it wasn't, it's
     basically a no-op, except for destroying the bucket_array_t itself. */
  ~bucket_array_t();

  /* Lacking a move constructor before C++11, we make a fake one. We use this
     to store the bucket_array_t on the stack in fill_in_buckets(), to remove
     one level of pointer dereferencing.
     This method copies all the fields of other, then sets them to 0 in other.
     Deconstructing the other bucket_array_t after move() is a no-op, as if
     other had been constructed, but never run allocate_memory(). */
  void move(bucket_array_t &other);

  /* Allocate enough memory to be able to store at least _n_bucket buckets,
     each of at least _bucket_size entries. If at least as much memory had
     already been allocated, does not resize it.  */
  void allocate_memory(const uint32_t _n_bucket, const double fill_ratio,
                       int logI,
                       const slice_index_t prealloc_slices = initial_slice_alloc);
  /* Return a begin iterator over the update_t entries in i_bucket-th
     bucket, generated by the i_slice-th slice */
  const update_t *begin(const size_t i_bucket, const slice_index_t i_slice) const {
    ASSERT_ALWAYS(i_slice < nr_slices);
    const update_t * const p = get_slice_pointers(i_slice)[i_bucket];
    /* The first slice we wrote must start at the bucket start */
    ASSERT_ALWAYS(i_slice != 0 || p == bucket_start[i_bucket]);
    return p;
  }
  /* Return an end iterator over the update_t entries in i_bucket-th
     bucket, generated by the i_slice-th slice */
  const update_t *end(const size_t i_bucket, const slice_index_t i_slice) const {
    ASSERT_ALWAYS(i_slice < nr_slices);
    return (i_slice + 1 < nr_slices) ? get_slice_pointers(i_slice + 1)[i_bucket] :
      bucket_write[i_bucket];
  }

  void reset_pointers();
  slice_index_t get_nr_slices() const {return nr_slices;}
  slice_index_t get_slice_index(const slice_index_t i_slice) const {
    ASSERT_ALWAYS(i_slice < nr_slices);
    return slice_index[i_slice];
  }
  void add_slice_index(const slice_index_t new_slice_index) {
    /* Write new set of pointers for the new factor base slice */
    ASSERT_ALWAYS(nr_slices <= alloc_slices);
    if (nr_slices == alloc_slices) {
      /* We're out of allocated space for the checkpoints. Realloc to bigger
         size. We add space for increase_slice_alloc additional entries. */
      realloc_slice_start(increase_slice_alloc);
    }
    aligned_medium_memcpy((uint8_t *)slice_start + size_b_align * nr_slices, bucket_write, size_b_align);
    slice_index[nr_slices++] = new_slice_index;
  }
  double max_full (unsigned int * fullest_index = NULL) const;
  double average_full () const;
  /* Push an update to the designated bucket. Also check for overflow, if
     SAFE_BUCKETS is defined. */
  void push_update(const int i, const update_t &update) {
#ifdef SAFE_BUCKETS
      if (bucket_start[i] + bucket_size == bucket_write[i]) {
          fprintf(stderr, "# Warning: hit end of bucket nb %d\n", i);
          ASSERT_ALWAYS(0);
          return;
      }
#endif
      *bucket_write[i]++ = update;
  }
  /* Create an update for a hit at location offset and push it to the
     coresponding bucket */
  void push_update(const uint64_t offset, const fbprime_t p,
      const slice_offset_t slice_offset, const slice_index_t slice_index,
      where_am_I& w MAYBE_UNUSED)
  {
      int logB = LOG_BUCKET_REGIONS[LEVEL];
    const uint64_t bucket_number = offset >> logB;
    ASSERT_EXPENSIVE(bucket_number < n_bucket);
    update_t update(offset & ((UINT64_C(1) << logB) - 1), p, slice_offset, slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);
#if defined(TRACE_K)
    log_this_update(update, offset, bucket_number, w);
#endif
    push_update(bucket_number, update);
  }
  void diagnosis(int side, int idx, fb_factorbase::slicing const & fbs) const;
};

/* Downsort sorts the updates in the bucket_index-th bucket of a level-n
   bucket array into a level n-1 bucket array. The update type of the level n
   bucket array can be short or long hint; the level n-1 bucket array is
   always longhint. */
template <int INPUT_LEVEL>
void
downsort(bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
         const bucket_array_t<INPUT_LEVEL, shorthint_t> &BA_in,
         uint32_t bucket_index, where_am_I& w);

template <int INPUT_LEVEL>
void
downsort(bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
         const bucket_array_t<INPUT_LEVEL, longhint_t> &BA_in,
         uint32_t bucket_index, where_am_I& w);

/* A class that stores updates in a single "bucket".
   It's really just a container class with pre-allocated array for storage,
   a persistent read and write pointer. */
template <int LEVEL, typename HINT>
class bucket_single {
  typedef bucket_update_t<LEVEL, HINT> update_t;
  update_t *start; /* start is a "strong" reference */
  update_t *read;  /* read and write are "weak" references into the allocated memory */
  update_t *write;
  size_t _size;
public:
  bucket_single (const size_t size) : _size(size)
  {
    start = electric_new<update_t>(size);
    // start = new update_t[size];
    read = start;
    write = start;
  }

  ~bucket_single() {
    electric_delete(start,_size);
    // delete[] start;
    start = read  = write = NULL;
    _size = 0;
  }

  /* A few of the standard container methods */
  const update_t * begin() const {return start;}
  const update_t * end() const {return write;}
  size_t size() const {return write - start;}

  /* Main writing function: appends update to bucket number i.
   * If SAFE_BUCKETS is not #defined, then there is no checking that there is
   * enough room for the update. This could lead to a segfault, with the
   * current implementation!
   */
  void push_update (const update_t &update)
  {
#ifdef SAFE_BUCKETS
      if (start + _size <= write) {
          fprintf(stderr, "# Warning: hit end of bucket\n");
          ASSERT_ALWAYS(0);
          write--;
      }
#endif
      *(write++) = update;
  }
  const update_t &get_next_update () {
#ifdef SAFE_BUCKETS
    ASSERT_ALWAYS (read < write);
#endif
    return *read++; 
  }
  void rewind_by_1() {if (read > start) read--;}
  bool is_end() const { return read == write; }

  void sort ();
};

/* Stores info containing the complete prime instead of only the low 16 bits */
class bucket_primes_t : public bucket_single<1, primehint_t> {
  typedef bucket_single<1, primehint_t> super;
public:  
  bucket_primes_t (const size_t size) : super(size){}
  ~bucket_primes_t(){}
  void purge (const bucket_array_t<1, shorthint_t> &BA, 
          int i, fb_factorbase::slicing const & fb, const unsigned char *S);
};

/* Stores info containing both slice index and offset instead of only the offset */
class bucket_array_complete : public bucket_single<1, longhint_t> {
    typedef bucket_single<1, longhint_t> super;
public:  
  bucket_array_complete (const size_t size) : super(size){}
  ~bucket_array_complete(){}
  template <typename HINT>
  void purge (const bucket_array_t<1, HINT> &BA, int i, const unsigned char *S);
  template <typename HINT>
  void purge (const bucket_array_t<1, HINT> &BA, int i,
              const unsigned char *S,
              const std::vector<typename bucket_update_t<1, HINT>::br_index_t> &survivors);
private:
#ifdef HAVE_SSE2
  template <typename HINT, int SIZE>
  void purge_1 (
      const bucket_array_t<1, HINT> &BA, const int i,
      const std::vector<typename bucket_update_t<1, HINT>::br_index_t> &survivors);
#endif
};


/* Compute a checksum over the bucket region.

   We import the bucket region into an mpz_t and take it modulo
   checksum_prime. The checksums for different bucket regions are added up,
   modulo checksum_prime. This makes the combined checksum independent of
   the order in which buckets are processed, but it is dependent on size of
   the bucket region. Note that the selection of the sieve region, i.e., of J
   depends somewhat on the number of threads, as we want an equal number of
   bucket regions per thread. Thus the checksums are not necessarily
   clonable between runs with different numbers of threads. */

class sieve_checksum {
  static const unsigned int checksum_prime = 4294967291u; /* < 2^32 */
  unsigned int checksum;
  void update(const unsigned int);

  public:
  sieve_checksum() : checksum(0) {}
  unsigned int get_checksum() {return checksum;}

  /* Combine two checksums */ 
  void update(const sieve_checksum &other) {
    /* Simply (checksum+checksum2) % checksum_prime, but using
       ularith_addmod_ul_ul() to handle sums >= 2^32 correctly. */
    this->update(other.checksum);
  }
  /* Update checksum with the pointed-to data */
  void update(const unsigned char *, size_t);
};

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
    bkmult_specifier(double x) : base(x) {}
    bkmult_specifier(const char * specifier);
    std::string print_all() const;
};

struct buckets_are_full : public clonable_exception {
    struct callback_base {
        virtual void diagnosis(int, fb_factorbase::slicing const &) const = 0;
    };
    callback_base const * base;
    int side;
    bkmult_specifier::key_type key;
    int bucket_number;
    int reached_size;
    int theoretical_max_size;
    std::string message;
    buckets_are_full(callback_base const *, int side, bkmult_specifier::key_type const&, int b, int r, int t);
    virtual const char * what() const noexcept { return message.c_str(); }
    bool operator<(buckets_are_full const& o) const {
        return (double) reached_size / theoretical_max_size < (double) o.reached_size / o.theoretical_max_size;
    }
    virtual clonable_exception * clone() const { return new buckets_are_full(*this); }
    void diagnosis(std::array<fb_factorbase::slicing const *, 2> fbs) const {
        base->diagnosis(side, *fbs[side]);
    }
};

#endif	/* BUCKET_HPP_ */
