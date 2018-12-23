#ifndef LAS_BASE_HPP_
#define LAS_BASE_HPP_

#include <string.h>
#include "utils_cxx.hpp"

/* C++ base classes that are not application-specific. */

// If a plain-old data type T inherits from _padded_pod<T>, it ensures that all
// the sizeof(T) bytes of any instance of T will be initialized, even if its
// data members do not occupy the whole memory region.
// For instance, this allows one to write `fwrite(&x, sizeof(T), 1, f)' without
// Valgrind complaining because of uninitialized memory reads.
//
// IT IS NOT VALID TO USE THIS IF T HAS NON-POD DATA MEMBERS !!!
template <typename T>
class _padded_pod {
  public:
    _padded_pod() {
      // see bug 21663
      memset((void*)this,  0, sizeof(T));
    }

    _padded_pod(const _padded_pod &x) {
      // see bug 21663
      memcpy((void*)this, (const void *) &x, sizeof(T));
    }

    const _padded_pod &operator=(const _padded_pod &x) {
      // see bug 21663
      memcpy((void*)this, &x, sizeof(T));
      return *this;
    }
};

#endif
