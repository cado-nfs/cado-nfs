#ifndef LAS_BASE_HPP_
#define LAS_BASE_HPP_

/* C++ base classes that are not application-specific. */

/* Base class with private copy-constructor and assignment operator.
   Classes which are not copy-constructible can inherit this with:
   private NonCopyable */
class NonCopyable {
 protected:
   NonCopyable() {}
   ~NonCopyable() {}
 private:
   NonCopyable(const NonCopyable&);
   NonCopyable& operator=(const NonCopyable&);
};

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
      memset(this,  0, sizeof(T));
    }

    _padded_pod(const _padded_pod &x) {
      memcpy(this, &x, sizeof(T));
    }

    const _padded_pod &operator=(const _padded_pod &x) {
      memcpy(this, &x, sizeof(T));
      return *this;
    }
};

#endif
