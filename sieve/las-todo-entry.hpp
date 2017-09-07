#ifndef LAS_TODO_ENTRY_HPP_
#define LAS_TODO_ENTRY_HPP_

#include <vector>
#include <stdint.h>
#include <gmp.h>
#include "las-base.hpp"
#include "cxx_mpz.hpp"

struct las_todo_entry {
  cxx_mpz p; /* this is the 'special-q', despite the 'p' name... */
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    cxx_mpz r;
    int side;
    bool prime_sq;
    std::vector<uint64_t> prime_factors;

    /* some fields which are specific to the descent */
    int depth;
    int iteration;      /* number of times we failed on this prime */

    /* Default ctor is just to enable sieve_info default ctor */
    las_todo_entry() : side(0), depth(0), iteration(0) {}
    /* Empty p, r is used for "closing brace" */
    las_todo_entry(const int _side, const int _depth) : side(_side), depth(_depth), iteration(0) { }
    las_todo_entry(const mpz_t _p, const mpz_t _r, const int _side, const int _depth = 0, const int _iteration = 0) {
      set(_p, _r, _side, _depth, _iteration);
    }

    // Given a *prime* ell, check whether ell is coprime to current
    // entry.
    bool is_coprime_to(unsigned long ell) const {
        if (prime_sq)
            return (mpz_cmp_ui(p, ell) != 0);
        else {
            for (unsigned int i = 0; i < prime_factors.size(); ++i) {
                if (prime_factors[i] == ell)
                    return false;
            }
            return true;
        }
    }
private:
  // TODO: this function becomes a bit long for a header file.
  // But for the moment we don't have a las-todo.cpp, so...
  void set(const mpz_t _p, const mpz_t _r, const int _side, const int _depth, const int _iteration);
};

#endif	/* LAS_TODO_ENTRY_HPP_ */

