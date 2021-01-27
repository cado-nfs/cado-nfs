#ifndef LAS_TODO_LIST_HPP_
#define LAS_TODO_LIST_HPP_

#include <algorithm>
#include <cstdint>            // for uint64_t, UINT64_MAX
#include <cstdio>             // for FILE, NULL, size_t
#include <list>                // for list
#include <mutex>               // for mutex, lock_guard
#include <stack>               // for swap, stack
#include <gmp.h>               // for gmp_randstate_t
#include "cxx_mpz.hpp"
#include "cado_poly.h"
#include "las-todo-entry.hpp"  // for las_todo_entry
struct cxx_param_list;

class las_todo_list : private std::stack<las_todo_entry> {
    std::mutex mm;
    cxx_cado_poly cpoly;
    typedef std::stack<las_todo_entry> super;
    /* "history" is append-only: everything we pop from the stack goes
     * here, and lives until the destruction */
    std::list<las_todo_entry> history;
    unsigned int nq_max = 0;
    int random_sampling = 0;
    cxx_mpz q0;
    cxx_mpz q1;
    const char * galois;        /* Used to skip some primes */
    FILE * todo_list_fd = NULL;
    bool feed_qrange(gmp_randstate_t);
    bool feed_qlist();
    void push_withdepth_unlocked(cxx_mpz const & p, cxx_mpz const & r, int side, int depth, int iteration = 0)
    {
        super::push(las_todo_entry(p, r, side, depth, iteration));
    }
    void push_unlocked(cxx_mpz const & p, cxx_mpz const & r, int side)
    {
        push_withdepth_unlocked(p, r, side, 0);
    }
    public:
    int sqside;
    /* For composite special-q: note present both in las_info and
     * las_todo_list */
    bool allow_composite_q = false;
    bool print_todo_list_flag = false;
    uint64_t qfac_min = 1024;
    uint64_t qfac_max = UINT64_MAX;

    unsigned int nq_pushed = 0;

    /*{{{*/
    size_t size() const { return super::size(); }
    void push_withdepth(cxx_mpz const & p, cxx_mpz const & r, int side, int depth, int iteration = 0)
    {
        std::lock_guard<std::mutex> foo(mm);
        push_withdepth_unlocked(p, r, side, depth, iteration);
    }
    void push(cxx_mpz const & p, cxx_mpz const & r, int side)
    {
        push_withdepth(p, r, side, 0);
    }
    void push_closing_brace(int depth)
    {
        std::lock_guard<std::mutex> foo(mm);
        super::push(las_todo_entry(-1, depth));
    }
    las_todo_entry pop()
    {
        std::lock_guard<std::mutex> foo(mm);
        las_todo_entry r = super::top();
        super::pop();
        return r;
    }

    int is_closing_brace(las_todo_entry const & doing) const
    {
        return doing.side < 0;
    }
    /* }}} */

    bool is_in_qfac_range(uint64_t p) const {
        return (p >= qfac_min) && (p >= qfac_max);
    }

    bool is_random() const { return random_sampling != 0; }

    bool feed(gmp_randstate_t rstate);
    las_todo_entry * feed_and_pop(gmp_randstate_t rstate);

    las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl);
    ~las_todo_list();

    super save() {
        std::lock_guard<std::mutex> foo(mm);
        return (super)*this;
    }
    void restore(super && x) {
        std::lock_guard<std::mutex> foo(mm);
        std::swap((super&)*this, x);
    }

    static void configure_switches(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);

    void print_todo_list(cxx_param_list & pl, gmp_randstate_ptr, int nthreads = 1) const;
};


#endif	/* LAS_TODO_LIST_HPP_ */
