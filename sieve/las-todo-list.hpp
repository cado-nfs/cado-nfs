#ifndef CADO_LAS_TODO_LIST_HPP
#define CADO_LAS_TODO_LIST_HPP

#include <cstdint>            // for uint64_t, UINT64_MAX
#include <cstdio>             // for FILE, NULL, size_t
#include <climits>

#include <algorithm>
#include <list>                // for list
#include <memory>
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
    unsigned int nq_max = UINT_MAX;
    int random_sampling = 0;
    cxx_mpz q0;
    cxx_mpz q1;
    const char * galois = nullptr;        /* Used to skip some primes */
    std::unique_ptr<std::ifstream> todo_list_fd;
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
    int sqside = -1;
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
        const std::lock_guard<std::mutex> foo(mm);
        push_withdepth_unlocked(p, r, side, depth, iteration);
    }
    void push(cxx_mpz const & p, cxx_mpz const & r, int side)
    {
        push_withdepth(p, r, side, 0);
    }
    void push_closing_brace(int depth)
    {
        const std::lock_guard<std::mutex> foo(mm);
        super::push(las_todo_entry(-1, depth));
    }
    las_todo_entry pop()
    {
        const std::lock_guard<std::mutex> foo(mm);
        las_todo_entry r = super::top();
        super::pop();
        return r;
    }

    static int is_closing_brace(las_todo_entry const & doing)
    {
        return doing.side < 0;
    }
    /* }}} */

    bool is_random() const { return random_sampling != 0; }

    bool feed(gmp_randstate_t rstate);
    las_todo_entry * feed_and_pop(gmp_randstate_t rstate);
    las_todo_list(las_todo_list const &) = delete;
    las_todo_list(las_todo_list &&) = delete;
    las_todo_list& operator=(las_todo_list const &) = delete;
    las_todo_list& operator=(las_todo_list &&) = delete;
    ~las_todo_list() = default;

    las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl);

    super save() {
        const std::lock_guard<std::mutex> foo(mm);
        return *this; /* NOLINT(cppcoreguidelines-slicing) */
    }
    void restore(super && x) {
        const std::lock_guard<std::mutex> foo(mm);
        std::swap((super&)*this, x);
    }

    static void configure_switches(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);

    void print_todo_list(cxx_param_list & pl, gmp_randstate_ptr, int nthreads = 1) const;
};


#endif	/* CADO_LAS_TODO_LIST_HPP */
