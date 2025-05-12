#ifndef CADO_LAS_TODO_LIST_HPP
#define CADO_LAS_TODO_LIST_HPP

#include <cstdint>
#include <cstdio>
#include <climits>

#include <condition_variable>
#include <fstream>
#include <algorithm>
#include <list>
#include <memory>
#include <mutex>
#include <stack>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "cado_poly.h"
#include "las-todo-entry.hpp"

struct cxx_param_list;

class las_todo_list : private std::stack<las_todo_entry> {
    std::mutex mm;
    std::condition_variable work_to_do;

    /* We use these flags to decide whether some extra work may appear or
     * not. (in the descent case)
     */
    public:
    size_t created = 0;
    size_t pulled = 0;
    size_t done = 0;
    size_t nq_max = SIZE_MAX;

    bool print_todo_list_flag = false;
    private:


    cxx_cado_poly cpoly;
    cxx_gmp_randstate rstate;

    typedef std::stack<las_todo_entry> super;

    int random_sampling = 0;
    cxx_mpz q0;
    cxx_mpz q1;
    const char * galois = nullptr;        /* Used to skip some primes */
    std::unique_ptr<std::ifstream> todo_list_fd;
    bool feed_qrange(gmp_randstate_t);
    bool feed_qlist();

    void push_withdepth_unlocked(cxx_mpz const & p, cxx_mpz const & r, int side, int depth, int iteration = 0)
    {
        super::emplace(p, r, side, depth, iteration);
        created++;
    }
    void push_unlocked(cxx_mpz const & p, cxx_mpz const & r, int side)
    {
        super::emplace(p, r, side);
        created++;
    }

    /* for the most part, the sqside shouldn't be considered relevant,
     * since the todo list can provide special-q's on either side.
     */
    int sqside = -1;
    public:
    /* For composite special-q: note present both in las_info and
     * las_todo_list */
    bool allow_composite_q = false;
    uint64_t qfac_min = 1024;
    uint64_t qfac_max = UINT64_MAX;

    public:
    /*{{{*/
    using super::size;
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
        /* it's very much unclear that we want to keep this.
         */
        const std::lock_guard<std::mutex> foo(mm);
        super::push(las_todo_entry::closing_brace(depth));
    }
    private:
    std::vector<uint64_t> next_legitimate_specialq(cxx_mpz & r, cxx_mpz const & s, unsigned long diff) const;
    public:

    /* }}} */

    las_todo_entry feed_and_pop();

    las_todo_list(las_todo_list const &) = delete;
    las_todo_list(las_todo_list &&) = delete;
    las_todo_list& operator=(las_todo_list const &) = delete;
    las_todo_list& operator=(las_todo_list &&) = delete;
    ~las_todo_list() = default;

    las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl);

    public:
    static void configure_switches(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);

    void print_todo_list(cxx_param_list & pl, int nthreads = 1) const;

    public:
    struct pulled_todo_entry : public las_todo_entry {
            las_todo_list * L = nullptr;
        explicit pulled_todo_entry(las_todo_list * L)
            : las_todo_entry(L->feed_and_pop())
            , L(L)
        { }
        ~pulled_todo_entry() {
            if (L) L->done++;
        }
        pulled_todo_entry() = delete;
        pulled_todo_entry(pulled_todo_entry const & E) = delete;
        pulled_todo_entry& operator=(pulled_todo_entry const & E) = delete;
        pulled_todo_entry(pulled_todo_entry && E) noexcept
            : las_todo_entry((las_todo_entry &&) E)
            , L(E.L)
        {
            E.L = nullptr;
        }
        pulled_todo_entry& operator=(pulled_todo_entry && E) noexcept
        {
            (las_todo_entry &) *this = (las_todo_entry const &) E;
            L = E.L;
            E.L = nullptr;
            return *this;
        }
    };
};


#endif	/* CADO_LAS_TODO_LIST_HPP */
