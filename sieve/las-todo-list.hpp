#ifndef CADO_LAS_TODO_LIST_HPP
#define CADO_LAS_TODO_LIST_HPP

#include <cstdint>
#include <cstdio>
#include <climits>

#include <fstream>
#include <vector>
#include <memory>
#include <mutex>
#include <stack>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "cado_poly.h"
#include "special-q.hpp"
#include "galois_action.hpp"

struct cxx_param_list;

class las_todo_list : private std::stack<special_q> {
    std::mutex mm;

    /* We use these flags to decide whether some extra work may appear or
     * not. (in the descent case)
     */
    public:
    /* the "created" data member is (almost) only used internally in
     * push_unlocked, and checked in feed_qrange and feed_qlist (to
     * compare against nq_max).
     */
    size_t created = 0;
    // size_t pulled = 0;
    // size_t done = 0;
    size_t nq_max = SIZE_MAX;

    bool print_todo_list_flag = false;
    private:


    cxx_cado_poly cpoly;
    cxx_gmp_randstate rstate;

    typedef std::stack<special_q> super;

    int random_sampling = 0;
    cxx_mpz q0;
    cxx_mpz q1;
    galois_action galois;        /* Used to skip some primes */
    std::unique_ptr<std::ifstream> todo_list_fd;
    bool feed_qrange(gmp_randstate_t);
    bool feed_qlist();

    void push_unlocked(special_q const & q, special_q const & = {})
    {
        super::emplace(q);
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
    void push(special_q const & q, special_q const & parent = {})
    {
        const std::lock_guard<std::mutex> foo(mm);
        push_unlocked(q, parent);
    }
    bool empty() {
        const std::lock_guard<std::mutex> foo(mm);
        return super::empty();
    }
    private:
    std::vector<uint64_t> next_legitimate_specialq(cxx_mpz & r, cxx_mpz const & s, unsigned long diff) const;
    public:

    /* }}} */

    special_q feed_and_pop();

    public:

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
};


#endif	/* CADO_LAS_TODO_LIST_HPP */
