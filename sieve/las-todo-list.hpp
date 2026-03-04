#ifndef CADO_LAS_TODO_LIST_HPP
#define CADO_LAS_TODO_LIST_HPP

#include <cstdint>
#include <cstdio>

#include <fstream>
#include <memory>
#include <mutex>
#include <stack>
#include <string>
#include <vector>

#include <gmp.h>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "cado_poly.h"
#include "special-q.hpp"
#include "galois_action.hpp"

struct cxx_param_list;

class todo_list_base : protected std::stack<special_q> {
    protected:
    using super = std::stack<special_q>;

    public:
    todo_list_base(cxx_cado_poly const & cpoly, cxx_param_list & pl);

    /* no copy nor move ctor and assigment operator */
    todo_list_base(todo_list_base const &) = delete;
    todo_list_base(todo_list_base &&) = delete;
    todo_list_base& operator=(todo_list_base const &) = delete;
    todo_list_base& operator=(todo_list_base &&) = delete;

    virtual ~todo_list_base() = default;

    /* create a unique_ptr to an object of a derived class */
    template<class derived, class... Args>
      requires std::derived_from<derived, todo_list_base>
               && std::is_constructible_v<derived, Args...>
    static std::unique_ptr<todo_list_base> create(Args&&... args)
    {
        todo_list_base *ptr = new derived(std::forward<Args>(args)...);
        return std::unique_ptr<todo_list_base>(ptr);
    }

    static void configure_switches(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);

    special_q feed_and_pop();
    using super::size;
    void push(special_q const & q, special_q const & parent = {})
    {
        const std::lock_guard<std::mutex> foo(mm);
        push_unlocked(q, parent);
    }
    bool empty()
    {
        const std::lock_guard<std::mutex> foo(mm);
        return super::empty();
    }

    void print_todo_list(cxx_param_list const & pl, int nthreads = 1) const;

    bool is_in_qfac_range(uint64_t p) const
    {
        return (p >= qfac_min) && (p <= qfac_max);
    }

    bool allow_composite_special_q() const
    {
        return allow_composite_q;
    }

    bool should_todo_list_be_printed() const
    {
        return print_todo_list_flag;
    }

    size_t ncreated() const
    {
        return created;
    }

    friend std::ostream& operator<<(
            std::ostream& os,
            todo_list_base const & todo)
    {
        return os << todo.string_repr();
    }

    protected:
    virtual bool feed_qrange(gmp_randstate_t) = 0;
    virtual std::unique_ptr<todo_list_base> create_sub_todo_list(
            unsigned int i,
            int nthreads,
            cxx_param_list const & pl) const = 0;
    virtual std::string string_repr() const = 0;
    bool feed_qlist();

    void push_unlocked(special_q const & q, special_q const & = {})
    {
        super::emplace(q);
        created++;
    }

    protected:
    cxx_cado_poly cpoly;
    cxx_gmp_randstate rstate;

    int random_sampling = 0;
    galois_action galois;        /* Used to skip some primes */
    size_t nq_max = SIZE_MAX;
    std::unique_ptr<std::ifstream> todo_list_fd;

    /* for the most part, the sqside shouldn't be considered relevant,
     * since the todo list can provide special-q's on either side.
     */
    int sqside = -1;

    /* The "created" data member is only used internally in push_unlocked, and
     * checked in feed_qrange and feed_qlist (to compare against nq_max).
     */
    size_t created = 0;
    bool print_todo_list_flag = false;

    bool allow_composite_q = false;
    uint64_t qfac_min = 1024;
    uint64_t qfac_max = UINT64_MAX;

    private:
    std::mutex mm;
};

namespace fmt {
    template <> struct formatter<todo_list_base>: ostream_formatter {};
}

class las_todo_list : public todo_list_base {
    private:
    cxx_mpz q0;
    cxx_mpz q1;

    std::vector<uint64_t> next_legitimate_specialq(cxx_mpz & r,
                                                   cxx_mpz const & s,
                                                   unsigned long diff) const;
    bool feed_qrange(gmp_randstate_t) override;
    std::unique_ptr<todo_list_base> create_sub_todo_list(unsigned int i,
                                    int nthreads,
                                    cxx_param_list const & pl) const override;
    std::string string_repr() const override;

    public:
    las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl);

    static void declare_usage(cxx_param_list & pl);
};


class siqs_todo_list : public todo_list_base {
    private:
    cxx_mpz qidx0;
    cxx_mpz qidx1;
    uint64_t qfac_nfac = 0U;
    std::vector<cxx_mpz> qfac_primes;

    uint64_t get_qfac_prime(uint64_t i);

    static void k_combination_from_index(std::vector<uint64_t>::reverse_iterator it,
                                         cxx_mpz & idx, uint64_t k);
    special_q special_q_from_index(cxx_mpz const & idx,
                                   cxx_mpz & scratch_q,
                                   std::vector<uint64_t> & scratch_vec);
    bool feed_qrange(gmp_randstate_t) override;
    std::unique_ptr<todo_list_base> create_sub_todo_list(unsigned int i,
                                    int nthreads,
                                    cxx_param_list const & pl) const override;
    std::string string_repr() const override;

    public:
    siqs_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl);

    static void declare_usage(cxx_param_list & pl);
};
#endif	/* CADO_LAS_TODO_LIST_HPP */
