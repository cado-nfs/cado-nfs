#ifndef SIEVE_SIEVE_METHODS_HPP_
#define SIEVE_SIEVE_METHODS_HPP_

#include <concepts>
#include <type_traits>

#include "las-plattice.hpp"
#include "las-qlattice.hpp"
#include "las-smallsieve.hpp"
#include "las-todo-list.hpp"
#include "smallsieve.hpp"
#include "siqs-smallsieve.hpp"
#include "siqs-largesieve.hpp"

template<typename T>
concept todo_list_class =
        std::derived_from<T, todo_list_base>
        && !std::is_abstract_v<T>
        && std::is_constructible_v<T, cxx_cado_poly const &, cxx_param_list &>;

template<typename T>
concept special_q_data_class =
        std::derived_from<T, special_q_data_base>
        && !std::is_abstract_v<T>
        && std::is_default_constructible_v<T>;

template<typename T>
concept small_sieve_data_class =
        std::derived_from<T, small_sieve_data>
        && !std::is_abstract_v<T>
        && std::is_default_constructible_v<T>;

template<typename T>
concept sieve_method = std::is_empty_v<T>
    && std::is_trivially_default_constructible_v<T> /* for tag dispatch */
    && todo_list_class<typename T::todo_list>
    && special_q_data_class<typename T::special_q_data>
    && small_sieve_data_class<typename T::smallsieve>
    && requires { typename T::largesieve; };

struct NFS {
  using todo_list = las_todo_list;
  using special_q_data = qlattice_basis;
  using largesieve = plattice_enumerator;
  using smallsieve = las_small_sieve_data;
};

struct SIQS {
  using todo_list = siqs_todo_list;
  using special_q_data = siqs_special_q_data;
  using largesieve = siqs_largesieve;
  using smallsieve = siqs_small_sieve_data;
};

static_assert(sieve_method<NFS>);
static_assert(sieve_method<SIQS>);

#endif /* SIEVE_SIEVE_METHODS_HPP_ */
