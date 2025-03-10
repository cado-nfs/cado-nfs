#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <array>
#include <iostream>
#include <utility>

#include "multityped_array.hpp"

template<int n> struct type_factory {
    typedef std::array<int, n> type;
};

struct print {
    template<typename T>
    void operator()(T const & x) {
        constexpr size_t n = std::tuple_size<T>::value;
        std::cout << "field has " << n << " elements";
        for(auto a : x) std::cout << " " << a;
        std::cout << "\n";
    }
};

struct print2 {
    int k;
    template<typename T>
    void operator()(T const & x) {
        constexpr size_t n = std::tuple_size<T>::value;
        std::cout << "field has " << n << " elements";
        for(auto a : x) std::cout << " " << a;
        std::cout << " [" << k << "]";
        std::cout << "\n";
    }
};

struct fill {
    template<typename T>
    void operator()(T & x) {
        for(unsigned int i = 0 ; i < x.size() ; ++i)
            x[i] = i*x.size();
    }
};

struct return_pointer_if_in_subrange {
    typedef int * type;
    typedef int key_type;
    template<typename T>
    type operator()(T & x, int & k) {
        if ((size_t) k < x.size()) {
            return &(x[k]);
        } else {
            k -= x.size();
            return nullptr;
        }
    }
    template<typename T>
    type operator()(T const & x, int & k) {
        if ((size_t) k < x.size()) {
            return &(x[k]);
        } else {
            k -= x.size();
            return nullptr;
        }
    }
};

struct accumulate_sizes {
    template<typename T>
    int operator()(int t0, T const & a) const {
        return t0 + a.size();
    }
};



int main()
{
    typedef multityped_array<type_factory, 0, 6> A_t;
    A_t A;

    multityped_array_foreach(fill(), A);

    type_factory<3>::type A3 = A.get<3>();
    std::cout << "A3[2] == " << A3[2] << "\n";

    multityped_array_foreach(print2 {1}, A);

    const int v = 5;
    multityped_array_locate<return_pointer_if_in_subrange>()(A, v);

    const int s = multityped_array_fold(accumulate_sizes(), 0, A);

    std::cout << "total " << s << "\n";
}
