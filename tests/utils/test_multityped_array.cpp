#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdio>

#include <array>
#include <iostream>
#include <type_traits>

#include "multityped_array.hpp"

template <int n> struct type_factory {
    using type = std::array<int, n>;
};

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
struct base {
    virtual size_t size() const = 0;
    virtual ~base() = default;
};
template <int n>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
struct type_hierarchy_factory
    : base
    , std::array<int, n> {
    using type = type_hierarchy_factory<n>;
    size_t size() const override { return n; }
    type_hierarchy_factory() { std::array<int, n>::fill(0); }
    ~type_hierarchy_factory() override = default;
};
template <int m, int n>
// NOLINTNEXTLINE(cert-dcl58-cpp)
struct std::common_type<type_hierarchy_factory<m>, type_hierarchy_factory<n>> {
    using type = base;
};

int main()
{
    using cado::multityped_array;

    using A_t = multityped_array<type_factory, 0, 6>;

    A_t A;

    /* testing the get<> accessors */
    auto A3 = A.get<3>();
    std::cout << "A3[2] == " << A3[2] << "\n";
    std::cout << "A3[2] == " << cado::get<3>(A)[2] << "\n";

    /* two tests of the foreach function */
    A.foreach([](auto & x) {
        for (unsigned int i = 0; i < x.size(); ++i)
            x[i] = i * x.size();
    });

    A.foreach([](auto const & x) {
        constexpr size_t n =
            std::tuple_size_v<std::remove_reference_t<decltype(x)>>;
        std::cout << "field has " << n << " elements";
        for (auto a: x)
            std::cout << " " << a;
        std::cout << "\n";
    });

    /* This additional test uses external state */
    size_t sum = 0;
    A.foreach([&](auto const & x) { sum += x.size(); });
    std::cout << "total " << sum << "\n";

    /* testing search. This only works if we have a common type among the
     * members!
     */
    using B_t = multityped_array<type_hierarchy_factory, 0, 6>;
    B_t B;

    B.foreach([](auto & x) {
        for (unsigned int i = 0; i < x.size(); ++i)
            x[i] = i * x.size();
    });

    auto item_containing_nth = [](auto & B, size_t v) {
        auto locate = [&](auto const & x) -> const base * {
                if (v < x.size())
                    return &x;
                v -= x.size();
                return nullptr;
                };
        return B.find(locate);
    };
    auto const * p = item_containing_nth(B, 5);

    if (p) {
        printf("Found match in the container of size %zu\n", p->size());
    } else {
        printf("no match\n");
    }

    A.foreach([](auto & x) {
        for (unsigned int i = 0; i < x.size(); ++i)
            x[i] = i * x.size();
    });

    auto nth = [](auto & B, size_t v) {
        auto locate = [&](auto & x) -> int const * {
            if (v < x.size())
                return &(x[v]);
            v -= x.size();
            return nullptr;
        };
        return B.find(locate);
    };

    for(size_t i = 0 ; i < sum ; i++) {
        auto const * q = nth(B, i);
        if (q) {
            printf("%zu-th entry in the container is %d\n", i, *q);
        } else {
            printf("no match\n");
        }
    }
}
