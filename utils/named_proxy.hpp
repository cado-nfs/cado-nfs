#ifndef UTILS_NAMED_PROXY_HPP_
#define UTILS_NAMED_PROXY_HPP_

#include <string>
#include <array>
#include <type_traits>
#include <utility>

namespace cado {
    template<typename T>
        class named_proxy {
            template<typename U> friend class named_proxy;
            static_assert(std::is_reference_v<T>, "T must be a reference");
            using V = std::remove_reference_t<T>;
            using nc = named_proxy<std::remove_const_t<V> &>;
            static constexpr int nvars = V::number_of_variables;
            public:
            T c;
            private:
            std::array<std::string, nvars> X;
            public:
            std::string x() const { return X[0]; }
            std::string y() const requires (nvars > 1) { return X[1]; }
            named_proxy(T c, std::string x) requires (nvars == 1)
                : c(c), X { std::move(x) }
            {}
            named_proxy(T c, std::string x, std::string y) requires (nvars == 2)
                : c(c), X { std::move(x), std::move(y) }
            {}
            /* Don't try to make this one explicit, it won't work. We really
             * need to be able to create a const named proxy from a non-const
             * one.
             *
             * And incidentally, we also need to have a copy ctor, therefore
             * it makes no sense to require T be a const reference: we want
             * this function in both cases.
             */
            named_proxy(nc const & c)
                : c(c.c), X(c.X) {}
        };
} /* namespace cado */

#endif	/* UTILS_NAMED_PROXY_HPP_ */
