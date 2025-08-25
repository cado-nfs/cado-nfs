#ifndef UTILS_NAMED_PROXY_HPP_
#define UTILS_NAMED_PROXY_HPP_

#include <string>
#include <type_traits>
#include <utility>

namespace cado {
    template<typename T>
        class named_proxy {
            static_assert(std::is_reference_v<T>, "T must be a reference");
            using V = std::remove_reference_t<T>;
            using nc = named_proxy<std::remove_const_t<V> &>;
            static constexpr const bool is_c = std::is_const_v<V>;
            public:
            T c;
            std::string x;
            named_proxy(T c, std::string x)
                : c(c), x(std::move(x))
            {}
            template<typename U = T>
                explicit named_proxy(U const & c)
                requires std::is_same_v<U, nc>
                : c(c.c), x(c.x) {}
        };
} /* namespace cado */


#endif	/* UTILS_NAMED_PROXY_HPP_ */
