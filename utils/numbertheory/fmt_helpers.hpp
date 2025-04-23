#ifndef CADO_UTILS_NUMBERTHEORY_FMT_HELPERS_HPP
#define CADO_UTILS_NUMBERTHEORY_FMT_HELPERS_HPP

#include "fmt/base.h"
#include "fmt/format.h"

namespace fmt {
    template<typename T>
    struct fmt_helper_sagemath {
        protected:
            enum { TEXT, SAGEMATH, MACHINE, MAGMA } custom_format = formatter<T>::custom_format_default;
            /* this can be overridden */
            static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        public:
            FMT_CONSTEXPR auto parse(basic_format_parse_context<char>& ctx)
                -> decltype(ctx.begin())
            {
                auto begin = ctx.begin(), end = ctx.end();
                if (begin != end && *begin == 'S') {
                    ++begin;
                    custom_format = SAGEMATH;
                } else if (begin != end && *begin == 'M') {
                    ++begin;
                    custom_format = MACHINE;
                } else if (begin != end && *begin == 'm') {
                    ++begin;
                    custom_format = MAGMA;
                } else if (begin != end && *begin == 'T') {
                    ++begin;
                    custom_format = TEXT;
                }
                return begin;
            }
    };
}
#endif	/* UTILS_NUMBERTHEORY_FMT_HELPERS_HPP_ */
