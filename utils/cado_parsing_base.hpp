#ifndef UTILS_CADO_PARSING_BASE_HPP_
#define UTILS_CADO_PARSING_BASE_HPP_

#include <exception>

namespace cado {

struct parse_error: public std::exception {
    const char * what() const noexcept override { return "parse error"; }
};
struct token_error: public std::exception {
    const char * what() const noexcept override { return "token error"; }
};

} /* namespace cado */


#endif	/* UTILS_CADO_PARSING_BASE_HPP_ */
