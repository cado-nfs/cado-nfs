#ifndef CADO_UTILS_NUMBERTHEORY_FWD_TYPES_HPP
#define CADO_UTILS_NUMBERTHEORY_FWD_TYPES_HPP

#include <exception>

class number_field;
class number_field_order;
class number_field_order_element;
class number_field_element;
class number_field_fractional_ideal;
class number_field_prime_ideal;

struct element_not_integral : public std::exception {};
struct number_field_inconsistency : public std::exception {};

#endif	/* UTILS_NUMBERTHEORY_FWD_TYPES_HPP_ */
