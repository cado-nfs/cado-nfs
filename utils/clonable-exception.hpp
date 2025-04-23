#ifndef CADO_CLONABLE_EXCEPTION_HPP
#define CADO_CLONABLE_EXCEPTION_HPP

#include <exception>

struct clonable_exception : public std::exception {
    virtual clonable_exception * clone() const = 0;
};


#endif	/* CADO_CLONABLE_EXCEPTION_HPP */
