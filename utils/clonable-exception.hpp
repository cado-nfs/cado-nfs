#ifndef CLONABLE_EXCEPTION_HPP_
#define CLONABLE_EXCEPTION_HPP_

#include <exception>

struct clonable_exception : public std::exception {
    virtual clonable_exception * clone() const = 0;
};


#endif	/* CLONABLE_EXCEPTION_HPP_ */
