#ifndef CADO_UTILS_ARITHXX_RESIDUE_STD_OP_HPP
#define CADO_UTILS_ARITHXX_RESIDUE_STD_OP_HPP

#include <cstddef>
#include <cstdint>

#include "macros.h"

/* Functions for residues. Slow and clean
 * This class extends the Residue type of an existing Modulus type by
 * adding standard C++ operators for the modular arithmetic. For this to
 * work, each residue needs a reference to the modulus w.r.t. which it
 * was declared, as there is no other way of passing the modulus to the
 * arithmetic operators.
 */

namespace arithxx_details {

template <typename layer>
class ResidueStdOp : public layer::Residue
{
    using Residue = layer::Residue;
    using Modulus = layer::Modulus;
    using Integer = layer::Integer;

  private:
    Modulus const & m;
    void assertValid(ResidueStdOp const & other MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(&m = &other.m);
    }

  public:
    ResidueStdOp(Modulus const & m) noexcept
        : Residue(m)
        , m(m)
    {
    }
    ResidueStdOp(Modulus const & m, Residue const & s)
        : Residue(m)
        , m(m)
    {
        m.set(*this, s);
    }
    ResidueStdOp(Modulus const & m, Integer const & s)
        : Residue(m)
        , m(m)
    {
        m.set(*this, s);
    }
    ResidueStdOp(Modulus const & m, uint64_t const & s)
        : Residue(m)
        , m(m)
    {
        m.set(*this, s);
    }
    ResidueStdOp(ResidueStdOp const & s)
        : Residue(s.m)
        , m(s.m)
    {
        *this = s;
    }
    ResidueStdOp(ResidueStdOp && s) noexcept
        : Residue(s.m)
        , m(s.m)
    {
        *this = s;
    }

    ResidueStdOp & operator=(ResidueStdOp const & s)
    {
        assertValid(s);
        if (this != &s) 
            m.set(*this, s);
        return *this;
    }
    ResidueStdOp & operator=(uint64_t const s)
    {
        m.set(*this, s);
        return *this;
    }
    ResidueStdOp & operator=(Integer const & s)
    {
        m.set(*this, s);
        return *this;
    }

    bool operator==(ResidueStdOp const & s) const
    {
        assertValid(s);
        return m.equal(*this, s);
    }
    bool operator!=(ResidueStdOp const & s) const
    {
        assertValid(s);
        return !m.equal(*this, s);
    }
    bool is0() const { return m.is0(*this); }
    bool is1() const { return m.is1(*this); }

    operator Integer()
    {
        return m.get(*this);
    }

    ResidueStdOp operator-() const
    {
        ResidueStdOp r = *this;
        m.neg(r, r);
        return r;
    }
    /* Prefix ++ and -- */
    ResidueStdOp & operator++()
    {
        m.add1(*this, *this);
        return *this;
    }
    ResidueStdOp & operator--()
    {
        m.sub1(*this, *this);
        return *this;
    }
    /* Postfix ++ and -- */
    ResidueStdOp operator++(int)
    {
        ResidueStdOp r = *this;
        m.add1(*this, *this);
        return r;
    }
    ResidueStdOp operator--(int)
    {
        ResidueStdOp r = *this;
        m.sub1(*this, *this);
        return r;
    }

    ResidueStdOp operator+(ResidueStdOp const & s) const
    {
        assertValid(s);
        ResidueStdOp r(m);
        m.add(r, *this, s);
        return r;
    }
    ResidueStdOp operator+(uint64_t const s) const
    {
        ResidueStdOp r(m, s);
        m.add(r, *this, r);
        return r;
    }
    ResidueStdOp operator+(Integer const & s) const
    {
        ResidueStdOp r(m, s);
        m.add(r, *this, r);
        return r;
    }
    ResidueStdOp & operator+=(ResidueStdOp const & s)
    {
        assertValid(s);
        m.add(*this, *this, s);
        return *this;
    }
    ResidueStdOp & operator+=(uint64_t const s)
    {
        ResidueStdOp t(m, s);
        m.add(*this, *this, t);
        return *this;
    }
    ResidueStdOp & operator+=(Integer const & s)
    {
        ResidueStdOp t(m, s);
        m.add(*this, *this, t);
        return *this;
    }
    friend ResidueStdOp operator+(uint64_t const s, ResidueStdOp const & t)
    {
        return t + s;
    }
    friend ResidueStdOp operator+(Integer const & s, ResidueStdOp const & t)
    {
        return t + s;
    }
    friend Integer & operator+=(Integer & s, ResidueStdOp const & t)
    {
        ResidueStdOp r(t.m, s);
        t.m.add(r, r, t);
        s = (Integer)r;
        return s;
    }

    ResidueStdOp operator-(ResidueStdOp const & s) const
    {
        assertValid(s);
        ResidueStdOp r(m);
        m.sub(r, *this, s);
        return r;
    }
    ResidueStdOp operator-(uint64_t const s) const
    {
        ResidueStdOp r(m, s);
        m.sub(r, *this, r);
        return r;
    }
    ResidueStdOp operator-(Integer const & s) const
    {
        ResidueStdOp r(m, s);
        m.sub(r, *this, r);
        return r;
    }
    ResidueStdOp & operator-=(ResidueStdOp const & s)
    {
        assertValid(s);
        m.sub(*this, *this, s);
        return *this;
    }
    ResidueStdOp & operator-=(uint64_t const s)
    {
        ResidueStdOp t(m, s);
        m.sub(*this, *this, t);
        return *this;
    }
    ResidueStdOp & operator-=(Integer const & s)
    {
        ResidueStdOp t(m, s);
        m.sub(*this, *this, t);
        return *this;
    }
    friend ResidueStdOp operator-(uint64_t const s, ResidueStdOp const & t)
    {
        return ResidueStdOp(t.m, s) - t;
    }
    friend ResidueStdOp operator-(Integer const & s, ResidueStdOp const & t)
    {
        return ResidueStdOp(t.m, s) - t;
    }
    friend Integer & operator-=(Integer & s, ResidueStdOp const & t)
    {
        return s = (Integer)(ResidueStdOp(t.m, s) - t);
    }

    ResidueStdOp operator*(ResidueStdOp const & s) const
    {
        assertValid(s);
        ResidueStdOp r(m);
        if (CONSTANT_P(this == &s) && this == &s) {
            m.sqr(r, *this);
        } else {
            m.mul(r, *this, s);
        }
        return r;
    }
    ResidueStdOp operator*(uint64_t const s) const
    {
        ResidueStdOp r(m, s);
        m.mul(r, *this, r);
        return r;
    }
    ResidueStdOp operator*(Integer const & s) const
    {
        ResidueStdOp r(m, s);
        m.mul(r, *this, r);
        return r;
    }
    ResidueStdOp & operator*=(ResidueStdOp const & s)
    {
        assertValid(s);
        m.mul(*this, *this, s);
        return *this;
    }
    ResidueStdOp & operator*=(uint64_t const s)
    {
        ResidueStdOp t(m, s);
        m.mul(*this, *this, t);
        return *this;
    }
    ResidueStdOp & operator*=(Integer const & s)
    {
        ResidueStdOp t(m, s);
        m.mul(*this, *this, t);
        return *this;
    }
    friend ResidueStdOp operator*(uint64_t const s, ResidueStdOp const & t)
    {
        ResidueStdOp r(t.m, s);
        t.m.mul(r, r, t);
        return r;
    }
    friend ResidueStdOp operator*(Integer const & s, ResidueStdOp const & t)
    {
        ResidueStdOp r(t.m, s);
        t.m.mul(r, r, t);
        return r;
    }
    friend Integer & operator*=(Integer & s, ResidueStdOp const & t)
    {
        ResidueStdOp r(t.m, s);
        t.m.mul(r, r, t);
        s = r;
        return s;
    }

    ResidueStdOp operator/(ResidueStdOp const & s) const
    {
        assertValid(s);
        ResidueStdOp r(m);
        m.inv(r, s);
        m.mul(r, *this, r);
        return r;
    }
    ResidueStdOp & operator/=(ResidueStdOp const & s)
    {
        assertValid(s);
        ResidueStdOp i(m);
        m.inv(i, s);
        m.mul(*this, *this, i);
        return *this;
    }
    ResidueStdOp operator/(uint64_t const s) const
    {
        ResidueStdOp r(m);
        if (CONSTANT_P(s) && s == 2) {
            m.div2(r, *this);
        } else if (CONSTANT_P(s) && s == 3) {
            m.div3(r, *this);
        } else if (CONSTANT_P(s) && s == 5) {
            m.div5(r, *this);
        } else if (CONSTANT_P(s) && s == 7) {
            m.div7(r, *this);
        } else if (CONSTANT_P(s) && s == 11) {
            m.div11(r, *this);
        } else if (CONSTANT_P(s) && s == 13) {
            m.div13(r, *this);
        } else {
            ResidueStdOp r(m);
            m.set(r, s);
            m.inv(r, r);
            m.mul(r, *this, r);
        }
        return r;
    }
    ResidueStdOp & operator/=(uint64_t const s)
    {
        return *this /= ResidueStdOp(m, s);
    }
    friend ResidueStdOp operator/(uint64_t const s, ResidueStdOp const & t)
    {
        if (CONSTANT_P(s) && s == 1) { /* For 1/t, don't multiply by 1 */
            ResidueStdOp r(t.m);
            t.m.inv(r, t);
            return r;
        } else {
            return ResidueStdOp(t.m, s) / t;
        }
    }
    ResidueStdOp operator/(Integer const & s) const
    {
        return *this / ResidueStdOp(m, s);
    }
    ResidueStdOp & operator/=(Integer const & s)
    {
        return *this /= ResidueStdOp(m, s);
    }
    friend ResidueStdOp operator/(Integer const & s, ResidueStdOp const & t)
    {
        return ResidueStdOp(t.m, s) / t;
    }
    friend Integer & operator/=(Integer & s, ResidueStdOp const & t)
    {
        return s = (Integer)(ResidueStdOp(t.m, s) / t);
    }

    ResidueStdOp pow(uint64_t const e)
    {
        ResidueStdOp r(m);
        m.pow_u64(r, *this, e);
        return r;
    }
    ResidueStdOp pow(Integer const & e)
    {
        ResidueStdOp r(m);
        m.pow(r, *this, e);
        return r;
    }
    ResidueStdOp pow(uint64_t const * e, size_t const l)
    {
        ResidueStdOp r(m);
        m.pow(r, *this, e, l);
        return r;
    }
    ResidueStdOp chebyshevV(uint64_t const e)
    {
        ResidueStdOp r(m);
        m.V(r, *this, e);
        return r;
    }
    ResidueStdOp chebyshevV(Integer const & e)
    {
        ResidueStdOp r(m);
        m.V(r, *this, e);
        return r;
    }
    ResidueStdOp chebyshevV(uint64_t const * e, size_t const l)
    {
        ResidueStdOp r(m);
        m.V(r, *this, e, l);
        return r;
    }

    /* Should we call it index or gcd? I can't decide */
    Integer index() const
    {
        Integer g;
        m.gcd(g, *this);
        return g;
    }
    Integer gcd() const
    {
        Integer g;
        m.gcd(g, *this);
        return g;
    }
    int jacobi() const { return m.jacobi(*this); }
};

} /* namespace arithxx_details */

#endif	/* UTILS_ARITHXX_RESIDUE_STD_OP_HPP_ */
