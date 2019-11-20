#ifndef LINGEN_MUL_SUBSTEPS_BASE_HPP_
#define LINGEN_MUL_SUBSTEPS_BASE_HPP_

#include <array>
#include <string>

struct op_mul_or_mp_base {
    enum op_type_t { OP_MP, OP_MUL } op_type;
    static const char * op_name(op_type_t op_type) {
        switch(op_type) {
            case OP_MP: return "MP";
            case OP_MUL: return "MUL";
        }
        throw std::runtime_error("bad op");
    }
    const char * op_name() const { return op_name(op_type); }
    virtual std::array<size_t, 3> get_alloc_sizes() const = 0;
    virtual std::string explain() const = 0;
    virtual ~op_mul_or_mp_base() = default;
};


#endif	/* LINGEN_MUL_SUBSTEPS_BASE_HPP_ */
