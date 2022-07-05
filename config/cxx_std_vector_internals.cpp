#include "cxx_std_vector_ugly_accessor.hpp"

struct fred : public std::vector<int> {
    typedef std::vector<int> Base;
    fred(size_type n, int a) : Base(n, a) {}
    bool foo() {
        auto & u = cado_nfs::ugly_accessor(*this);
        void * x = u._start();
        void * y = u._finish();
        void * z = u._end_of_storage();
        return x == y || x == z;
    }
    Base::allocator_type & bar() {
        auto & u = cado_nfs::ugly_accessor(*this);
        return u._allocator();
    }
};

int main() {
    fred F(12, 1);
}
