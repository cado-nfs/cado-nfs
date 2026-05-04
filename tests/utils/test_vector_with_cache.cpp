#include "cado.h"       // IWYU pragma: keep

#include <cstddef>

#include "fmt/base.h"

#include "vector_with_cache.hpp"
#include "cxx_mpz.hpp"

int main()
{
    {
        vector_with_cache<cxx_mpz, 4> A;

        for(size_t i = 0 ; i < 20 ; i++) {
            A.emplace_back(i);
            ASSERT_ALWAYS(A.last() == i);
        }

        for(size_t i = 0 ; i < A.size() ; i++) {
            ASSERT_ALWAYS(A[i] == i);
            A[i] *= 2;
        }

        for(auto const & v : A)
            fmt::print("{}\n", v);
    }
}

