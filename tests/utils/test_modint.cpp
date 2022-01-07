#include "cado.h" // IWYU pragma: keep
#include <cstdlib>
#include <sstream>
#include <iostream> // cerr
#include <string>
#include "cxx_mpz.hpp"
#include "modint.hpp"
#include "tests_common.h"

template <typename INTEGER>
class Tests {
    typedef INTEGER Integer;

    bool test_operator_assign_cxx_mpz() {
        Integer i;
        cxx_mpz z;
        
        i = 1;
        i = z;
        if (i != 0) {
            std::cerr << "Setting to 0 failed" << std::endl;
            return false;
        }
        z = 5;
        i = z;
        if (i != 5) {
            std::cerr << "Setting to 5 failed" << std::endl;
            return false;
        }
        return true;
    }
    
    bool test_stream_operator() {
        cxx_mpz V;
        Integer v;
        bool ok = true;
        for (int i = 0; i < 256; i++) {
            V = 1;
            V <<= i;
            if (!V.fits<Integer>())
                break;
            v = V;
            std::stringstream V_str, v_str;
            V_str << V;
            v_str << v;
            if (V_str.str() != v_str.str()) {
                std::cerr << "GMP output: " << V << ", Integer output: " << v << std::endl;
                ok = false;;
            }
        }
        return ok;
    }
public:
    bool runTests() {
        bool ok = true;
        ok &= test_operator_assign_cxx_mpz();
        ok &= test_stream_operator();
        return ok;
    }
};

int main(int argc, const char **argv) {
    bool ok = true;
  
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);

    try {
        Tests<Integer64> test1;
        ok &= test1.runTests();
        
        Tests<Integer128> test2;
        ok &= test2.runTests();
    } catch (std::runtime_error const & e) {
        std::cerr << e.what();
        ok = false;
    }
      
    tests_common_clear();
    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
