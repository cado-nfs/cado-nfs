#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>                                // for strlen

#include <typeinfo>
#include <iostream>                                // for operator<<, endl
#include <memory>

#include <gmp.h>                                   // for mpz_ptr, gmp_sscanf

#include "tests_common.h"
#include "cxx_mpz.hpp"                             // for cxx_mpz, operator==
#include "arithxx/mod64.hpp"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/mod_mpz_new.hpp"
#include "sieve/ecm/ec_arith_Weierstrass_new.hpp"
#include "macros.h"

template <typename layer>
class TestWeierstrass {
    int verbose;
public:
    using Modulus = typename layer::Modulus;
    using Residue = typename layer::Residue;
    using Integer = typename layer::Integer;
    using Curve = ECWeierstrass<layer>;
    using AffinePoint = typename ECWeierstrass<layer>::AffinePoint;
    using ProjectivePoint = typename ECWeierstrass<layer>::ProjectivePoint;
    
    explicit TestWeierstrass(int verbose) : verbose(verbose) {}
    
    std::unique_ptr<Modulus> initMod(const cxx_mpz &s) const {
        Integer i;
        if (!s.fits<Integer>())
            return {};
        i = s;
        if (!Modulus::valid(i))
            return {};
        return std::unique_ptr<Modulus>(new Modulus(i));
    }

    bool
    setResidue ( Residue &r, const Modulus &m, const cxx_mpz &v ) const
    {
        Integer i;
        if (!v.fits<Integer> ())
            return false;
        i = v;
        m.set ( r, i );
        return true;
    }

    bool
    setPoint ( AffinePoint &p, const Modulus &m, const cxx_mpz &_x, const cxx_mpz &_y ) const
    {
        Residue x ( m ), y ( m );
        if (_x == 0 && _y == 0) {
            p.set0();
            return true;
        }
        if (!setResidue ( x, m, _x ) || !setResidue ( y, m, _y ))
            return false;
        p.set ( x, y );
        return true;
    }

    bool
    setPoint ( ProjectivePoint &p, const Modulus &m, const cxx_mpz &_x, const cxx_mpz &_y, const cxx_mpz &_z ) const
    {
        Residue x ( m ), y ( m ), z ( m );

        if (!setResidue ( x, m, _x ) || !setResidue ( y, m, _y ) || !setResidue ( z, m, _z ) )
            return false;
        p.set ( x, y, z );
        return true;
    }

    template <typename Point>
    bool oneTest ( const Curve &c,
        int operation, const Point &p1, const Point &p2,
        uint64_t mult, const Point &pReference,
        uint64_t expectedOrder ) const;

    bool parseLine(const char *line) const;

};

template<typename layer>
template <typename Point>
bool
TestWeierstrass<layer>::oneTest ( const Curve &c,
        const int operation,  const Point &p1, const Point &p2,
        const uint64_t mult, const Point &pReference,
        const uint64_t expectedOrder ) const
{
    const char *operationNames[4] = {"add", "dbl", "mul", "ord"};
    bool ok = true;
    Integer order;
    Point pResult ( c );

    ASSERT_ALWAYS(0 <= operation && operation < 4);

    if (operation == 0) {
        pResult = p1 + p2;
    } else if (operation == 1) {
        p1.dbl ( pResult );
    } else if (operation == 2) {
        pResult = p1 * mult;
    } else if (operation == 3) {
        order = p1.point_order (1, 0, 1);
        if (order != expectedOrder) {
            std::cerr << "Computed order is " << order << " but expected " << expectedOrder << "\n";
            return false;
        }
        return true;
    } else {
        abort();
    }

    if ( !pResult.isValid() || pReference != pResult ) {
        ok = false;
        std::cerr << "Computed point is wrong\n";
        std::cerr << "Testing " << typeid(*this).name() << "::" << operationNames[operation] << ":\n";
        std::cerr << c << "\n";
        std::cerr << "p1 = " << p1 << "\n";
        if (operation == 0)
            std::cerr << "p2 = " << p2 << "\n";
        if (operation == 2)
            std::cerr << "mult = " << mult << "\n";
        std::cerr << "pReference = " << pReference << "\n";
        std::cerr << "pResult = " << pResult << "\n";
    }
    return ok;
}

template<typename layer>
bool
TestWeierstrass<layer>::parseLine(const char *line) const {
    cxx_mpz M, A, P1x, P1y, P1z, P2x, P2y, P2z, PReferencex, PReferencey, PReferencez;
    uint64_t mult = 0, order = 0;
    int operation;
    bool affine;

    if (gmp_sscanf(line, "add M %Zd A %Zd P1 %Zd:%Zd P2 %Zd:%Zd RESULT %Zd:%Zd",
                   (mpz_ptr) M, (mpz_ptr) A, (mpz_ptr) P1x, (mpz_ptr) P1y, 
                   (mpz_ptr) P2x, (mpz_ptr) P2y, (mpz_ptr) PReferencex, (mpz_ptr) PReferencey) == 8) {
        operation = 0;
        affine = true;
    } else if (gmp_sscanf(line, "add M %Zd A %Zd P1 %Zd:%Zd:%Zd P2 %Zd:%Zd:%Zd RESULT %Zd:%Zd:%Zd",
                   (mpz_ptr) M, (mpz_ptr) A, (mpz_ptr) P1x, (mpz_ptr) P1y, (mpz_ptr) P1z,
                   (mpz_ptr) P2x, (mpz_ptr) P2y, (mpz_ptr) P2z, 
                   (mpz_ptr) PReferencex, (mpz_ptr) PReferencey, (mpz_ptr) PReferencez) == 11) {
        operation = 0;
        affine = false;
    } else if (gmp_sscanf(line, "dbl M %Zd A %Zd P1 %Zd:%Zd RESULT %Zd:%Zd",
                          &M, &A, &P1x, &P1y, &PReferencex, &PReferencey) == 6) {
        operation = 1;
        affine = true;
    } else if (gmp_sscanf(line, "dbl M %Zd A %Zd P1 %Zd:%Zd:%Zd RESULT %Zd:%Zd:%Zd",
                          &M, &A, &P1x, &P1y, &P1z, &PReferencex, &PReferencey, &PReferencez) == 8) {
        operation = 1;
        affine = false;
    } else if (gmp_sscanf(line, "mul M %Zd A %Zd P1 %Zd:%Zd MULT %" SCNu64 " RESULT %Zd:%Zd",
                          &M, &A, &P1x, &P1y, &mult, &PReferencex, &PReferencey) == 7) {
        operation = 2;
        affine = true;
    } else if (gmp_sscanf(line, "mul M %Zd A %Zd P1 %Zd:%Zd:%Zd MULT %" SCNu64 " RESULT %Zd:%Zd:%Zd",
                          &M, &A, &P1x, &P1y, &P1z, &mult, &PReferencex, &PReferencey, &PReferencez) == 9) {
        operation = 2;
        affine = false;
    } else if (gmp_sscanf(line, "ord M %Zd A %Zd P1 %Zd:%Zd ORDER %" SCNu64,
                          &M, &A, &P1x, &P1y, &order) == 5) {
        operation = 3;
        affine = true;
    } else {
        std::cerr << "Did not understand line:\n" << line << "\n";
        return false;
    }

    auto pm = initMod(M);
    if (!pm) {
        if (verbose) {
            std::cout << "Could not process modulus " << M << "\n";
        }
        return true; /* This Modulus type can't test this value. This is not an error. */
    }

    Modulus & m = *pm;

    Residue a ( m );
    bool ok = true;
    if (!setResidue(a, m, A)) {
        std::cerr << "Could not set point\n";
        ok = false;
    }
    Curve const c ( m, a );

    if (affine) {
        AffinePoint p1 ( c ), p2 ( c ), pReference ( c );

        if (!setPoint ( p1, m, P1x, P1y ) ||
            !setPoint ( p2, m, P2x, P2y ) ||
            !setPoint ( pReference, m, PReferencex, PReferencey )) {
            std::cerr << "Could not set point\n";
            ok = false;
        } else {
            ok &= oneTest(c, operation, p1, p2, mult, pReference, order);
        }
    } else {
        ProjectivePoint p1 ( c ), p2 ( c ), pReference ( c );

        if (!setPoint ( p1, m, P1x, P1y, P1z ) ||
            !setPoint ( p2, m, P2x, P2y, P2z ) ||
            !setPoint ( pReference, m, PReferencex, PReferencey, PReferencez )) {
            std::cerr << "Could not set point\n";
            ok = false;
        } else {
            ok &= oneTest(c, operation, p1, p2, mult, pReference, order);
        }
    }
    return ok;
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    bool ok = true;
    tests_common_cmdline (&argc, &argv, PARSE_VERBOSE);
    int const verbose = tests_common_get_verbose ();

    if (argc < 2) {
        fprintf(stderr, "Input file missing\n");
        exit(EXIT_FAILURE);
    }
    FILE *inputfile = fopen(argv[1], "r");
    DIE_ERRNO_DIAG(!inputfile, "fopen(%s)", argv[1]);
    
    TestWeierstrass<arithxx_mod64> const test1(verbose);
    TestWeierstrass<arithxx_modredc64> const test2(verbose);
    TestWeierstrass<arithxx_modredc126> const test3(verbose);
    TestWeierstrass<arithxx_mod_mpz_new> const test4(verbose);

    constexpr size_t buflen = 1024;
    char line[buflen];
    while (!feof(inputfile)) {
        if (fgets(line, sizeof(line), inputfile) == nullptr)
            break;
        size_t len = strlen(line);
        if (len > 0) {
            if (line[len - 1] == '\n')
                len--;
            else if (len == buflen - 1) {
                std::cerr << "No newline in input. Buffer too small?\n";
                abort();
            }
        }
        if (len == 0)
            continue;
        ok &= test1.parseLine(line);
        ok &= test2.parseLine(line);
        ok &= test3.parseLine(line);
        ok &= test4.parseLine(line);
    }
    
    fclose(inputfile);
    tests_common_clear ();
    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
