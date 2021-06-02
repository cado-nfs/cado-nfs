#include "cado.h" // IWYU pragma: keep
#include "tests_common.h"
#include <cstdint>
#include <cinttypes>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>                                // for strlen
#include <gmp.h>                                   // for mpz_ptr, gmp_sscanf
#include <iostream>                                // for operator<<, endl
#include "cxx_mpz.hpp"                             // for cxx_mpz, operator==
#include "modint.hpp"                              // for operator<<
#include "mod64.hpp"
#include "modredc64.hpp"
#include "modredc126.hpp"
#include "mod_mpz_new.hpp"
#include "sieve/ecm/ec_arith_Weierstrass_new.hpp"
#include "macros.h"

template <typename MODULUS>
class TestWeierstrass {
    int verbose;
public:
    typedef MODULUS Modulus;
    typedef typename MODULUS::Residue Residue;
    typedef typename MODULUS::Integer Integer;
    typedef ECWeierstrass<MODULUS> Curve;
    typedef typename ECWeierstrass<MODULUS>::AffinePoint AffinePoint;
    typedef typename ECWeierstrass<MODULUS>::ProjectivePoint ProjectivePoint;
    
    TestWeierstrass(int verbose) : verbose(verbose) {}
    
    Modulus *initMod(const cxx_mpz &s) const {
        Integer i;
        if (!s.fits<Integer>())
            return NULL;
        i = s;
        if (!Modulus::valid(i))
            return NULL;
        Modulus *m = new Modulus(i);
        return m;
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
        const int operation, const Point &p1, const Point &p2,
        const uint64_t mult, const Point &pReference,
        const uint64_t expectedOrder ) const;

    bool parseLine(const char *line) const;

};

template<typename MODULUS>
template <typename Point>
bool
TestWeierstrass<MODULUS>::oneTest ( const Curve &c,
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
            std::cerr << "Computed order is " << order << " but expected " << expectedOrder << std::endl;
            return false;
        }
        return true;
    } else {
        abort();
    }

    if ( !pResult.isValid() || pReference != pResult ) {
        ok = false;
        std::cerr << "Computed point is wrong" << std::endl;
        std::cerr << "Testing " << typeid(*this).name() << "::" << operationNames[operation] << ":" << std::endl;
        std::cerr << c << std::endl;
        std::cerr << "p1 = " << p1 << std::endl;
        if (operation == 0)
            std::cerr << "p2 = " << p2 << std::endl;
        if (operation == 2)
            std::cerr << "mult = " << mult << std::endl;
        std::cerr << "pReference = " << pReference << std::endl;
        std::cerr << "pResult = " << pResult << std::endl;
    }
    return ok;
}

template<typename MODULUS>
bool
TestWeierstrass<MODULUS>::parseLine(const char *line) const {
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
        std::cerr << "Did not understand line:" << std::endl << line << std::endl;
        return false;
    }

    Modulus *m = initMod(M);
    if (m == NULL) {
        if (verbose) {
            std::cout << "Could not process modulus " << M << std::endl;
        }
        return true; /* This Modulus type can't test this value. This is not an error. */
    }
    Residue a ( *m );
    bool ok = true;
    if (!setResidue(a, *m, A)) {
        std::cerr << "Could not set point" << std::endl;
        ok = false;
    }
    Curve c ( *m, a );

    if (affine) {
        AffinePoint p1 ( c ), p2 ( c ), pReference ( c );

        if (!setPoint ( p1, *m, P1x, P1y ) ||
            !setPoint ( p2, *m, P2x, P2y ) ||
            !setPoint ( pReference, *m, PReferencex, PReferencey )) {
            std::cerr << "Could not set point" << std::endl;
            ok = false;
        } else {
            ok &= oneTest(c, operation, p1, p2, mult, pReference, order);
        }
    } else {
        ProjectivePoint p1 ( c ), p2 ( c ), pReference ( c );

        if (!setPoint ( p1, *m, P1x, P1y, P1z ) ||
            !setPoint ( p2, *m, P2x, P2y, P2z ) ||
            !setPoint ( pReference, *m, PReferencex, PReferencey, PReferencez )) {
            std::cerr << "Could not set point" << std::endl;
            ok = false;
        } else {
            ok &= oneTest(c, operation, p1, p2, mult, pReference, order);
        }
    }
    delete m;
    return ok;
}

// coverity[root_function]
int main(int argc, const char *argv[]) {
    bool ok = true;
    tests_common_cmdline (&argc, &argv, PARSE_VERBOSE);
    int verbose = tests_common_get_verbose ();

    if (argc < 2) {
        fprintf(stderr, "Input file missing\n");
        exit(EXIT_FAILURE);
    }
    FILE *inputfile = fopen(argv[1], "r");
    if (inputfile == NULL) {
        perror("Could not open file");
        exit(EXIT_FAILURE);
    }
    
    TestWeierstrass<Modulus64> test1(verbose);
    TestWeierstrass<ModulusREDC64> test2(verbose);
    TestWeierstrass<ModulusREDC126> test3(verbose);
    TestWeierstrass<ModulusMPZ> test4(verbose);

    constexpr size_t buflen = 1024;
    char line[buflen];
    while (!feof(inputfile)) {
        if (fgets(line, sizeof(line), inputfile) == NULL)
            break;
        size_t len = strlen(line);
        if (len > 0) {
            if (line[len - 1] == '\n')
                len--;
            else if (len == buflen - 1) {
                std::cerr << "No newline in input. Buffer too small?" << std::endl;
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
