#include "cado.h"
#include "tests_common.h"
#include <cstdint>
#include <cinttypes>
#include "mod64.hpp"
#include "modredc64.hpp"
#include "modredc126.hpp"
#include "mod_mpz_new.hpp"
#include "sieve/ecm/ec_arith_Weierstrass_new.hpp"

template <typename MODULUS>
class TestWeierstrassAffine {
    int verbose;
public:
    typedef MODULUS Modulus;
    typedef typename MODULUS::Residue Residue;
    typedef typename MODULUS::Integer Integer;
    typedef ECWeierstrass<MODULUS> Curve;
    typedef typename ECWeierstrass<MODULUS>::AffinePoint Point;
    
    TestWeierstrassAffine(int verbose) : verbose(verbose) {}
    
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
    setPoint ( Point &p, const Modulus &m, const cxx_mpz &_x, const cxx_mpz &_y ) const
    {
        Residue x ( m ), y ( m );

        if (!setResidue ( x, m, _x ) || !setResidue ( y, m, _y ))
            return false;
        p.set ( x, y );
        return true;
    }

    bool oneTest ( const cxx_mpz &M, const cxx_mpz &A,
        const int operation, const cxx_mpz &x1, const cxx_mpz &y1, const cxx_mpz &x2,
        const cxx_mpz &y2, const uint64_t mult, const cxx_mpz &xReference, const cxx_mpz &yReference,
        const bool expectedIsFinite, const uint64_t expectedOrderexpectedIsFinite ) const;

    bool parseLine(const char *line) const;

};

template<typename MODULUS>
bool
TestWeierstrassAffine<MODULUS>::oneTest ( const cxx_mpz &M, const cxx_mpz &A,
        const int operation, const cxx_mpz &x1, const cxx_mpz &y1, const cxx_mpz &x2,
        const cxx_mpz &y2, const uint64_t mult, const cxx_mpz &xReference, const cxx_mpz &yReference,
        const bool expectedIsFinite, const uint64_t expectedOrder ) const
{
    const char *operationNames[4] = {"add", "dbl", "mul", "ord"};
    bool ok = true;
    Integer order;
    Modulus *m = initMod(M);
    if (m == NULL) {
        if (verbose) {
            std::cout << "Could not process modulus " << M << std::endl;
        }
        return true; /* This Modulus type can't test this value. This is not an error. */
    }
    Residue a ( *m ), x ( *m ), y ( *m );
    setResidue(a, *m, A);
    Curve c ( *m, a );
    Point p1 ( c ), p2 ( c ), pReference ( c ), pResult ( c );

    setPoint ( p1, *m, x1, y1 );
    setPoint ( p2, *m, x2, y2 );
    setPoint ( pReference, *m, xReference, yReference );

    if (operation < 0 || operation > 3)
        abort();

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
    
    if ( (!pResult.is0()) != expectedIsFinite ) {
        std::cerr << "Expected point " << ( expectedIsFinite ? "not " : "" ) << "at infinity" << std::endl;
        std::cerr << "but result is " << (pResult.is0() ? "" : "not ") << "point at infinity" << std::endl;
        ok = false;
    } else if ( (!pResult.is0()) && pReference != pResult ) {
        std::cerr << "Computed point is wrong" << std::endl;
        ok = false;
    }
    if ( !ok ) {
        std::cerr << "Testing TestWeierstrassAffine<MODULUS>::" << operationNames[operation] << ":" << std::endl;
        std::cerr << c << std::endl;
        std::cerr << "p1 = " << p1 << std::endl;
        if (operation == 0)
            std::cerr << "p2 = " << p2 << std::endl;
        if (operation == 3)
            std::cerr << "mult = " << mult << std::endl;
        std::cerr << "pReference = " << pReference << std::endl;
        std::cerr << "pResult = " << pResult << std::endl;
    }
    delete m;
    return ok;
}

template<typename MODULUS>
bool
TestWeierstrassAffine<MODULUS>::parseLine(const char *line) const {
    cxx_mpz M, A, P1x, P1y, P2x, P2y, PReferencex, PReferencey;
    uint64_t mult = 0, order = 0;
    int operation;

    if (gmp_sscanf(line, "add M %Zd A %Zd P1 %Zd:%Zd P2 %Zd:%Zd RESULT %Zd:%Zd",
                   (mpz_ptr) M, (mpz_ptr) A, (mpz_ptr) P1x, (mpz_ptr) P1y, (mpz_ptr) P2x, (mpz_ptr) P2y, (mpz_ptr) PReferencex, (mpz_ptr) PReferencey) == 8) {
        operation = 0;
    } else if (gmp_sscanf(line, "dbl M %Zd A %Zd P1 %Zd:%Zd RESULT %Zd:%Zd",
                          &M, &A, &P1x, &P1y, &PReferencex, &PReferencey) == 6) {
        operation = 1;
    } else if (gmp_sscanf(line, "mul M %Zd A %Zd P1 %Zd:%Zd MULT %" SCNu64 " RESULT %Zd:%Zd",
                          &M, &A, &P1x, &P1y, &mult, &PReferencex, &PReferencey) == 7) {
        operation = 2;
    } else if (gmp_sscanf(line, "ord M %Zd A %Zd P1 %Zd:%Zd ORDER %" SCNu64,
                          &M, &A, &P1x, &P1y, &order) == 5) {
        operation = 3;
    } else {
        std::cerr << "Did not understand line:" << std::endl << line << std::endl;
        return false;
    }
    const bool expectedIsFinite = (PReferencex != 0 || PReferencey != 0);
    return oneTest(M, A, operation, P1x, P1y, P2x, P2y, mult, PReferencex, PReferencey, expectedIsFinite, order);
}


int main(int argc, const char *argv[]) {
    bool ok = true;
    int verbose = 0;
    tests_common_cmdline (&argc, &argv, PARSE_VERBOSE);
    verbose = tests_common_get_verbose ();

    if (argc < 2) {
        fprintf(stderr, "Input file missing\n");
        exit(EXIT_FAILURE);
    }
    FILE *inputfile = fopen(argv[1], "r");
    if (inputfile == NULL) {
        perror("Could not open file");
        exit(EXIT_FAILURE);
    }
    
    TestWeierstrassAffine<Modulus64> test1(verbose);
    TestWeierstrassAffine<ModulusREDC64> test2(verbose);
    TestWeierstrassAffine<ModulusREDC126> test3(verbose);
    TestWeierstrassAffine<ModulusMPZ> test4(verbose);

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
