#include "cado.h"
#include <gmp.h>
#include <string.h>
#include "minimum_spanning_tree.hpp"
#include <iostream>
#include <limits>
#include <climits>
#include <string>
#include <sstream>
#include "macros.h"
#include "utils.h"

using namespace std;


string do_random_test(gmp_randstate_t rstate, int maxsize)
{
    int n;
    for( ; (n = gmp_urandomm_ui(rstate, maxsize)) <= 2 ; );
    matrix<int> A(n, n, INT_MAX);
    int m = 2 * n + gmp_urandomm_ui(rstate, n * n - 2 * n);
    for(int k = 0 ; k < m ; k++) {
        int i = gmp_urandomm_ui(rstate, n);
        int j = gmp_urandomm_ui(rstate, n);
        int w = gmp_urandomm_ui(rstate, 1000);
        A(i,j) = w;
    }
    auto mst = minimum_spanning_tree(A);

    stringstream os;
    os << n << " " << m << " ";
    if (mst.first == numeric_limits<decltype(mst.first)>::max())
        os << -1;
    else
        os << mst.first;
    return os.str();
}


int wikipedia_example()
{
    // this is the graph from
    // https://en.wikipedia.org/wiki/File:Minimum_spanning_tree.svg
    matrix<int> A(10,10,INT_MAX);
    A(0,1)=4;
     A(0,2)=1;
      A(0,3)=4;
     A(1,2)=5;
       A(1,4)=9;
        A(1,5)=9;
          A(1,7)=7;
      A(2,3)=3;
          A(2,7)=9;
          A(3,7)=10;
            A(3,9)=18;
        A(4,5)=2;
         A(4,6)=4;
           A(4,8)=6;
         A(5,6)=2;
          A(5,7)=8;
          A(6,7)=9;
           A(6,8)=3;
            A(6,9)=9;
            A(7,9)=8;
            A(8,9)=9;
    for(int i = 0 ; i < 10 ; i++) {
        for(int j = i + 1 ; j < 10 ; j++) {
            A(j,i)=A(i,j);
        }
    }
    auto mst = minimum_spanning_tree(A);
    cout << "MST has weight " << mst.first << "\n";
    for(auto const & x : mst.second) {
        cout << " " << x.first << "--" << x.second << "\n";
    }
    ASSERT_ALWAYS(mst.first == 38);

    /*
     * we expect an MST of weight 38, which the code would be welcome to
     * print nicely as follows
     0--2--3
     |
     \--1--7--9
           |
           \--5--6--8
              |
              \--4
     */

    return 0;
}


int main(int argc, char * argv[])
{
    int seed = 0;
    int ntrials = 10;
    int timings = 0;
    int maxsize = 100;
    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        if (strcmp(argv[0], "--seed") == 0) {
            seed = atoi(argv[1]);
            argc--,argv++;
            continue;
        }
        if (strcmp(argv[0], "--maxsize") == 0) {
            maxsize = atoi(argv[1]);
            argc--,argv++;
            continue;
        }
        if (strcmp(argv[0], "--ntrials") == 0) {
            ntrials = atoi(argv[1]);
            argc--,argv++;
            continue;
        }
        if (strcmp(argv[0], "--timings") == 0) {
            timings=1;
            continue;
        }
    }

    wikipedia_example();

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);

    for(int i = 0 ; i < ntrials ; i++) {
        if (timings) {
            double tt = seconds();
            int nrun = 0;
            string s;
            for( ; (seconds()-tt) < 1 ; nrun++) {
                s = do_random_test(rstate, maxsize);
            }
            printf("%s # %g\n", s.c_str(), (seconds()-tt) / nrun);
        } else {
            cout << do_random_test(rstate, maxsize) << "\n";
        }
    }

    gmp_randclear(rstate);
}

