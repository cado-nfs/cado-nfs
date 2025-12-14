#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <algorithm>
#include <string>

#include "fmt/base.h"
#include <gmp.h>

#include "point.hpp"
#include "tab_point.hpp"
#include "convex_hull.hpp"
#include "params.h"
#include "random_distributions.hpp"
#include "gmp_aux.h"
#include "macros.h"

struct coarse_cmp_points_coordinates {
    static int cmp(double a, double b)
    {
        const double diff = a - b;
        constexpr const double precision = 1e-10;

        if (diff <= -precision) return -1;
        if (diff >= precision)  return 1;
        return 0;
    }

    bool operator()(point const & a, point const & b) const
    {
        if (a.x < b.x) return true;
        if (b.x < a.x) return false;
        return a.y < b.y;
    }
};

struct coarse_cmp_table_points {
    bool operator()(tabular_point const & a, tabular_point const & b) const
    {
        return std::ranges::lexicographical_compare(
                a, b, coarse_cmp_points_coordinates());
    }
};

static bool operator<(tabular_point const & a, tabular_point const & b)
{
    return coarse_cmp_table_points()(a, b);
}
static bool operator==(tabular_point const & a, tabular_point const & b)
{
    return !(a < b) && !(b < a);
}

static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "seed", "seed");
    param_list_decl_usage(pl, "N", "number of points");
    param_list_decl_usage(pl, "x0", "lower bound");
    param_list_decl_usage(pl, "x1", "lower bound");
}


int main(int argc, char const * argv[])
{
    const char * argv0 = argv[0];
    cxx_param_list pl;
    cxx_gmp_randstate state;

    declare_usage(pl);
    argv++,argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);


    unsigned int N = 10;
    double x0 = -2;
    double x1 = 10;
    unsigned long seed = 0;

    param_list_parse(pl, "N", N);
    param_list_parse(pl, "x0", x0);
    param_list_parse(pl, "x1", x1);
    param_list_parse(pl, "seed", seed);

    ASSERT_ALWAYS(x0 < x1);
    gmp_randseed_ui(state, seed);

    /* a reference convex hull */
    tabular_point H;
    for(unsigned int i = 0 ; i < N ; i++) {
        const double x = x0 + (x1 - x0) * random_uniform(state);
        /* take any strictly convex function, really. Here we take the
         * square */
        H.emplace_back(point { .number=i, .x=x, .y=x * x });
    }
    std::ranges::sort(H, coarse_cmp_points_coordinates());

    ASSERT_ALWAYS(N >= 3);

    for(const double s : { 1, -1 }) {
        /* s == 1: take our convex hull, and add many interior points.
         * We should get the exact same convex hull
         *
         * s == -1: same, but now with _exterior_ points. We should _not_ get
         * the same convex hull
         */
        tabular_point C = H;

        for(unsigned int i = 0 ; i < 4 * N ; i++) {
            const int i0 = gmp_urandomm_ui(state, N);
            int i1; for(i1 = i0 ; i1 == i0 ; i1 = gmp_urandomm_ui(state, N));
            int i2; for(i2 = i0 ; i2 == i0 || i2 == i1 ; i2 = gmp_urandomm_ui(state, N));
            const point P0 = H[i0];
            const point P1 = H[i1];
            const point P2 = H[i2];
            double t0 = random_uniform(state);
            double t1 = random_uniform(state);
            double t2 = random_uniform(state);
            const double T = t0 + t1 + t2;
            t0 /= T;
            t1 /= T;
            t2 /= T;
            if (s == -1) {
                /* This still sums to one, but we now get an exterior
                 * point P' = P0 - (t1+t2)*(P12-P0) (where P12 is the
                 * weighted sum of P1 and P2)
                 */
                t0 = 2 - t0;
                t1 = - t1;
                t2 = - t2;
            }

            const point Pi { N + i,
                    t0 * P0.x + t1 * P1.x + t2 * P2.x,
                    t0 * P0.y + t1 * P1.y + t2 * P2.y };

            /*
            fmt::print("i={} s={} i0={} i1={} i2={} t0={} t1={} t2={}  -> ({},{})\n",
                    i, s, i0, i1, i2, t0, t1, t2, Pi.x, Pi.y);
                    */
            C.emplace_back(Pi);


        }

        const auto H2 = convex_hull(C);

        fmt::print("Test s=={}\n", s);
        for(auto const & P : C)
            fmt::print("({},{}),\n", P.x, P.y);
        fmt::print("Hull:\n");
        for(auto const & P : H2)
            fmt::print("({},{}),\n", P.x, P.y);

        if ((s == 1) != (H == H2)) {
            fmt::print(stderr, "Error in convex_hull\n");
            exit(EXIT_FAILURE);
        }
    }

    return EXIT_SUCCESS;
}
