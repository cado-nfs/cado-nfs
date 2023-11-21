#include "cado.h"

#include "facul_method.hpp"
#include "facul_ecm.h"  // for ecm_plan_t, ecm_make_plan, ecm_clear_plan
#include "macros.h"     // for ASSERT_ALWAYS, ASSERT, MAYBE_UNUSED
#include "pm1.h"        // for pm1_plan_t, pm1_clear_plan, pm1_make_plan
#include "pp1.h"        // for pp1_plan_t, pp1_clear_plan, pp1_make_plan

facul_method::~facul_method()
{
    if (method == PM1_METHOD)
        pm1_clear_plan ((pm1_plan_t*) plan);
    else if (method == PP1_27_METHOD)
        pp1_clear_plan ((pp1_plan_t*) plan);
    else if (method == PP1_65_METHOD)
        pp1_clear_plan ((pp1_plan_t*) plan);
    else if (method == EC_METHOD)
        ecm_clear_plan ((ecm_plan_t*) plan);
    free(plan);
}


facul_method::facul_method(parameters const & p, const int verbose)
{
    method = p.method;
    plan = NULL;

    /* we don't need to store B1 and B2 in the method object. They all go
     * in the bytecode anyway
     */
    int B1 = p.B1;
    int B2 = p.B2;

    switch(method) {
        case PP1_27_METHOD:
        case PP1_65_METHOD:
            plan = malloc (sizeof (pp1_plan_t));
            pp1_make_plan ((pp1_plan_t *) plan, B1, B2, verbose);
            break;
        case PM1_METHOD:
            plan = malloc (sizeof (pm1_plan_t));
            pm1_make_plan ((pm1_plan_t*) plan, B1, B2, verbose);
            break;
        case EC_METHOD:
            plan = malloc (sizeof (ecm_plan_t));
            if (!ec_parameter_is_valid (p.parameterization, p.parameter)) {
                fprintf (stderr,
                        "Parameter %lu is not valid with parametrization %d\n",
                        p. parameter, p.parameterization);
                exit (EXIT_FAILURE);
            }

            ecm_make_plan ((ecm_plan_t *) plan, B1, B2, p.parameterization, p.parameter, p.extra_primes, verbose);
            break;
    }

    ASSERT_ALWAYS(plan != NULL);
}

void facul_method_side::fix_is_last(std::vector<facul_method_side> & v)
{
    int is_last[2] = {1,1};
    for(auto it = v.rbegin() ; it != v.rend() ; ++it) {
        it->is_last = is_last[it->side];
        is_last[it->side] = 0;
    }
}
