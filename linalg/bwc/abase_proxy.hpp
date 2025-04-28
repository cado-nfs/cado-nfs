#ifndef CADO_LINALG_BWC_ABASE_PROXY_HPP
#define CADO_LINALG_BWC_ABASE_PROXY_HPP

#include <map>
#include <memory>

#include <gmp.h>

#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "bw-common.h"
#include "parallelizing_info.hpp"

struct abase_proxy {
    parallelizing_info_ptr pi;
    std::unique_ptr<arith_generic> A;
    pi_datatype_ptr A_pi;

    abase_proxy(parallelizing_info_ptr pi, unsigned int width)
        : pi(pi)
        , A(arith_generic::instance(bw->p, width))
        , A_pi(pi_alloc_arith_datatype(pi, A.get()))
    {
    }
    abase_proxy(abase_proxy &) = delete;
    abase_proxy & operator=(abase_proxy const &) = delete;
    abase_proxy(abase_proxy &&) = default;
    abase_proxy & operator=(abase_proxy &&) = default;
    static abase_proxy most_natural(parallelizing_info_ptr pi)
    {
        const unsigned int width = mpz_cmp_ui(bw->p, 2) == 0 ? 64 : 1;
        return { pi, width };
    }
    std::map<arith_generic *, std::shared_ptr<arith_cross_generic>> tdict;
    arith_cross_generic * templates(arith_generic * A1)
    {
        auto it = tdict.find(A1);
        if (it == tdict.end())
            tdict[A1] = std::shared_ptr<arith_cross_generic>(
                arith_cross_generic::instance(A.get(), A1));
        return tdict[A1].get();
    }
    ~abase_proxy() { pi_free_arith_datatype(pi, A_pi); }
};

#endif /* LINALG_BWC_ABASE_PROXY_HPP_ */
