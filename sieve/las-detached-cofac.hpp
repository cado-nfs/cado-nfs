#ifndef CADO_LAS_DETACHED_COFAC_HPP
#define CADO_LAS_DETACHED_COFAC_HPP

#include <memory>
#include <utility>

#include "las-cofac-standalone.hpp"
#include "threadpool.hpp"

class nfs_aux;
class nfs_work_cofac;
struct relation;

struct detached_cofac_parameters
    : public cofac_standalone
    , public task_parameters
{
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    detached_cofac_parameters(
            std::shared_ptr<nfs_work_cofac> wc_p,
            std::shared_ptr<nfs_aux> aux_p,
            cofac_standalone&& C)
    : cofac_standalone(C)
    , wc_p(std::move(wc_p))
    , aux_p(std::move(aux_p))
    { }
};

struct detached_cofac_result : public task_result {
    std::shared_ptr<relation> rel_p;
};

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int);

#endif	/* CADO_LAS_DETACHED_COFAC_HPP */
