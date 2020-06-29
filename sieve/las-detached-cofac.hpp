#ifndef LAS_DETACHED_COFAC_HPP_
#define LAS_DETACHED_COFAC_HPP_

#include <memory>                    // for shared_ptr
#include "las-cofac-standalone.hpp"  // for cofac_standalone
#include "threadpool.hpp"            // for task_parameters, task_result
class nfs_aux;
class nfs_work_cofac;
struct relation;

struct detached_cofac_parameters : public cofac_standalone, public task_parameters {
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    detached_cofac_parameters(std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, cofac_standalone&& C) : cofac_standalone(C), wc_p(wc_p), aux_p(aux_p) {}
};

struct detached_cofac_result : public task_result {
    std::shared_ptr<relation> rel_p;
};

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int);

#endif	/* LAS_DETACHED_COFAC_HPP_ */
