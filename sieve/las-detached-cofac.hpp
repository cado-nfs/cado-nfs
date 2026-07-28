#ifndef CADO_LAS_DETACHED_COFAC_HPP
#define CADO_LAS_DETACHED_COFAC_HPP

#include "las-cofac-standalone.hpp"

class nfs_aux;
class nfs_work_cofac;
class worker_thread;
struct relation;

extern relation detached_cofac(worker_thread * worker,
        nfs_work_cofac & wc,
        nfs_aux & aux,
        cofac_standalone && cur);

#endif	/* CADO_LAS_DETACHED_COFAC_HPP */
