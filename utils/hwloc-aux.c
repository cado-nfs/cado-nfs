#include "cado.h" // IWYU pragma: keep
#include <stddef.h>  // for NULL
#include <hwloc.h>   // for hwloc_type_sscanf_as_depth, HWLOC_API_VERSION
#include "macros.h"     // ASSERT_ALWAYS // IWYU pragma: keep
#include "hwloc-aux.h"

/* missing api calls in hwloc */
/* Part of this code is derived from hwloc-calc.h in hwloc-1.11.0,
 * which carries the following copyright notice.
 * {{{
 *   Copyright © 2004-2006 The Trustees of Indiana University and Indiana
 *   University Research and Technology Corporation.  All rights
 *   reserved.
 *
 *   Copyright © 2004-2005 The University of Tennessee and The University
 *   of Tennessee Research Foundation.  All rights reserved.
 *
 *   Copyright © 2004-2005 High Performance Computing Center Stuttgart,
 *   University of Stuttgart.  All rights reserved.
 *
 *   Copyright © 2004-2005 The Regents of the University of California.
 *   All rights reserved.
 *
 *   Copyright © 2009      CNRS
 *
 *   Copyright © 2009-2016 Inria.  All rights reserved.
 *
 *   Copyright © 2009-2015 Université Bordeaux
 *
 *   Copyright © 2009-2015 Cisco Systems, Inc.  All rights reserved.
 *
 *   Copyright © 2009-2012 Oracle and/or its affiliates.  All rights
 *   reserved.
 *
 *   Copyright © 2010      IBM
 *
 *   Copyright © 2010      Jirka Hladky
 *
 *   Copyright © 2012      Aleksej Saushev, The NetBSD Foundation
 *
 *   Copyright © 2012      Blue Brain Project, EPFL. All rights reserved.
 *
 *   Copyright © 2015      Research Organization for Information Science
 *   and Technology (RIST). All rights reserved.
 *
 *   Copyright © 2015-2016 Intel, Inc.  All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3. The name of the author may not be used to endorse or promote
 *      products derived from this software without specific prior
 *      written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 *   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *   ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 *   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 *   IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 *   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 *   IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *   }}}
 */
int hwloc_aux_get_depth_from_string(hwloc_topology_t topology, const char * desc)
{
#if HWLOC_API_VERSION >= 0x00020000
    int depth = 0;
    hwloc_type_sscanf_as_depth(desc, NULL, topology, &depth);
    return depth;
#else
    hwloc_obj_type_t type;
    int depthattr;
    hwloc_obj_cache_type_t cachetypeattr;

    int rc = hwloc_obj_type_sscanf(desc,
            &type,
            &depthattr,
            &cachetypeattr,
            sizeof(cachetypeattr));

    if (rc < 0)
        return HWLOC_TYPE_DEPTH_UNKNOWN;

    if (depthattr == -1) {
        /* matched a type without depth attribute, try to get the
         * depth from the type if it exists and is unique */
        return hwloc_get_type_or_above_depth(topology, type);
        /* Error cases are:
         * 
         * HWLOC_TYPE_DEPTH_MULTIPLE
         *   as provided, this type exists at multiple depths in the
         *   hierarchy.
         *
         * HWLOC_TYPE_DEPTH_UNKNOWN
         *   this does not exist. (e.g. we asked for a numanode and
         *   there isn't one in this topology).
         *
         */

    } else if (type == HWLOC_OBJ_GROUP) {
        /* matched a type with a depth attribute, look at the
         * first object of each level to find the depth */
        for(int i=0; i < (int) hwloc_topology_get_depth(topology) ; i++) {
            hwloc_obj_t obj = hwloc_get_obj_by_depth(topology, i, 0);
            ASSERT_ALWAYS(obj != NULL);
            if (obj->type == type
                    && (unsigned) depthattr == obj->attr->group.depth)
                return i;
        }
    } else if (type == HWLOC_OBJ_CACHE) {
        return hwloc_get_cache_type_depth(topology, depthattr, cachetypeattr);
        /* Same error cases as above */
    }
    return HWLOC_TYPE_DEPTH_UNKNOWN;
#endif
}
