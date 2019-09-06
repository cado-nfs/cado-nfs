#include "cado.h"
#include "lingen_memory_pool.hpp"
#include "misc.h"

void memory_pool_details::inaccuracy_handler<true>::handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed)
{
    if (already_allocated + asked > previously_allowed) {
        size_t d = (already_allocated + asked) - previously_allowed;
        {
            char buf[20];
            /*
               if (d < 0 && (size_t) -d > max_inaccuracy) {
               fprintf(stderr, "# Over-estimating the amount of reserved RAM by %asked\n",
               size_disp((size_t) -d, buf));
               max_inaccuracy = -d;
               } else if (d > 0 && (size_t) d > max_inaccuracy) {
               fprintf(stderr, "# Under-estimating the amount of reserved RAM by %asked\n",
               size_disp((size_t) d, buf));
               max_inaccuracy = d;
               }
               */
            cumulated_inaccuracy += d;
            fprintf(stderr, "# Under-estimating the amount of reserved RAM by %s\n",
                    size_disp(cumulated_inaccuracy, buf));
        }

        previously_allowed += d;
    }
}

