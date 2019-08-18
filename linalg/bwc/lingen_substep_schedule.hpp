#ifndef LINGEN_SUBSTEP_SCHEDULE_HPP_
#define LINGEN_SUBSTEP_SCHEDULE_HPP_

#include <array>

struct lingen_substep_schedule {
    /* output characteristics -- the ones we have to choose */

    /* batch, shrink0, shrink2. Used to be static parameters, now they're
     * dynamic.
     *
     * batch: with batch=b, schedule b times more simultaneous
     *  transforms, so that we can have more parallelism. Costs b times
     *  more RAM
     * shrink0: divide the number of transforms on dimension 0 (m/n for
     *  MP, (m+n)/r for MUL) so as to save memory
     * shrink2: divide the number of transforms on dimension 2 (m+n)/r
     *  for both MP and MUL) so as to save memory
     */

    /* shrink0 and shrink2 are important. Here, we restrict our operation
     * to shrink0*shrink2 multiplications of matrices of size
     * (n0/shrink0) times (n2/shrink2).
     */
    unsigned int shrink0;
    unsigned int shrink2;

    /* We may batch in all three dimensions. At most, batch is {n0, n1,
     * n2}. More generally, if it is {b0,b1,b2}, then the amount of RAM
     * spent for transforms is r * b1 * (b0 + b2) transforms. This brings
     * the benefit of computing b0*b1 transforms simultaneously on A, and
     * b1*b2 transforms simultaneously on B.
     *
     * The outer loop, of size n1/r, is processed b1 steps at a time.
     */
    std::array<unsigned int, 3> batch;

    lingen_substep_schedule() : shrink0(1), shrink2(1), batch {1,1,1}  {}
    lingen_substep_schedule(lingen_substep_schedule const&) = default;

    bool operator<(lingen_substep_schedule const & o) const {
        if (shrink0 < o.shrink0) return true;
        if (shrink0 > o.shrink0) return false;
        if (shrink2 < o.shrink2) return true;
        if (shrink2 > o.shrink2) return false;
        if (batch < o.batch) return true;
        if (batch > o.batch) return false;
        return false;
    }
    bool operator==(lingen_substep_schedule const & o) const {
        return shrink0 == o.shrink0
            && shrink2 == o.shrink2
            && batch == o.batch;
    }
};

#endif	/* LINGEN_SUBSTEP_SCHEDULE_HPP_ */
