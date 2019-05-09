#ifndef LINGEN_SUBSTEP_SCHEDULE_H_
#define LINGEN_SUBSTEP_SCHEDULE_H_

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

    /* batch is not used much, because it goes in the direction of
     * spending _more_ memory, which is rarely what we're interested in.
     * What this does is that several steps of the outer loop (off size
     * n1/r) are done simultaneously: more transforms are kept live.
     */
    unsigned int batch;

#ifdef __cplusplus
    lingen_substep_schedule() : shrink0(1), shrink2(1), batch(1) {}
    lingen_substep_schedule(lingen_substep_schedule const&) = default;
#endif
};

#endif	/* LINGEN_SUBSTEP_SCHEDULE_H_ */
