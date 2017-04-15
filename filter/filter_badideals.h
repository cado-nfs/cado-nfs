#ifndef __FILTER_BADIDEALS_H__
#define __FILTER_BADIDEALS_H__

/* The information from this structure is the one found in the
 * .badidealinfo filename
 */

// Let's assume we will never have polynomials of degree > 10.
#define MAX_DEG 10

typedef struct {
    p_r_values_t p;
    p_r_values_t r;
    unsigned int k;
    int side;
    p_r_values_t pk;
    p_r_values_t rk;
    unsigned int ncol;
    int val[MAX_DEG];
} badid_info_struct_t;

typedef struct {
    int n;
    badid_info_struct_t * badid_info;
} allbad_info_struct_t;

typedef allbad_info_struct_t allbad_info_t[1];

void read_bad_ideals_info(const char *filename, allbad_info_t info);
void handle_bad_ideals (int *exp_above, int64_t a, int64_t b, p_r_values_t p,
                        int e, int side, allbad_info_t info);

#endif   /* __FILTER_BADIDEALS_H__ */
