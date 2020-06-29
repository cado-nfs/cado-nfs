#ifndef AREA_H_
#define AREA_H_

/* default parameters for Murphy's E-value */
#define BOUND_F 1e7
#define BOUND_G 5e6
#define AREA    1e16

#ifdef __cplusplus
extern "C" {
#endif

extern double area, bound_f, bound_g;

#ifdef __cplusplus
}
#endif


/* default rootsieve effort */
#define DEFAULT_ROPTEFFORT 5.0

#endif	/* AREA_H_ */
