#ifndef POLYSELECT_SPECIAL_Q_H_
#define POLYSELECT_SPECIAL_Q_H_

#ifdef __cplusplus
extern "C" {
#endif


#if ULONG_MAX == 4294967295UL
#define LEN_SPECIAL_Q 57
#else
#define LEN_SPECIAL_Q 59
#endif
//#define DEBUG_HASH_TABLE
extern const unsigned int SPECIAL_Q[LEN_SPECIAL_Q];

/* declarations */

extern const unsigned int SPECIAL_Q[];

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_SPECIAL_Q_H_ */

