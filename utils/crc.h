#ifndef CRC_H_
#define CRC_H_

#include <stdint.h>

struct cado_crc_lfsr_s {
    uint32_t t[32];
    int i; 
    int r;
};

typedef struct cado_crc_lfsr_s cado_crc_lfsr[1];
typedef struct cado_crc_lfsr_s * cado_crc_lfsr_ptr;

#ifdef __cplusplus
extern "C" {
#endif

/* This computes the crc32 value of the n bytes pointed to by c.  */
extern uint32_t crc32(const void * c, size_t n);

/* These provide the possibility of computing a checkum in several parts */
extern uint32_t cado_crc_lfsr_turn1(cado_crc_lfsr_ptr, uint32_t);
extern void cado_crc_lfsr_init(cado_crc_lfsr_ptr);
extern void cado_crc_lfsr_clear(cado_crc_lfsr_ptr);
extern uint32_t cado_crc_lfsr_turn(cado_crc_lfsr_ptr, const void *, unsigned int);
extern uint32_t cado_crc_lfsr_turn32_little(cado_crc_lfsr_ptr, const uint32_t *, unsigned int);

#ifdef __cplusplus
}
#endif

#endif	/* CRC_H_ */
