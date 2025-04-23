#ifndef CADO_MEMUSAGE_H
#define CADO_MEMUSAGE_H

#include <stddef.h>

#ifdef	__cplusplus
extern "C" {
#endif

extern size_t Memusage (void);
extern size_t Memusage2 (void);
extern size_t PeakMemusage (void);

#ifdef	__cplusplus
}
#endif

#endif /* CADO_MEMUSAGE_H */
