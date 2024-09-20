#ifndef MEMUSAGE_H_
#define MEMUSAGE_H_

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

#endif /* MEMUSAGE_H_ */
