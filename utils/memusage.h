#ifndef MEMUSAGE_H_
#define MEMUSAGE_H_

// IWYU pragma: private, include "utils.h"

#ifdef	__cplusplus
extern "C" {
#endif

extern long Memusage (void);
extern long Memusage2 (void);
extern long PeakMemusage (void);

#ifdef	__cplusplus
}
#endif

#endif /* MEMUSAGE_H_ */
