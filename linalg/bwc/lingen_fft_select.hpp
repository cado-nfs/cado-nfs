#ifndef LINGEN_FFT_SELECT_HPP_
#define LINGEN_FFT_SELECT_HPP_

#ifdef SELECT_MPFQ_LAYER_u64k1
#include "gf2x-fft.h"
#else
#include "flint-fft/fft.h"
#endif

#endif	/* LINGEN_FFT_SELECT_HPP_ */
