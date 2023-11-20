#ifndef LINGEN_FFT_SELECT_HPP_
#define LINGEN_FFT_SELECT_HPP_

#ifdef LINGEN_BINARY
#include "gf2x-fft.h"   // IWYU pragma: export
#include "gf2x-fake-fft.h"   // IWYU pragma: export
#include "gf2x-cantor-fft.h"   // IWYU pragma: export
#include "gf2x-ternary-fft.h"   // IWYU pragma: export
#else
#include "flint-fft/fft.h"      // IWYU pragma: export
#include "flint-fft/transform_interface.h" // IWYU pragma: export
#endif

#endif	/* LINGEN_FFT_SELECT_HPP_ */
