#ifndef CADO_ARITH_MODP_CARRY_SAVE_HPP
#define CADO_ARITH_MODP_CARRY_SAVE_HPP

#include "arith-modp-main.hpp"

#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
#include <x86intrin.h>
#endif

namespace arith_modp::details {
#if 0
    /* This is old AVX/SSE carry-save code. It hasn't been tested in a
     * long while, and it has never been put in production as far as I
     * can tell.
     */
    /* 
     * We have
     * little use for this code at the moment since it is suboptimal.
     */
#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
    template<> struct fast_type<gfp<3, 1> > {
        typedef gfp<3, 1> super;
        struct elt;
        typedef elt elt_ur_for_add;
        struct elt {
#ifdef HAVE_AVX2
        static const int alignment = 32;
#else
        static const int alignment = 16;
#endif
            typedef elt self;
#ifdef HAVE_AVX2
            __m256i data[1];
#else
            __m128i data[2];
#endif
            elt() { zero(); }
            elt(elt const& a) = default;
            elt& operator=(elt const& a) = default;

            /* we do not construct (nor affect) from mpz, because we're not
             * positional */
            void zero() {
#ifdef HAVE_AVX2
                data[0] = _mm256_setzero_si256();
#else
                data[0] = _mm_setzero_si128();
                data[1] = _mm_setzero_si128();
#endif
            }
            static void zero(elt * x, int N) {
                // see bug 21663
                memset(x->data, 0, N * sizeof(data));
            }
            static void copy(elt * y, const elt * x, int N) {
                // see bug 21663
                memcpy(y->data, x->data, N * sizeof(data));
            }
            bool operator==(elt const& a) {
                return memcmp(data, a.data, sizeof(data)) == 0;
            }
            elt(super::elt const& a) {
                convert(*this, a);
            }

            operator super::elt_ur_for_add() const {
                super::elt_ur_for_add carries(conv_backend_get_carries(*this));
                super::add(carries, conv_backend_get_main(*this));
                return carries;
            }

            /* same, but we assume carry is zero */
            operator super::elt() const {
                ASSERT(conv_backend_get_carries(*this).is_zero());
                return conv_backend_get_main(*this);
            }
        };

        static inline void stream_store(elt * dst, elt const& src) {
            /* Do we want to stream that or not ? In fact it's slower
             * when streaming... */
#if 0
#ifdef HAVE_AVX2
            _mm256_stream_si256(dst->data+0, src.data[0]);
#else
            _mm_stream_si128(dst->data+0, src.data[0]);
            _mm_stream_si128(dst->data+1, src.data[1]);
#endif
#else
#ifdef HAVE_AVX2
            _mm256_storeu_si256(dst->data+0, src.data[0]);
#else
            _mm_storeu_si128(dst->data+0, src.data[0]);
            _mm_storeu_si128(dst->data+1, src.data[1]);
#endif
#endif
        }
        static inline void add(elt & dst, elt const & src)
        {
#ifdef HAVE_AVX2
            dst.data[0] = _mm256_add_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_add_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_add_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub(elt & dst, elt const & src)
        {
#ifdef HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub_ur(elt_ur_for_add & dst, elt_ur_for_add const & src)
        {
#ifdef HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        /* conversions are done as a combination of blend & shuffle */

#ifdef HAVE_AVX2
        /* We grok only values for w_i which are integer immediates
         * within {-1} \cup {0..15}
         */
#define shuffle_16bit_words(                                                   \
  out, in, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf)     \
    do {                                                                       \
        *out = _mm256_xor_si256(                                               \
          _mm256_shuffle_epi8(                                                 \
            in,                                                                \
            _mm256_setr_epi8(                                                  \
              (w0 < 0) ? -1 : ((w0 < 8) ? 2 * (w0 & 7) + 0 : -1),              \
              (w0 < 0) ? -1 : ((w0 < 8) ? 2 * (w0 & 7) + 1 : -1),              \
              (w1 < 0) ? -1 : ((w1 < 8) ? 2 * (w1 & 7) + 0 : -1),              \
              (w1 < 0) ? -1 : ((w1 < 8) ? 2 * (w1 & 7) + 1 : -1),              \
              (w2 < 0) ? -1 : ((w2 < 8) ? 2 * (w2 & 7) + 0 : -1),              \
              (w2 < 0) ? -1 : ((w2 < 8) ? 2 * (w2 & 7) + 1 : -1),              \
              (w3 < 0) ? -1 : ((w3 < 8) ? 2 * (w3 & 7) + 0 : -1),              \
              (w3 < 0) ? -1 : ((w3 < 8) ? 2 * (w3 & 7) + 1 : -1),              \
              (w4 < 0) ? -1 : ((w4 < 8) ? 2 * (w4 & 7) + 0 : -1),              \
              (w4 < 0) ? -1 : ((w4 < 8) ? 2 * (w4 & 7) + 1 : -1),              \
              (w5 < 0) ? -1 : ((w5 < 8) ? 2 * (w5 & 7) + 0 : -1),              \
              (w5 < 0) ? -1 : ((w5 < 8) ? 2 * (w5 & 7) + 1 : -1),              \
              (w6 < 0) ? -1 : ((w6 < 8) ? 2 * (w6 & 7) + 0 : -1),              \
              (w6 < 0) ? -1 : ((w6 < 8) ? 2 * (w6 & 7) + 1 : -1),              \
              (w7 < 0) ? -1 : ((w7 < 8) ? 2 * (w7 & 7) + 0 : -1),              \
              (w7 < 0) ? -1 : ((w7 < 8) ? 2 * (w7 & 7) + 1 : -1),              \
              (w8 < 0) ? -1 : ((w8 >= 8) ? 2 * (w8 & 7) + 0 : -1),             \
              (w8 < 0) ? -1 : ((w8 >= 8) ? 2 * (w8 & 7) + 1 : -1),             \
              (w9 < 0) ? -1 : ((w9 >= 8) ? 2 * (w9 & 7) + 0 : -1),             \
              (w9 < 0) ? -1 : ((w9 >= 8) ? 2 * (w9 & 7) + 1 : -1),             \
              (wa < 0) ? -1 : ((wa >= 8) ? 2 * (wa & 7) + 0 : -1),             \
              (wa < 0) ? -1 : ((wa >= 8) ? 2 * (wa & 7) + 1 : -1),             \
              (wb < 0) ? -1 : ((wb >= 8) ? 2 * (wb & 7) + 0 : -1),             \
              (wb < 0) ? -1 : ((wb >= 8) ? 2 * (wb & 7) + 1 : -1),             \
              (wc < 0) ? -1 : ((wc >= 8) ? 2 * (wc & 7) + 0 : -1),             \
              (wc < 0) ? -1 : ((wc >= 8) ? 2 * (wc & 7) + 1 : -1),             \
              (wd < 0) ? -1 : ((wd >= 8) ? 2 * (wd & 7) + 0 : -1),             \
              (wd < 0) ? -1 : ((wd >= 8) ? 2 * (wd & 7) + 1 : -1),             \
              (we < 0) ? -1 : ((we >= 8) ? 2 * (we & 7) + 0 : -1),             \
              (we < 0) ? -1 : ((we >= 8) ? 2 * (we & 7) + 1 : -1),             \
              (wf < 0) ? -1 : ((wf >= 8) ? 2 * (wf & 7) + 0 : -1),             \
              (wf < 0) ? -1 : ((wf >= 8) ? 2 * (wf & 7) + 1 : -1))),           \
          _mm256_shuffle_epi8(/* 0x4E is 0b01001110 aka (1,0,3,2) */           \
                              _mm256_permute4x64_epi64(                        \
                                in, _MM_SHUFFLE(1, 0, 3, 2)),                  \
                              _mm256_setr_epi8(                                \
                                (w0 < 0)                                       \
                                  ? -1                                         \
                                  : ((w0 >= 8) ? 2 * (w0 & 7) + 0 : -1),       \
                                (w0 < 0)                                       \
                                  ? -1                                         \
                                  : ((w0 >= 8) ? 2 * (w0 & 7) + 1 : -1),       \
                                (w1 < 0)                                       \
                                  ? -1                                         \
                                  : ((w1 >= 8) ? 2 * (w1 & 7) + 0 : -1),       \
                                (w1 < 0)                                       \
                                  ? -1                                         \
                                  : ((w1 >= 8) ? 2 * (w1 & 7) + 1 : -1),       \
                                (w2 < 0)                                       \
                                  ? -1                                         \
                                  : ((w2 >= 8) ? 2 * (w2 & 7) + 0 : -1),       \
                                (w2 < 0)                                       \
                                  ? -1                                         \
                                  : ((w2 >= 8) ? 2 * (w2 & 7) + 1 : -1),       \
                                (w3 < 0)                                       \
                                  ? -1                                         \
                                  : ((w3 >= 8) ? 2 * (w3 & 7) + 0 : -1),       \
                                (w3 < 0)                                       \
                                  ? -1                                         \
                                  : ((w3 >= 8) ? 2 * (w3 & 7) + 1 : -1),       \
                                (w4 < 0)                                       \
                                  ? -1                                         \
                                  : ((w4 >= 8) ? 2 * (w4 & 7) + 0 : -1),       \
                                (w4 < 0)                                       \
                                  ? -1                                         \
                                  : ((w4 >= 8) ? 2 * (w4 & 7) + 1 : -1),       \
                                (w5 < 0)                                       \
                                  ? -1                                         \
                                  : ((w5 >= 8) ? 2 * (w5 & 7) + 0 : -1),       \
                                (w5 < 0)                                       \
                                  ? -1                                         \
                                  : ((w5 >= 8) ? 2 * (w5 & 7) + 1 : -1),       \
                                (w6 < 0)                                       \
                                  ? -1                                         \
                                  : ((w6 >= 8) ? 2 * (w6 & 7) + 0 : -1),       \
                                (w6 < 0)                                       \
                                  ? -1                                         \
                                  : ((w6 >= 8) ? 2 * (w6 & 7) + 1 : -1),       \
                                (w7 < 0)                                       \
                                  ? -1                                         \
                                  : ((w7 >= 8) ? 2 * (w7 & 7) + 0 : -1),       \
                                (w7 < 0)                                       \
                                  ? -1                                         \
                                  : ((w7 >= 8) ? 2 * (w7 & 7) + 1 : -1),       \
                                (w8 < 0) ? -1                                  \
                                         : ((w8 < 8) ? 2 * (w8 & 7) + 0 : -1), \
                                (w8 < 0) ? -1                                  \
                                         : ((w8 < 8) ? 2 * (w8 & 7) + 1 : -1), \
                                (w9 < 0) ? -1                                  \
                                         : ((w9 < 8) ? 2 * (w9 & 7) + 0 : -1), \
                                (w9 < 0) ? -1                                  \
                                         : ((w9 < 8) ? 2 * (w9 & 7) + 1 : -1), \
                                (wa < 0) ? -1                                  \
                                         : ((wa < 8) ? 2 * (wa & 7) + 0 : -1), \
                                (wa < 0) ? -1                                  \
                                         : ((wa < 8) ? 2 * (wa & 7) + 1 : -1), \
                                (wb < 0) ? -1                                  \
                                         : ((wb < 8) ? 2 * (wb & 7) + 0 : -1), \
                                (wb < 0) ? -1                                  \
                                         : ((wb < 8) ? 2 * (wb & 7) + 1 : -1), \
                                (wc < 0) ? -1                                  \
                                         : ((wc < 8) ? 2 * (wc & 7) + 0 : -1), \
                                (wc < 0) ? -1                                  \
                                         : ((wc < 8) ? 2 * (wc & 7) + 1 : -1), \
                                (wd < 0) ? -1                                  \
                                         : ((wd < 8) ? 2 * (wd & 7) + 0 : -1), \
                                (wd < 0) ? -1                                  \
                                         : ((wd < 8) ? 2 * (wd & 7) + 1 : -1), \
                                (we < 0) ? -1                                  \
                                         : ((we < 8) ? 2 * (we & 7) + 0 : -1), \
                                (we < 0) ? -1                                  \
                                         : ((we < 8) ? 2 * (we & 7) + 1 : -1), \
                                (wf < 0) ? -1                                  \
                                         : ((wf < 8) ? 2 * (wf & 7) + 0 : -1), \
                                (wf < 0)                                       \
                                  ? -1                                         \
                                  : ((wf < 8) ? 2 * (wf & 7) + 1 : -1))));     \
    } while (0)
#else
#define shuffle_16bit_words(                                                   \
  out, lo, hi, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf) \
    do {                                                                       \
        out[0] = _mm_xor_si128(                                                \
          _mm_shuffle_epi8(                                                    \
            lo,                                                                \
            _mm_setr_epi8((w0 < 0) ? -1 : ((w0 < 8) ? 2 * (w0 & 7) + 0 : -1),  \
                          (w0 < 0) ? -1 : ((w0 < 8) ? 2 * (w0 & 7) + 1 : -1),  \
                          (w1 < 0) ? -1 : ((w1 < 8) ? 2 * (w1 & 7) + 0 : -1),  \
                          (w1 < 0) ? -1 : ((w1 < 8) ? 2 * (w1 & 7) + 1 : -1),  \
                          (w2 < 0) ? -1 : ((w2 < 8) ? 2 * (w2 & 7) + 0 : -1),  \
                          (w2 < 0) ? -1 : ((w2 < 8) ? 2 * (w2 & 7) + 1 : -1),  \
                          (w3 < 0) ? -1 : ((w3 < 8) ? 2 * (w3 & 7) + 0 : -1),  \
                          (w3 < 0) ? -1 : ((w3 < 8) ? 2 * (w3 & 7) + 1 : -1),  \
                          (w4 < 0) ? -1 : ((w4 < 8) ? 2 * (w4 & 7) + 0 : -1),  \
                          (w4 < 0) ? -1 : ((w4 < 8) ? 2 * (w4 & 7) + 1 : -1),  \
                          (w5 < 0) ? -1 : ((w5 < 8) ? 2 * (w5 & 7) + 0 : -1),  \
                          (w5 < 0) ? -1 : ((w5 < 8) ? 2 * (w5 & 7) + 1 : -1),  \
                          (w6 < 0) ? -1 : ((w6 < 8) ? 2 * (w6 & 7) + 0 : -1),  \
                          (w6 < 0) ? -1 : ((w6 < 8) ? 2 * (w6 & 7) + 1 : -1),  \
                          (w7 < 0) ? -1 : ((w7 < 8) ? 2 * (w7 & 7) + 0 : -1),  \
                          (w7 < 0) ? -1                                        \
                                   : ((w7 < 8) ? 2 * (w7 & 7) + 1 : -1))),     \
          _mm_shuffle_epi8(                                                    \
            hi,                                                                \
            _mm_setr_epi8((w0 < 0) ? -1 : ((w0 >= 8) ? 2 * (w0 & 7) + 0 : -1), \
                          (w0 < 0) ? -1 : ((w0 >= 8) ? 2 * (w0 & 7) + 1 : -1), \
                          (w1 < 0) ? -1 : ((w1 >= 8) ? 2 * (w1 & 7) + 0 : -1), \
                          (w1 < 0) ? -1 : ((w1 >= 8) ? 2 * (w1 & 7) + 1 : -1), \
                          (w2 < 0) ? -1 : ((w2 >= 8) ? 2 * (w2 & 7) + 0 : -1), \
                          (w2 < 0) ? -1 : ((w2 >= 8) ? 2 * (w2 & 7) + 1 : -1), \
                          (w3 < 0) ? -1 : ((w3 >= 8) ? 2 * (w3 & 7) + 0 : -1), \
                          (w3 < 0) ? -1 : ((w3 >= 8) ? 2 * (w3 & 7) + 1 : -1), \
                          (w4 < 0) ? -1 : ((w4 >= 8) ? 2 * (w4 & 7) + 0 : -1), \
                          (w4 < 0) ? -1 : ((w4 >= 8) ? 2 * (w4 & 7) + 1 : -1), \
                          (w5 < 0) ? -1 : ((w5 >= 8) ? 2 * (w5 & 7) + 0 : -1), \
                          (w5 < 0) ? -1 : ((w5 >= 8) ? 2 * (w5 & 7) + 1 : -1), \
                          (w6 < 0) ? -1 : ((w6 >= 8) ? 2 * (w6 & 7) + 0 : -1), \
                          (w6 < 0) ? -1 : ((w6 >= 8) ? 2 * (w6 & 7) + 1 : -1), \
                          (w7 < 0) ? -1 : ((w7 >= 8) ? 2 * (w7 & 7) + 0 : -1), \
                          (w7 < 0) ? -1                                        \
                                   : ((w7 >= 8) ? 2 * (w7 & 7) + 1 : -1))));   \
        out[1] = _mm_xor_si128(                                                \
          _mm_shuffle_epi8(                                                    \
            lo,                                                                \
            _mm_setr_epi8((w8 < 0) ? -1 : ((w8 < 8) ? 2 * (w8 & 7) + 0 : -1),  \
                          (w8 < 0) ? -1 : ((w8 < 8) ? 2 * (w8 & 7) + 1 : -1),  \
                          (w9 < 0) ? -1 : ((w9 < 8) ? 2 * (w9 & 7) + 0 : -1),  \
                          (w9 < 0) ? -1 : ((w9 < 8) ? 2 * (w9 & 7) + 1 : -1),  \
                          (wa < 0) ? -1 : ((wa < 8) ? 2 * (wa & 7) + 0 : -1),  \
                          (wa < 0) ? -1 : ((wa < 8) ? 2 * (wa & 7) + 1 : -1),  \
                          (wb < 0) ? -1 : ((wb < 8) ? 2 * (wb & 7) + 0 : -1),  \
                          (wb < 0) ? -1 : ((wb < 8) ? 2 * (wb & 7) + 1 : -1),  \
                          (wc < 0) ? -1 : ((wc < 8) ? 2 * (wc & 7) + 0 : -1),  \
                          (wc < 0) ? -1 : ((wc < 8) ? 2 * (wc & 7) + 1 : -1),  \
                          (wd < 0) ? -1 : ((wd < 8) ? 2 * (wd & 7) + 0 : -1),  \
                          (wd < 0) ? -1 : ((wd < 8) ? 2 * (wd & 7) + 1 : -1),  \
                          (we < 0) ? -1 : ((we < 8) ? 2 * (we & 7) + 0 : -1),  \
                          (we < 0) ? -1 : ((we < 8) ? 2 * (we & 7) + 1 : -1),  \
                          (wf < 0) ? -1 : ((wf < 8) ? 2 * (wf & 7) + 0 : -1),  \
                          (wf < 0) ? -1                                        \
                                   : ((wf < 8) ? 2 * (wf & 7) + 1 : -1))),     \
          _mm_shuffle_epi8(                                                    \
            hi,                                                                \
            _mm_setr_epi8((w8 < 0) ? -1 : ((w8 >= 8) ? 2 * (w8 & 7) + 0 : -1), \
                          (w8 < 0) ? -1 : ((w8 >= 8) ? 2 * (w8 & 7) + 1 : -1), \
                          (w9 < 0) ? -1 : ((w9 >= 8) ? 2 * (w9 & 7) + 0 : -1), \
                          (w9 < 0) ? -1 : ((w9 >= 8) ? 2 * (w9 & 7) + 1 : -1), \
                          (wa < 0) ? -1 : ((wa >= 8) ? 2 * (wa & 7) + 0 : -1), \
                          (wa < 0) ? -1 : ((wa >= 8) ? 2 * (wa & 7) + 1 : -1), \
                          (wb < 0) ? -1 : ((wb >= 8) ? 2 * (wb & 7) + 0 : -1), \
                          (wb < 0) ? -1 : ((wb >= 8) ? 2 * (wb & 7) + 1 : -1), \
                          (wc < 0) ? -1 : ((wc >= 8) ? 2 * (wc & 7) + 0 : -1), \
                          (wc < 0) ? -1 : ((wc >= 8) ? 2 * (wc & 7) + 1 : -1), \
                          (wd < 0) ? -1 : ((wd >= 8) ? 2 * (wd & 7) + 0 : -1), \
                          (wd < 0) ? -1 : ((wd >= 8) ? 2 * (wd & 7) + 1 : -1), \
                          (we < 0) ? -1 : ((we >= 8) ? 2 * (we & 7) + 0 : -1), \
                          (we < 0) ? -1 : ((we >= 8) ? 2 * (we & 7) + 1 : -1), \
                          (wf < 0) ? -1 : ((wf >= 8) ? 2 * (wf & 7) + 0 : -1), \
                          (wf < 0) ? -1                                        \
                                   : ((wf >= 8) ? 2 * (wf & 7) + 1 : -1))));   \
    } while (0)
#endif

        /* case of 192 bits within 256 bits. Three 64-bit words
         * split into four 48-bit words.
         */
        static void convert(elt& dst, const super::elt& a) {
            /* index of 16-bit word in destination, fetched from
             * which index of 16-bit word in the gfp::elt. This is
             * given for the 256-bit registers
             *
             * 0    0
             * 1    1
             * 2    2
             * 3    <empty>
             * 4    3
             * 5    4
             * 6    5
             * 7    <empty>
             * 8    6
             * 9    7
             * 10   8
             * 11   <empty>
             * 12   9
             * 13   10
             * 14   11
             * 15   <empty>
             */
#ifdef HAVE_AVX2
            /* I'm really upset here. _mm256_shuffle_epi8, aka VPSHUFB,
             * reads only 4-byte immediates (and discards the rest). As a
             * consequence, the following does not work: the indices
             * 12,13,14,15 read off bounds, while the 16,17, etc actually
             * do what they want, but based on the fact that they're
             * reduced mod 16 + implicitly considered wrt the high part
             * of the operand...
            dst.data[0] = _mm256_shuffle_epi8(
                    _mm256_loadu_si256((__m256i*) a.x),
                    _mm256_setr_epi8( 
                        0,1,2,3,4,5,-1,-1,
                        6,7,8,9,10,11,-1,-1,
                        12,13,14,15,16,17,-1,-1,
                        18,19,20,21,22,23,-1,-1));
            */

#if 0
            __m256i in = _mm256_loadu_si256((__m256i*) a.x);
            dst.data[0] =
                    _mm256_xor_si256(
                        _mm256_shuffle_epi8(
                        in,
                        _mm256_setr_epi8( 
                            0,1,2,3,4,5,-1,-1,
                            6,7,8,9,10,11,-1,-1,
                            -1,-1,-1,-1,0,1,-1,-1,
                            2,3,4,5,6,7,-1,-1)),
                        _mm256_shuffle_epi8(
                            /* 0x4E is 0b01001110 aka (1,0,3,2) */
                            _mm256_permute4x64_epi64 (in, _MM_SHUFFLE(1,0,3,2)), 
                        _mm256_setr_epi8( 
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            12,13,14,15,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1)));
#endif
            __m256i in = _mm256_loadu_si256((__m256i*) a.x);
            shuffle_16bit_words(dst.data, in,
                    0,1,2,-1, 3,4,5,-1, 6,7,8,-1, 9,10,11,-1);
#else /* SSSE3 !! */
            __m128i lo = _mm_loadu_si128((__m128i*) a.x);
            __m128i hi = _mm_loadu_si128((__m128i*) (a.x + 2));
            shuffle_16bit_words(dst.data, lo, hi,
                    0,1,2,-1, 3,4,5,-1, 6,7,8,-1, 9,10,11,-1);
            /* note that 16bit-wide shuffles use an 8-bit immediate,
             * but do not offer the option to selectively insert
             * zeroes. So we're probably better off shuffling bytes.
             */
#endif
        }

        static super::elt conv_backend_get_main(elt const& src) {
            /* This is needed because we knowingly write off bounds */
            union station {
                super::elt e;
#ifdef HAVE_AVX2
                __m256i v[1];
#else
                __m128i v[2];
#endif
                station() {};
            } main;
#ifdef HAVE_AVX2
            shuffle_16bit_words(main.v, src.data[0],
                    0,1,2, 4,5,6, 8,9,10, 12,13,14, -1,-1,-1,-1);
#else
            shuffle_16bit_words(main.v, src.data[0], src.data[1],
                    0,1,2, 4,5,6, 8,9,10, 12,13,14, -1,-1,-1,-1);
#endif
            return main.e;
        }
        static super::elt_ur_for_add conv_backend_get_carries(elt const& src) {
            union station {
                super::elt_ur_for_add e;
#ifdef HAVE_AVX2
                __m256i v[1];
#else
                __m128i v[2];
#endif
                station() {};
            } carries, ncarries;

            /* It's slightly more complicated than it seems. The carry
             * words may be negative. So we must sign-extend them to the
             * full unreduced element size.
             */
#ifdef HAVE_AVX2
            shuffle_16bit_words(carries.v, src.data[0],
                    -1,-1,-1,3,
                    -1,-1,7,-1,
                    -1,11,-1,-1,
                    15,-1,-1,-1);
            __m256i zero = _mm256_setzero_si256();
            shuffle_16bit_words(ncarries.v,
                    _mm256_sub_epi16(zero,
                        _mm256_cmpgt_epi16(zero, carries.v[0])),
                    -1,-1,-1,-1,
                    3,-1,-1,6,
                    -1,-1,9,-1,
                    -1,12,-1,-1);
#else
            shuffle_16bit_words(carries.v, src.data[0], src.data[1],
                    -1,-1,-1,3,
                    -1,-1,7,-1,
                    -1,11,-1,-1,
                    15,-1,-1,-1);
            __m128i zero = _mm_setzero_si128();
            shuffle_16bit_words(ncarries.v,
                    _mm_sub_epi16(zero, _mm_cmpgt_epi16(zero, carries.v[0])),
                    _mm_sub_epi16(zero, _mm_cmpgt_epi16(zero, carries.v[1])),
                    -1,-1,-1,-1,
                    3,-1,-1,6,
                    -1,-1,9,-1,
                    -1,12,-1,-1);
#endif
            super::sub_ur(carries.e, ncarries.e);
            return carries.e;
        }



        /* (add|sub)mul_ui go through convert, do naively and convert
         * back. Yes, it's slightly painful. Here we assume that src
         * has undergone little to no accumulated additions, so that
         * it can basically be converted lossless to a gfp::elt
         *
         * This prototype was changed recently in the "main"
         * implementation. Ideally, the way to go would be to make the
         * version below a member function rather than a static function
         * that receives the prime and preinv.
         */
        static inline void addmul_ui(elt & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            ARITH_MODP_TEMPORARY_ALLOC(T, elt, zr);
            ARITH_MODP_TEMPORARY_ALLOC(T, elt_ur_for_add, z);
            z = dst;
            super::addmul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }
        static inline void submul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            ARITH_MODP_TEMPORARY_ALLOC(T, elt, zr);
            ARITH_MODP_TEMPORARY_ALLOC(T, elt_ur_for_add, z);
            z = dst;
            super::submul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }

        /* we have *TWO* reduction functions here. One which assigns to a
         * standard gfp::elt, and one which assigns to fast_type::elt */
        static void reduce(super::elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            ARITH_MODP_TEMPORARY_ALLOC(T, elt_ur_for_add, z);
            z = a;
            super::reduce(r, z, p, j);
        }
        static void reduce(elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            ARITH_MODP_TEMPORARY_ALLOC(T, elt, zr);
            reduce(zr, a, p, j);
            r = zr;
        }
    };
#endif /* defined(HAVE_AVX2) || defined(HAVE_SSSE3) */
#endif
}

#endif /* ARITH_MODP_CARRY_SAVE_HPP_ */
