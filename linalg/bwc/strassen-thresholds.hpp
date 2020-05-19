#ifndef STRASSEN_THRESHOLDS_H_
#define STRASSEN_THRESHOLDS_H_

#ifdef  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS
#define C128_STRASSEN_THRESHOLDS_D 3
#define C128_STRASSEN_THRESHOLDS_I 4
#define C128_STRASSEN_THRESHOLDS {      \
         {1312, -1}, {1313, 768}, {1314, 256}, {1315, 256},     \
         {1316, 256}, {1317, 256}, {1318, 256}, {1319, 256},    \
         {2336, -1}, {2337, 768}, {2338, 256}, {2339, 256},     \
         {2340, 256}, {2341, 256}, {2342, 256}, {2343, 256},    \
        }
#else   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
template<> unsigned int foo<gf2x_cantor_fft>::default_selector_data[] = {
        /* D = 3 */
        /* I = 4 */
         [1312] = -1, [1313] = 768, [1314] = 256, [1315] = 256, \
         [1316] = 256, [1317] = 256, [1318] = 256, [1319] = 256,        \
         [2336] = -1, [2337] = 768, [2338] = 256, [2339] = 256, \
         [2340] = 256, [2341] = 256, [2342] = 256, [2343] = 256,        \
        }
#endif   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
#ifdef  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS
#define FAKE_STRASSEN_THRESHOLDS_D 3
#define FAKE_STRASSEN_THRESHOLDS_I 4
#define FAKE_STRASSEN_THRESHOLDS {      \
         {1312, -1}, {1313, 1152}, {1314, 384}, {1315, 128},    \
         {1316, 128}, {1317, 128}, {1318, 128}, {1319, 128},    \
         {2336, -1}, {2337, 640}, {2338, 384}, {2339, 128},     \
         {2340, 128}, {2341, 128}, {2342, 128}, {2343, 128},    \
        }
#else   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
template<> unsigned int foo<gf2x_fake_fft>::default_selector_data[] = {
        /* D = 3 */
        /* I = 4 */
         [1312] = -1, [1313] = 1152, [1314] = 384, [1315] = 128,        \
         [1316] = 128, [1317] = 128, [1318] = 128, [1319] = 128,        \
         [2336] = -1, [2337] = 640, [2338] = 384, [2339] = 128, \
         [2340] = 128, [2341] = 128, [2342] = 128, [2343] = 128,        \
        }
#endif   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
#ifdef  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS
#define GF2X_TFFT_STRASSEN_THRESHOLDS_D 3
#define GF2X_TFFT_STRASSEN_THRESHOLDS_I 4
#define GF2X_TFFT_STRASSEN_THRESHOLDS { \
         {1312, -1}, {1313, 20736}, {1314, 10368}, {1315, 128}, \
         {1316, 128}, {1317, 128}, {1318, 128}, {1319, 128},    \
         {2336, -1}, {2337, 31104}, {2338, 10368}, {2339, 128}, \
         {2340, 128}, {2341, 128}, {2342, 128}, {2343, 128},    \
        }
#else   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
template<> unsigned int foo<gf2x_ternary_fft>::default_selector_data[] = {
        /* D = 3 */
        /* I = 4 */
         [1312] = -1, [1313] = 20736, [1314] = 10368, [1315] = 128,     \
         [1316] = 128, [1317] = 128, [1318] = 128, [1319] = 128,        \
         [2336] = -1, [2337] = 31104, [2338] = 10368, [2339] = 128,     \
         [2340] = 128, [2341] = 128, [2342] = 128, [2343] = 128,        \
        }
#endif   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */
#endif	/* STRASSEN_THRESHOLDS_H_ */
