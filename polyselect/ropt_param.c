/**
 * @file ropt_param.c
 * Holds extern vars. The customizable parameters 
 * are in ropt_param.h.
 */


#include "cado.h" // IWYU pragma: keep
#include <math.h>       // sqrt log
#include "ropt_arith.h" // ROPT_NPRIMES
#include "ropt_param.h"


/**
 * Default values for L1 cache size and tune sieve length (which should be half of
 * the L1 cache size). The code will also try to auto-detect them.
 */
unsigned int L1_cachesize = 16384;
unsigned int size_tune_sievearray = 8192;


/* ------------------
   possible to change
   ------------------ */


/**
 * Total number of sublattices in the tuning and sieving steps.
 * With the default parameters (SIZE_SIEVEARRAY_V_MAX 4194304),
 * each root sieve takes about 2-4 seconds. The first column is `
 * ranked by the digits of integers to be factored and the right
 * column is the number (actually number+1).
 *
 * Usually, there is one or more tuning steps before the final root
 * sieve. In that case, more sublattices are checked.
 * 
 * The number is linearly scaled by param->effort, but a larger value
 * is not always necessary since the sieving is done in order with
 * the best sublattices first.
 *
 * It is organized as {digits, A, B, C}
 *
 * C is number of sublattices used in tuning the best lognorm bound;
 * A is number of sublattices used in tuning;
 * B number of sublattices used in final sieve;
 * In general, A is large; while B, C are small 
 *  (perhaps C could also be large).
 * They will be scaled by ropt_effort linearly.
 * Note the total running-time is scaled approx. linearly since 
 * each lats in B take more time but there are less lats in B
 * than others.
 */
const unsigned int size_total_sublattices[NUM_DEFAULT_DIGITS][4] = {
  { 80,    8,    4,    2},  /* up to 79 digits */
  {100,   16,    8,    3},  /* up to 99 digits */
  {120,   32,   12,    4},  /* up to 119 digits */
  {140,   64,   16,    5},  /* up to 139 digits */
  {150,   96,   20,    6},  /* up to 149 digits */
  {160,  128,   24,    7},  /* up to 159 digits */
  {170,  160,   28,    8},  /* up to 169 digits */
  {180,  192,   32,    9},  /* up to 179 digits */
  {190,  224,   36,    10},  /* up to 189 digits */
  {200,  256,   40,    11},  /* up to 199 digits */
  {210,  288,   44,    12},  /* up to 209 digits */
  {220,  320,   48,    13},  /* up to 219 digits */
  {230,  368,   52,    14},  /* up to 229 digits */
  {240,  416,   56,    15},  /* up to 239 digits */
  {250,  464,   60,    16},  /* up to 249 digits */
  {260,  512,   64,    17},  /* up to 259 digits */
  {270,  560,   68,    18},  /* up to 269 digits */
  {280,  608,   72,    19},  /* up to 279 digits */
  {290,  656,   76,    20},  /* up to 289 digits */
  {300,  704,   80,    21}   /* up to 299 digits */
};


/**
 * Number of top sublattice for individual primes[i] in stage 1,
 * where i < NUM_SUBLATTICE_PRIMES. The construction should depend
 * on the total num of primes in s1param->e_sl[]. 
 * The main purpose is to prevent too much crt computations in stage 1.
 * They will be passed to s1param->individual_nbest_sl[] later.
 */
const unsigned int
s1_size_each_sublattice[NUM_SUBLATTICE_PRIMES][NUM_SUBLATTICE_PRIMES] = {
  { 64,  0,  0,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 1
  { 64, 64,  0,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 2
  { 64, 64, 64,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 3
  { 64, 32, 32, 32,  0,  0,  0,  0,  0 }, // tlen_e_sl = 4
  { 32, 32, 16,  8,  8,  0,  0,  0,  0 }, // tlen_e_sl = 5
  { 32, 16, 16,  8,  8,  4,  0,  0,  0 }, // tlen_e_sl = 6 (current bound)
  { 32, 16,  8,  8,  4,  4,  4,  0,  0 }, // tlen_e_sl = 7
  { 32, 16,  8,  8,  4,  4,  4,  2,  0 }, // tlen_e_sl = 8
  { 32, 16,  8,  8,  4,  4,  4,  2,  2 }, // tlen_e_sl = 9
};


/**
 * As above, but it's only used for tuning good w in 'ropt_quadratic.c'.
 * Therefore, the values are much smaller than above.
 */
const unsigned int
s1_size_each_sublattice_tune[NUM_SUBLATTICE_PRIMES] = {
  16,  8,  8,  4,  4,  2,  2,  2,  2
};


/* -------------------------
   perhaps no need to change
   ------------------------- */


/**
 * Default parameters for sublattice p_i^e_i.
 * Non-decreasing due to ropt_s1param_setup() in ropt_str.c
 */
const unsigned int 
default_sublattice_pe[NUM_DEFAULT_SUBLATTICE][NUM_SUBLATTICE_PRIMES] = {
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23 */
  { 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // 2
  { 2, 0, 0, 0, 0, 0, 0, 0, 0 }, // 4
  { 2, 1, 0, 0, 0, 0, 0, 0, 0 }, // 12
  { 3, 1, 0, 0, 0, 0, 0, 0, 0 }, // 24
  { 4, 1, 0, 0, 0, 0, 0, 0, 0 }, // 48
  { 3, 2, 0, 0, 0, 0, 0, 0, 0 }, // 72
  { 3, 1, 1, 0, 0, 0, 0, 0, 0 }, // 120
  { 4, 2, 0, 0, 0, 0, 0, 0, 0 }, // 144
  { 4, 1, 1, 0, 0, 0, 0, 0, 0 }, // 240
  { 5, 2, 0, 0, 0, 0, 0, 0, 0 }, // 288
  { 3, 2, 1, 0, 0, 0, 0, 0, 0 }, // 360
  { 5, 1, 1, 0, 0, 0, 0, 0, 0 }, // 480
  { 4, 2, 1, 0, 0, 0, 0, 0, 0 }, // 720
  { 5, 2, 1, 0, 0, 0, 0, 0, 0 }, // 1440
  { 5, 3, 1, 0, 0, 0, 0, 0, 0 }, // 4320
  { 5, 2, 1, 1, 0, 0, 0, 0, 0 }, // 10080
  { 6, 2, 1, 1, 0, 0, 0, 0, 0 }, // 20160
  { 5, 3, 2, 0, 0, 0, 0, 0, 0 }, // 21600
  { 5, 2, 2, 1, 0, 0, 0, 0, 0 }, // 50400
  { 6, 2, 2, 1, 0, 0, 0, 0, 0 }, // 100800
  { 6, 3, 2, 1, 0, 0, 0, 0, 0 }, // 302400
  { 6, 3, 2, 1, 1, 0, 0, 0, 0 }, // 3326400
  { 7, 3, 2, 1, 1, 0, 0, 0, 0 }, // 6652800
  { 8, 3, 2, 1, 1, 0, 0, 0, 0 }, // 13305600
  { 7, 4, 2, 1, 1, 0, 0, 0, 0 }, // 19958400
  { 6, 3, 2, 1, 1, 1, 0, 0, 0 }, // 43243200
  { 7, 3, 2, 1, 1, 1, 0, 0, 0 }, // 86486400
  { 7, 4, 2, 1, 1, 1, 0, 0, 0 }, // 259459200
};


/**
 * Hard-coded product for above default_sublattice_pe.
 */
const unsigned long
default_sublattice_prod[NUM_DEFAULT_SUBLATTICE] = {
  2, 4, 12, 24, 48, 72, 120, 144, 240, 288, 360, 480,
  720, 1440, 4320, 10080, 20160, 21600, 50400, 100800,
  302400, 3326400, 6652800, 13305600, 19958400, 43243200,
  86486400, 259459200 };


/**
 * Primes.
 */
const unsigned int primes[ROPT_NPRIMES] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199
};


/**
 * next_prime_idx[3] = ind(5) = 2 in above table.
 */
const unsigned char next_prime_idx[] = {
  0, // 0
  0, 1, 2, 2, 3, 3, 4, 4, 4, 4, // 1 - 10
  5, 5, 6, 6, 6, 6, 7, 7, 8, 8,
  8, 8, 9, 9, 9, 9, 9, 9, 10, 10,
  11, 11, 11, 11, 11, 11, 12, 12, 12, 12,
  13, 13, 14, 14, 14, 14, 15, 15, 15, 15,
  15, 15, 16, 16, 16, 16, 16, 16, 17, 17,
  18, 18, 18, 18, 18, 18, 19, 19, 19, 19,
  20, 20, 21, 21, 21, 21, 21, 21, 22, 22,
  22, 22, 23, 23, 23, 23, 23, 23, 24, 24,
  24, 24, 24, 24, 24, 24, 25, 25, 25, 25,
  26, 26, 27, 27, 27, 27, 28, 28, 29, 29,
  29, 29, 30, 30, 30, 30, 30, 30, 30, 30,
  30, 30, 30, 30, 30, 30, 31, 31, 31, 31,
  32, 32, 32, 32, 32, 32, 33, 33, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 35, 35,
  36, 36, 36, 36, 36, 36, 37, 37, 37, 37,
  37, 37, 38, 38, 38, 38, 39, 39, 39, 39,
  39, 39, 40, 40, 40, 40, 40, 40, 41, 41,
  42, 42, 42, 42, 42, 42, 42, 42, 42, 42,
  43, 43, 44, 44, 44, 44, 45, 45 // 191 - 198
};



