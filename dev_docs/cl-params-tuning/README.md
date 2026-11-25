# Tuning the SIQS parameters for class group computations (d80 – d160)

These are the scripts and the model behind `parameters/cl/params.qs.d80` …
`params.qs.d160` (November 2026). They are research scripts, not part of
the build; they are kept here so that the numbers can be reproduced or
redone on another machine.

## What is being optimised

The cost of `cado-nfs.py -qs -cl` is sieving + linear algebra, and the linear
algebra is not the one of factoring:

* each of the `nmatrices = 2` square matrices needs its determinant over Z,
  obtained by CRT from one block-Wiedemann run per modulus;
* the determinant is not "class number times something small": it has about
  **1.9 bits per row** (measured 1.87 – 1.95 on d60/d80 matrices, by CRT and
  by floating-point LU), because the minors of an (N+1)×N relation matrix
  are the class number times the coordinates of its kernel vector;
* so the number of moduli is `1.9 N / moduli_nbits`, each run costs
  `1.5 N` products by a matrix with `170 N` coefficients, and the total is
  **cubic in N**.

Consequences seen in all the experiments:

* the matrix size is governed by the factor base bound `lim0`, hardly by
  `lpb0`: a small `lim0` and heavy oversieving are right (as in the existing
  files), and the best `lim0` doubles every 10 digits;
* `lpb0` is flat between `log2(lim0) + 6` and `+ 9`; three large primes do
  not pay (tested at d80);
* wide moduli are cheaper per bit of determinant (704 bits: −31 % compared to
  256 bits). `moduli_nbits = 960` does not work: the prime is part of the
  workunit file name, which becomes longer than 255 characters;
* `m`, `n` of block Wiedemann and `target_density` (100 – 400) make no
  difference; the sieve threshold and `ncurves0 ≥ 6` lose no relations; the
  bucket thresholds are best left at their defaults.

## Model

1. **Sieve speed** (`timing2.py`): core-seconds per special-q, 2 processes ×
   2 threads on 4 cores, start-up excluded. Per sieved position it is
   `a + b sqrt(lim0)`, with `b` decreasing with `qfac_nfac`.
2. **Yield** (`dickman.py`, `ymodel.py`, `alpha.py`): relations per special-q
   = number of positions × P(a number of size `X` is `lim0`-smooth up to two
   primes below `2^lpb0`), with `X = ln(2^(I-1) sqrt|D|) + 1.83 + 1.17 alpha(D)`
   and `alpha(D) = - sum_p chi_D(p) ln p / (p-1)` (2 treated separately).
   The two constants were fitted on the d80 and d100 samples; the model then
   predicted the d90, d110, d120 measurements within 0 – 15 %.
   **The discriminant matters**: over random prime discriminants alpha is
   spread over about [−1.6, +1.0] (10 % – 90 %), and the sieving time is
   roughly proportional to `exp(0.3 alpha)`. All estimates below are for
   alpha = 0.
3. **Filtering** (`fakeqs.cpp`, `simfilter.py`, `proxy.py`): synthetic
   relations are fed to the real `purge` and `merge-dl`. A synthetic relation
   is a real one whose primes are moved to nearby prime ideals; for sizes
   where real relations are too expensive, real relations of a smaller
   discriminant are padded with factor base primes up to the right size, and
   the number of large primes is drawn from the yield model.
   Checks: d80 real data, 12042 purged rows / 3298 matrix rows, simulated
   12294 / 3357; d80 → d100 proxy within 4 % of the d100 matrix sizes;
   d100 → d120 proxy within 12 %.
4. **Linear algebra** (`opt.py`, `labench*.py`): 10.4 ns per matrix
   coefficient per iteration with 704-bit moduli (5.5 ns with 256 bits), flat
   from N = 4000 to N = 40000 on the test machine.

`final.py` minimises sieving + linear algebra over the number of relations
for a grid of `lim0`, `lpb0`; `reeval.py` re-evaluates the stored simulations
for other sieve geometries and prints the table; `writeparams.py` writes the
parameter files from `choices.json`.

## Results

Core time on the test machine (4 cores of a 2.1 GHz Xeon, Emerald Rapids
generation), typical discriminant (alpha = 0):

| digits | I | qfac_nfac | lim0 | lpb0 | rels_wanted | special-q | matrix rows | sieving | lin. alg. | total (core-days) |
|---|---|---|---|---|---|---|---|---|---|---|
| 80 | 13 | 13 | 200000 | 26 | 640000 | 8e+03 | 4200 | 0.014 | 0.012 | 0.026 |
| 90 | 13 | 14 | 400000 | 26 | 806000 | 2.9e+04 | 7100 | 0.128 | 0.072 | 0.2 |
| 100 | 13 | 16 | 800000 | 27 | 1540000 | 5.6e+04 | 12000 | 0.99 | 0.33 | 1.3 |
| 110 | 13 | 16 | 1600000 | 28 | 2920000 | 3.9e+05 | 23500 | 8.4 | 2.3 | 10.7 |
| 120 | 15 | 16 | 3200000 | 29 | 5660000 | 8.7e+05 | 44000 | 54 | 15 | 69 |
| 130 | 15 | 16 | 6400000 | 30 | 9000000 | 4.2e+06 | 90000 | 350 | 140 | 490 |
| 140 | 15 | 16 | 12800000 | 31 | 14700000 | 1.9e+07 | 190000 | 1 900 | 1 200 | 3 100 |
| 150 | 15 | 16 | 25600000 | 32 | 28100000 | 9.8e+07 | 340000 | 14 700 | 6 700 | 21 400 |
| 160 | 15 | 16 | 51200000 | 33 | 53900000 | 4.6e+08 | 620000 | 97 600 | 41 400 | 139 000 |

d80 and d90 were run to completion with these files (see below); for d100 –
d120 the sieving statistics were measured and the filtering simulated;
d130 – d160 are extrapolations (sieve speed measured, yield and filtering
from the model). Expect ±30 % on the measured sizes and a factor 2 on the
extrapolated ones, before the effect of alpha.

## Validation runs (4 cores, complete computations)

| run | CPU seconds | sieving | matrix rows | moduli per matrix | lin. alg. |
|---|---|---|---|---|---|
| d80, previous file (lim0 200000, lpb0 22, 256-bit moduli) | 4588 | 3844 | 3298 | 27 | 740 |
| d80, new file | 3024 | 2019 | 4269 | 14 | 966 |
| d90, new file | 15979 | 10430 | 7698 | 23 | 5390 |

The d80 test discriminant has alpha = +2.1 (only 31 splits among the small
primes), the d90 one alpha = −1.0. For the d80 run the model had predicted
2126 core-seconds of sieving, 17264 purged rows and 4177 matrix rows
(measured: 2223 core-seconds occupied, 17557, 4212 – 4269). For d90 it had
predicted 0.128 core-days of sieving and 7100 rows (measured 0.133 and 7698);
that run found 27.8 relations per special-q where 33.8 were expected, because
the special-q's were 3 – 4 times smaller than the optimum for this
discriminant (see below).

In the same model the previous files cost 0.032, 0.235 and 1.47 core-days
at 80, 90 and 100 digits for a typical discriminant, against 0.026, 0.20
and 1.3 for the new ones: they were already close to the optimum.

## Things that would change the picture

* The special-q prime range (`qfac_min`, `qfac_max`) is fixed per parameter
  file, whereas the right range depends on the discriminant: the primes
  actually used should be centred on `(sqrt|D| / 2^(I-1))^(1/qfac_nfac)`.
  With 14 – 16 factors a 5 % shift of the primes is a factor 2 on the
  special-q; the d90 validation run lost about 18 % of yield this way.
  Computing the range from the discriminant in `cado-nfs.py` would fix it.
* The cubic linear algebra. Anything that gives the class number from fewer
  than `1.9 N / moduli_nbits` solves per matrix would move the optimum
  towards larger factor bases and cheaper sieving.
* For N above about 10^5 the vectors no longer fit in cache on CPUs with a
  small L3 per core; the linear algebra figures for d140 – d160 were not
  measured at that size.

## Reproducing

    export CADO_BUILD=/path/to/build/dir
    g++ -O2 -o fakeqs fakeqs.cpp
    python3 timing2.py 120 15 16 29 3200000        # sieve speed, d120 file
    python3 study.py sample d=100 I=13 nfac=16 lim=51200000 lpb=35 nq=24   # templates
    python3 final.py 120 13 18 samples/<template>.json lims=1600000,3200000,6400000 lpbs=29
    python3 reeval.py

`timing2.py` on the d120 line above gives 5.4 core-seconds per special-q on
the test machine; the ratio on another machine is the factor to apply to the
sieving column.
