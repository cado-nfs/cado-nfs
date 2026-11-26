# Handoff: cl SIQS parameters, to be continued from a machine with cluster access

Written 2026-11-26 at the end of a cloud session that could not reach
Grid'5000 (no credentials, no route). Everything below was measured on that
session's machine: 4 cores of a 2.1 GHz Xeon (Emerald Rapids class), 15 GB.
Nothing has run on grdix or grvingt yet.

## State

Branch `claude/bold-bell-zgnx30`, on top of `692ecb7`:

* `parameters/cl/params.qs.d80` … `d160` (d80, d90, d100 retuned; d110 – d160
  new), registered in `files.dist`.
* this directory: scripts, model, `README.md` (method, table, validation).
* `final/*.jsonl`: the filtering simulations behind the table (matrix size
  against number of relations, per size and `lim0`). `reeval.py` reads them.
* `samples/d120…rels.gz`: 3513 real 120-digit relations (the only sample
  that is slow to regenerate: about 45 core-minutes).
* `timings.json`: measured core-seconds per special-q, per geometry and `lim0`.

How far each size was checked:

| size | status |
|---|---|
| d80 | complete run with the committed file: 3024 CPU s (old file: 4588), same class group |
| d90 | complete run: 15979 CPU s, model said 0.20 core-days |
| d100, d110 | yield and sieve speed measured, filtering simulated; never run |
| d120 | same, plus a 5-minute start of the real file (5.5 CPU s per special-q, 6.9 relations each, 2.7 GB) |
| d130 – d160 | sieve speed measured, yield and filtering from the model; never started |

## What to do first on the cluster side

1. **Calibrate** on one grdix node and one grvingt node (a few minutes each;
   the script loads every core, 1.4 GB per core):

       export CADO_BUILD=/path/to/build/<host>
       python3 calibrate.py

   It prints one factor for the sieving column and one for the linear
   algebra column, relative to the reference machine (where it prints
   about 1.0 and 1.0), and writes `calibration.<host>.json`.

2. **Recompute the table** with the measured factors:

       python3 reeval.py sscale=<sieve factor> lascale=<lin. alg. factor>

   This also re-optimises the number of relations; if the two factors differ
   a lot the best `rels_wanted` moves (cheaper linear algebra → fewer
   relations), and `choices.json` + `writeparams.py` regenerate the files.
   Until then the cluster figures quoted with the table assume 1.0 for a
   grdix core and 1.5 for a grvingt core; these were guesses.

3. **Run d100, then d110, to completion** on grdix (about 1.3 and 11
   reference core-days). These are the first sizes that were never run.
   Compare with `README.md` / `choices.json`:

   | quantity | d100 expected | d110 expected |
   |---|---|---|
   | relations per special-q (alpha = 0) | 27 | 7.4 |
   | core-s per special-q | 1.5 | 1.8 |
   | duplicates | 8 % assumed | 8 % assumed |
   | matrix rows | 12 000 | 23 500 |
   | moduli per matrix | 36 | 68 |

   A matrix more than about 15 % off, or a duplicate rate far from the
   assumption, is the signal to redo the simulation for that size
   (`final.py`, see `README.md`).

## For a specific discriminant

* `python3 -c "from alpha import alpha; print(alpha(D))"` with D negative.
  The sieving time scales roughly like `exp(0.27 alpha)`; the table is for
  alpha = 0.
* Centre the special-q primes on `g = (sqrt|D| / 2^(I-1))^(1/qfac_nfac)`:
  the files assume |D| ≈ 10^(d − 1/2) and an average density of split
  primes. `mkparams.qfac_range` shows the rule used. The d90 run lost 18 % of
  yield from a range that sat too low for its discriminant.
* Parameters can be overridden as `key=value` right after the discriminant
  on the `cado-nfs.py` command line. With `--parameters <file>` no local
  client is started and the run just waits.

## Known weak points

* **Duplicates**: 6 – 7 % measured at d80/d90; 8 – 10 % assumed above. It
  grows with the number of special-q's and was not measurable in small
  samples.
* **`rels_wanted` at d140 – d160** is the first point of a coarse ladder
  (steps of 1.4) with positive excess, i.e. barely enough; `cado-nfs.py` tops
  up by 1 % steps if short. The best `lim0` at d160 was at the edge of the
  grid (51.2M); 102.4M was not simulated.
* **Linear algebra above 40 000 rows** was not measured. The d140 – d160
  matrices (190 000 – 620 000 rows, 17 – 55 MB per vector) will not fit the
  per-core L3 of an EPYC 9754.
* **Proxy relations** (smaller-discriminant relations padded up) gave matrix
  sizes within 4 % at d80 → d100 and about 12 % low at d100 → d120.
* **Start-up costs for large `lpb0`**: `freerel` took 74 CPU s for
  `lpb0 = 29`; the renumber table for `lpb0 = 33` is 16 times larger and
  was not tried.
* **`moduli_nbits`**: 704 is the largest value that works; 960 fails when
  the client uploads a result ("File name too long", the prime is in the
  workunit name).
* `replay-dl` hit `addRowsUpdateIndex: Assertion e1 != 0 && e2 != 0` on a
  merge history produced from *synthetic* relations (single-threaded too).
  It never happened on real relations, and was not investigated; noted in
  case it shows up on a real run.
* Not tried: three large primes above 80 digits; the papers of Kleinjung
  (Math. Comp. 2016) and CSI-FiSh were unreachable, so nothing is compared
  with their parameters.

## Build notes

Ubuntu 24.04 needed `libgmp-dev libhwloc-dev libfmt-dev libecm-dev` and the
Python `flask` module before `make cmake` passed; gmp-ecm is required for
the `factor` binary used at the end of a class group run.
