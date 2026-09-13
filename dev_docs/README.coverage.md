Coverage tests
==============

Two things are covered here.
 - how to do coverage tests by hand locally
 - how to visualize the coverage reports that are run automatically by
   the cado-nfs pipelines on the gitlab server.

How to do coverage tests by hand locally
----------------------------------------

I add this to my `local.sh` file:
```
if [ "$COV" ] ; then
    build_tree="${build_tree}.cov"
    DEBUG=1
    CFLAGS="-O0 -g --coverage -fprofile-update=atomic"
    CXXFLAGS="-O0 -g --coverage -fprofile-update=atomic"
    LDFLAGS="--coverage"
fi
```

Two programs are needed: `grcov`, which reads the counter files that the
instrumented binaries leave behind, and `genhtml`, which comes with
`lcov`. On debian, `apt install grcov lcov`; elsewhere, grcov is a single
statically linked binary, see https://github.com/mozilla/grcov/releases.

Then:

```
    COV=1 make -j8 all all_test_dependencies && COV=1 scripts/cov.sh -R numbertheory
```

In the line above, `-R numbertheory` can be replaced by anything you see
fit and that restricts the test breadth in order to focus on what you're
working on exactly. Anything you pass is handed over to ctest, so `-E`
and `-j8` work too. You can also run `cov.sh` alone in order to run all
tests, but be warned that this takes a while.

The report lands in `$build_tree/coverage/index.html` and `cov.sh` prints
the path (and opens it, if `gio` is around). Note that the percentages
are always relative to the whole project, even when you only ran fifteen
tests -- a file that no test touched is reported at 0%, not omitted.

`scripts/cov.sh` sources `ci/ci.bash` and calls the very same function
that the pipeline calls, so what you see locally is what the pipeline
will show.

How the pipeline does it
------------------------

Worth knowing if you ever have to touch it, because the split of
responsibilities is deliberate:

 - the `coverage tests N/M` jobs run their share of the test suite and
   then do nothing but `tar` up the `.gcda` files that the compiler
   runtime wrote. That costs a fraction of a second and a couple of
   megabytes of artifacts. One of them, `1/M`, also ships the `.gcno`
   files, which describe the shape of the code and are identical in every
   job of the pipeline.
 - every coverage job compiles with the same `-frandom-seed` (the
   commit, see `ci/ci/005-build-environment.sh`). Without it gcc stamps
   each build with its own timestamp, the `.gcno` of one job pairs with
   the `.gcda` of no other, and gcov reports 0% rather than complaining.
 - `merge coverage tests` combines the counters of every job with
   `gcov-tool merge` -- which sums the counters themselves, so the result
   is exactly what one job running the whole suite would have produced --
   and renders the report once, with grcov and genhtml.

The list of third-party paths that are kept out of the reports is
`coverage_excluded_paths` in `ci/ci.bash`, and it is the only place where
that list lives.

Where to go from there
----------------------

 - four strides is a guess, not a measurement. The fixed per-job cost
   that the old thirteen-way split was working around is gone, so the
   split is now about test time alone, and one job may well be enough.
   Read the timings of a pipeline or two before settling.
 - grcov is the one thing the pipeline fetches from outside (a github
   release binary, pinned and checksummed in `ci/ci.sh`). Mirroring it
   as a gitlab.inria.fr project upload, the way `ecm-7.0.6.tar.gz`
   already is, would remove that.
 - do not swap grcov for lcov at the capture step without reinstating a
   zero-coverage baseline pass (`lcov -c -i`). `lcov --capture` only
   reports a compilation unit that has a `.gcda`, so files that no test
   loaded vanish from the report instead of showing at 0%; on our tree
   that is about 190 files and it turns a truthful 47% into a
   meaningless 62%. gcovr and grcov do not have that problem.
 - clang source-based coverage (`-fprofile-instr-generate
   -fcoverage-mapping`) would give region and branch coverage rather
   than line coverage, and would drop the `-fprofile-update=atomic`
   slowdown that `TIMEOUT_SCALE=2` currently papers over. The catch is
   that `llvm-cov` needs the binaries at report time, and the coverage
   build tree is 1.5G.

How to visualize the automatic coverage reports
-----------------------------------------------

Simple: follow this link: [![coverage report](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/coverage.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/file/coverage/index.html?job=merge+coverage+tests) (for the master branch).

Coverage of the lines that a merge request actually touches is shown in
the merge request diff itself, which gitlab builds from the cobertura
report that the same job publishes.
