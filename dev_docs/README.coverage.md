
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
    CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
    CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
fi
```

And then the following commands produce an html report, provided the required tools are installed (notably, `lcov`).

```
    COV=1 make cmake
    COV=1 make -j8
    COV=1 ./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
    eval $(COV=1 make show)
    geninfo --ignore-errors gcov,source -q --output-filename /tmp/app.info -b . $build_tree --no-external
    ./ci/utilities/coverage_local_infofile_modifications.pl -d $build_tree /tmp/app.info
    ./ci/utilities/coverage-postprocess.sh -o /tmp/coverage /tmp/app.info
    xdg-open /tmp/coverage/index.html
```

How to visualize the automatic coverage reports
-----------------------------------------------

Simple: follow this link: [![coverage report](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/coverage.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/file/coverage/index.html?job=merge+coverage+tests) (for the master branch).
