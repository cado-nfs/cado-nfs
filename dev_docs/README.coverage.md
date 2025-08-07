
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

And then the following commands produce an html report in $PWD/coverage,
provided gcovr is installed.

```
    COV=1 make -j8 all all_test_dependencies && COV=1 scripts/cov.sh -R numbertheory
```

In the line above, `-R numbertheory` can be replaced by anything you see
fit and that restricts the test breadth in order to focus on what you're
working on exactly. You can also run `cov.sh` alone, or perhaps `cov.sh
-j8`, in order to run all tests.

How to visualize the automatic coverage reports
-----------------------------------------------

Simple: follow this link: [![coverage report](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/coverage.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/file/coverage/lcov/index.html?job=merge+coverage+tests) (for the master branch).
