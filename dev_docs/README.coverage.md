
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
    COV=1 make cmake
    COV=1 make -j8
    eval $(COV=1 make show)
    C=coverage-report
    commit=$(git rev-parse --short=8 HEAD)
    commit_ref="https://gitlab.inria.fr/cado-nfs/cado-nfs/-/commit/$commit"
    (cd $build_tree ; time gcovr --merge-mode-functions=separate -r $src_tree --json $src_tree/${C}-base.json)
    COV=1 ./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
    (cd $build_tree ; time gcovr --merge-mode-functions=separate -r $src_tree --json $src_tree/${C}-app.json)
    rm -rf coverage 2>/dev/null || : ; mkdir coverage
    gcovr --merge-mode-functions=separate -a ${C}-base.json -a ${C}-app.json --html-title "Coverage for commit <a href=\"$commit_ref\">$commit</a>" --html-nested coverage/coverage.html --print-summary
    xdg-open coverage/coverage.html
```

How to visualize the automatic coverage reports
-----------------------------------------------

Simple: follow this link: [![coverage report](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/coverage.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/file/coverage/index.html?job=merge+coverage+tests) (for the master branch).
