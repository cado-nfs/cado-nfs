
Short of a hitchhiker's guide to cado && cmake.

cmake is a source tree configurator and makefile generator. The two most
important things to know are:

- cmake's configuration files are CMakeLists.txt files.
- cmake builds _out of the source tree_. This means that the built object
  files and binaries don't go in the same directory as the source files.

Tweaking CMakeLists.txt is relatively easy using a follow-by-example
strategy. Lacking sufficient information in the already existing
CMakeLists.txt files, the usual sources of information with respect to
cmake are given by
  * cmake --help-commands | less
  * cmake --help-variables | less
  * cmake --help-properties | less
  * all the ``See Also'' items given in the output of the commands above
  * the CMake book in office A213

This being said, here are the peculiarities of the cmake+cado combination:

* Makefile ; call_cmake.sh

In order to avoid breaking users' old habits, it is possible to live
almost without noticing cmake is behind the scenes. From any source
directory within cado, it is possible to run "make <some target name>",
and there is some black magic supposed to do Just The Right Thing (tm).
In particular, there is a Makefile hanging around whose job is to hand
over the responsibility of building the target to cmake. In particular,
this initiates the build tree (see item *local.sh below), and configures
it by calling cmake. The magic script is called scripts/call_cmake.sh.

Note that in case cmake is not available on the host computer, the script
tries to automatically wget and install it.

This way of building requires some default settings to be honoured as to
where the build tree goes. The default setting is to use a directory
called build/<hostname> from the top of the cado source tree. However
this can bee customized (see item *local.sh below)

Passing variables to makefiles is no longer possible with command-line
overrides, or more precisely it may or may not to what you want. It is
more portable to use environment variables overrides, like:

CFLAGS=-g make

Note also that a tree which has been built partially with one set of flags may
require being wiped out before using another set of flags. To do so, use the
"make tidy" command (see item *special targets below). A more precise way to
act calls the "make cmake" command (see item *special targets below)

The set of environment variables understood by cmake is restricted. See
the top-level CMakeLists.txt file for the current list.

See local.sh.example for a documentation on the recognized flags.

* local.sh

At the top of the source tree, if a file is called local.sh exists, it is
sourced by scripts/call_cmake.sh ; therefore, any environment variable
which is set by this script is passed to cmake.

In addition to the environment variables looked up by cmake in
CMakeLists.txt, it is also possible to override the default setting for
the build tree location. For this, set the build_tree environment
variable. It must be a relative path with respect to the current
directory, which may be anywhere within the cado source tree. Mimick
scripts/call_cmake.sh, and use the variables defined there in order to
set something sensible.

* CMakeLists-nodist.txt

Several sub-directories contain a CMakeLists-nodist.txt file. This file
lists the targets which are only built for the development tree, but not
relevant for the distributed tree. Therefore, these CMakeLists-nodist.txt
files are not distributed either (see the item *make dist below).

* special targets

Several new meta-targets are defined.

- make tidy ; this wipes out the build tree
- make show ; this is primarily of use to some scripts (the output can be
  sourced by a shell).
- make dist ; see item *make dist below.
- make install ; do not use it. It is not really relevant, and not
  reasonably tested.
- make cmake ; force cmake to regenerate the makefiles, possibly gaining
  knowledge of newly set environment variables.
- make clean ; cleans the objects and binaries in the subdir of the build
  tree which corresponds to the current directory.

* make dist

This target expectably build a cado-nfs-<version>.tar.gz file. This
archive name may be forced to an arbitrary user-defined value by setting
the environment variable CADO_DIST_ARCHIVE_NAME.

What goes in the distribution, what doesn't ?

This is controlled by the files called files.dist, files.nodist, and
files.unknown. These three files must form a set of patterns matching
each file of the distribution exactly once. make dist will abort if it is
not the case. ``pattern'' here doesn't mean anything rich here. It is
either a full path to a file, or a directory name ending with a slash ;
in the latter case, it means ``all files below this directory''. 

If ``make dist'' does not work for you because someone else has committed
a new file, and you don't know whether it should go into the distribution or
not, then update the files.unknown file, and blame the person responsible
for the mess. The normal situation is when files.unknown is empty.

Note: the tools required for building a distribution are not distributed.

* nightly-test

This script is used for checking the distribution.

* run_example.sh

This script supersedes the old new_run.c59 script and friends. An
invocation of new_run.c59 thus becomes run_example.sh c59

# vim: set ft=conf:
