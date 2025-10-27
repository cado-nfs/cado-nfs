
The main page of the Cado-NFS source code is
[https://gitlab.inria.fr/cado-nfs/cado-nfs](https://gitlab.inria.fr/cado-nfs/cado-nfs).
If you're accessing the cado-nfs source from a different link, it may be
an outdated fork. (This being said, all commits to the
`master` branch are automatically mirrored to the [cado-nfs project on
GitHub](https://github.com/cado-nfs/cado-nfs), so the latter should be
up-to-date as well.)

[![pipeline status](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/pipeline.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/pipelines?ref=master)
[![coverage report](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/coverage.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/file/coverage/lcov/index.html?job=merge+coverage+tests)
[![coverity scan](https://scan.coverity.com/projects/23184/badge.svg)](https://scan.coverity.com/projects/cado-nfs)

Quick install:
==============

(see also the section about [running cado-nfs in a
container](#containers))

in most cases the following should work to factor the number 903...693
using all cores on the local machine

1. `make`
2. `./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693`

More details follow.

Important note: For a larger factorization (distributed on several
machines), the command line to be used is a priori more involved. Please
refer to [`scripts/cadofactor/README.md`](scripts/cadofactor/README.md).
Documented example parameter files are in
[`parameters/factor/params.c90`](parameters/factor/params.c90) and
[`scripts/cadofactor/parameters*`](scripts/cadofactor/).

Supported platforms
===================

The primary development platform is `x86_64` linux with gcc 10 or newer,
the most common processor being Intel core2-like or more recent.

Other architectures are checked regularly, and should work. Please refer
to the gitlab-ci page for the list of regularly tested platforms, and their
current status. The overall pipeline status for the master branch is
[![pipeline
status](https://gitlab.inria.fr/cado-nfs/cado-nfs/badges/master/pipeline.svg)](https://gitlab.inria.fr/cado-nfs/cado-nfs/-/pipelines?ref=master),
and details can be obtained by clicking on the badges.  Note however that
a failing pipeline might mean that a bug affects only one platform in
particular, or could be caused by one runner encountering difficulties.
Such things do happen.

Since we use gitlab-ci pipelines, the authoritative source as to what
gets tested on a routine basis is of course the [`.gitlab-ci.yml`
file](.gitlab-ci.yml). In plain english, the configurations that we test,
or that we at least used to test at some point, are as follows.
Anything beyond the set of regularly tested machines perhaps works,
perhaps does not.

 * The primary development platform is `x86_64` Debian GNU/Linux, latest
   version, with gcc. If it doesn't work, we have a problem.

 * The current version, as well as a few other versions of Debian, Alpine
   Linux, Fedora, Centos Stream, and OpenSuse are also tested. CentOS
   distributions (or derivatives) that have been EOL'd for some time, or
   are deemed to be shortly because they have vastly out-of-date
   software, are not tested.

 * Current FreeBSD is also routinely tested.

 * `x86_64` with icc used to work, but we switched to icx later on. Now
   that cado-nfs requires C++20, icx 2023 or later is required.  Routine
   checks use the icx version that is provided by Intel's
   [`intel/oneapi-hpckit:latest` docker
   image](https://hub.docker.com/r/intel/oneapi-hpckit).

 * Mac OS X is used for CI routine compilation checks, using the
   default compiler from XCode. All version from 10.5 onwards were part
   of the CI routine checks at some point in time, and worked
   successfully. However we obviously do not commit to continued support
   for old versions. Since C++20 is required, you need a fairly recent
   XCode installation to compile cado-nfs.

 * Windows used to be partly supported, but this has been abandoned for
   some time now (see a longer note
   [there](#using-cado-nfs-under-windows) at the end of this file).

These configurations are run within specific containers, or on specific
machines. Compared to the base install, these are equipped with the
necessary dependencies (see below). The console outputs for the different
builds contain information related to the compiler versions being used.


Required software tools
=======================

 * [GMP](https://gmplib.org/), version 5 or newer: usually installed in
   most Linux distributions (on some Linux distributions you need to
   install the `libgmp-dev` or `gmp-devel` package that includes
   `gmp.h`. It is often not installed by default). Note: make sure to
   configure GMP with `--enable-shared` so that a shared library is
   installed (`libgmp.so` under Linux) otherwise CADO-NFS might not
   compile.
 * As of `cado-nfs-3.0.0`, a C/C++ compiler and C/C++ standard library that
   conform to the C99 and C++20 standards are required.
   * GCC: the minimal required version is >= 10
   * LLVM Clang: the minimal required version is >= 12.0.0
   * Apple Clang: the minimal required version is >= 16.0.0
   * Intel ICX: the minimal required version is >= 2023
 * GNU make and CMake (`cmake 3.18` or later) for building.
 * Python 3.8 or later is required, as well as a few fairly common
   packages such as python3-requests and python3-flask. These are
   packaged with most software distributions, or alternatively you can
   install them via pip.  As a convenience means, we also provide the script
   [`./scripts/setup-venv.sh`](scripts/setup-venv.sh), which you can use
   to install the python requirements of cado-nfs in a venv.  Running
   this script is probably the quickest way to get you going on the
   Python front.
 * Several very common unix tools are used in many places. Among them,
   `gzip`, and under some circumstances, `ssh` and `rsync`.
 * On MacOS X, the cado-nfs client script needs an alternative to the
   system-provided curl binary, which it doesn't like. The simplest way
   to deal with this issue is to install the wget downloader (e.g. via
   homebrew).
 * For a large computation, it is recommended to use a MySQL backend for
   the database.

Optionally, cado-nfs can use the following additional software.

* Support for OpenMP (at least version 3.0)
* Support for MPI (see [`local.sh.example`](local.sh.example) and [`linalg/bwc/README`](linalg/bwc/README))
* Support for hwloc (see [`parameters/misc/cpubinding.conf`](parameters/misc/cpubinding.conf))
  Under Debian the command to install HWLOC is:
  	apt-get install libhwloc-dev
* Support for GMP-ECM. Define the environment variable GMPECM if it is
  installed in a non-standard place.
* A formatting library [`fmt`](https://fmt.dev/) is used if found, otherwise a
  snapshot is embedded in the cado-nfs code anyway.

Configure
=========

Normally you have nothing to do to configure cado-nfs.

However if your site needs tweaks, set such tweaks using environment variables,
or via the shell script local.sh ; you may start with

```
cp local.sh.example local.sh
```

Edit according to your local settings and your taste: local.sh.example
gives documentation on the recognized environment variables and their effect.

Note that tweaks in local.sh are global (relevant to all sub-directories
of the cado-nfs tree, not to a particular directory).

As a rule of thumb, whenever you happen to modify the `local.sh` script, it
is advised to trigger re-configuration of the build system, by the
special command `make cmake`. Another way to achieve the same goal is
to remove the build tree, which is below `build/` (normally named as the
hostname of the current machine): `make tidy` should do it.
Then, of course, you must recompile with `make`, since `make cmake`
is just the equivalent of `./configure` in a classical build system.

Optional (alternative): configure using cmake directly
======================================================

cado-nfs includes a top-level `Makefile` which builds the binary objects in
a special build directory which is below `build/` (normally named as the
hostname of the current machine). This way, parallel builds for different
machines may coexist in one shared directory. This is sort of a magic
`out-of-source` build.

Another way to do `out-of-source` build is to create a separate build
directory and build from there, by calling cmake directly for the
configuration. This proceeds as follows:

```
mkdir /some/build/directory
cd /some/build/directory
cmake /path/to/cado-nfs/source/tree
```

Note however that in this case, the `local.sh` file found in the source
tree is not read (but you may elect to do so before calling cmake).

Compile:
========

```
make
```

Install:
========

The relevance of the `make install` step depends on the platform.
Cado-nfs binaries link to shared libraries, and some do so dynamically.
For this to work, we rely on some control logic by cmake and cooperation
with the operating system's dynamic linker.

* if `make` is directly called from the source directory `$SRCDIR`,
  then `make install` installs all programs and binaries in `$SRCDIR/installed`.

* otherwise programs and binaries will be installed in
  `/usr/local/share/cado-nfs-x.y.z`, and this default installation prefix
  can be changed by one of the following commands:

```
cmake .../cado-nfs -DCMAKE_INSTALL_PREFIX=/tmp/install
export PREFIX=/tmp/install ; cmake .../cado-nfs
```

There are several ways to call cado-nfs scripts (e.g., `cado-nfs.py`).
Here we assume that `$SRCDIR` is the source directory, that `$BUILDDIR`
is the build tree for the local machine (typically `$SRCDIR/build/$(hostname)`),
and that `$PREFIX` is the installation prefix (see above).  We refer to
these different ways, and later discuss how they work on different
systems (which is mostly impacted by the shared library mechanism).

* `$SRCDIR/cado-nfs.py`
  This deduces `$BUILDDIR` from the machine hostname, and amounts to
  calling binaries from there. Parameter files are obtained from
  `$SRCDIR/parameters/`

* `$PREFIX/bin/cado-nfs.py`
  This calls binaries from `$PREFIX/bin`, and loads parameter files from
  `$PREFIX/share/cado-nfs-x.y.z/`

* `$BUILDDIR/cado-nfs.py`
  This calls binaries within `$BUILDDIR`. This is useful when several
  processors with the same architecture share a file-system, you have
  to compile CADO-NFS on one only.

Linux, BSD: the first two choices above should work ok.
MacOS X:
    For any invocation to work, the `LD_LIBRARY_PATH` or `DYLD_LIBRARY_PATH`
    variable must be set up properly. The easiest method is to do make
    install, and include in these environment variables the directory
    `$PREFIX/lib/cado-nfs-x.y.z`.

Run a factorization on the current machine:
===========================================

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 -t 2
```

where the option `-t 2` tells how many cores (via threads) to use on the
current machine (for polynomial selection, sieving, linear algebra, among
others).  It is possible to set `-t all` (which, in fact, is the default)
to use all threads on the current machine.

CADO-NFS is optimized only for numbers above 85 digits, and no support will
be provided for numbers of less than 60 digits (for very large numbers,
no support is promised). Note that it is a good idea to remove small prime
factors using special-purpose algorithms such as trial division, P-1, P+1,
or ECM, and use CADO-NFS only for the remaining composite factor.

Parts of the Number Field Sieve computation are massively distributed. In
this mode, client scripts (namely, `cado-nfs-client.py`) are run on many
nodes, connect to a central server, and run programs according to which
computations need to be done.  The programs (for the polynomial selection
and sieving steps) can run multithreaded. It is better to have them
run with a capped number of threads (say, 2), and run several clients per
node. By default, programs in this step are limited to 2 threads. When
running the computation on the local machine, the number of clients is
set so that the number of cores specified by `-t` are kept busy.

Run a larger factorization on several machines:
===============================================

CADO-NFS has several ways of operation, which can be roughly split into
three modes as follows.

 * For small computations, or for tests where it is important to have a
   single command line, the strategy
   [above](#run-a-factorization-on-the-current-machine) works. The
   `cado-nfs.py` script can arrange so that the binaries for all steps of
   the computation are run on the current machine, or possibly on other
   machines, via SSH. Some of the documentation here is specific to this
   mode of operation (see in particular
   [there](#check-that-your-network-configuration-is-correct) or
   [there](#check-that-your-ssh-configuration-is-correct)). If you get it
   right, you may manage to run factorizations as follows. However be
   aware that this mode of operation is fragile, and we advise not to use
   it beyond trivial testing.

   ```
   ./cado-nfs.py 353493749731236273014678071260920590602836471854359705356610427214806564110716801866803409 slaves.hostnames=hostname1,hostname2,hostname3 --slaves 4 --client-threads 2
   ```

   This starts 4 clients per host on the hosts `hostname1`, `hostname2`,
   and `hostname3`, and each client uses two cpus (threads). For
   hostnames that are not `localhost`, ssh is used to connect to the host
   and start a client there.  To configure ssh, see the [next
   section](#check-that-your-ssh-configuration-is-correct). For tasks
   that use the local machine only (not massively distributed tasks), the
   number of threads used is the one given by `-t` (which defaults to all
   threads on the local machine).

 * For larger computations where work distribution is an important point
   (distribution on several machines, possibly with different parameters
   for different machines), it is considerably more flexible to let the
   server be _just_ a server, and start clients separately, that will be
   used to offload the distributed tasks (polynomial selection and
   relation collection). Clients can come and go. The server has plenty
   of timeout provisions to deal with such events.

   This is called the `--server` mode (see
   [`scripts/cadofactor/README.md`](scripts/cadofactor/README.md) and
   [`scripts/cadofactor/parameters`](scripts/cadofactor/parameters)).
   See also [this thread](https://sympa.inria.fr/sympa/arc/cado-nfs/2020-03/msg00001.html)
   on the `cado-nfs` list. If you want to use cado-nfs even to a little
   extent, we recomment that you familiarize with this mode of operation.

 * For much larger computations, the `cado-nfs.py` is only of moderate
   use. The individual cado-nfs binaries and internal scripts are the
   most flexible entry points, and should be used in order to adapt to
   the specificities of the platform being used (e.g. to deal with
   various requirements such as memory for filtering, interconnect for
   linear algebra, and so on).

Check that your network configuration is correct:
=================================================

In case you run a factorization on the local machine, the clients should be able
to connect to the server. Under Linux, CADO-NFS uses 'localhost' to identify the
server, thus you should have the following line in your /etc/hosts file:

```
127.0.0.1   localhost
```

Check that your SSH configuration is correct:
=============================================

The master script (unless in `--server` mode) uses SSH to connect to
available computing resources.  In order to avoid the script asking your
password or passphrase, you must have public-key authentication and an
agent.

The SSH keys are usually installed in the files `~/.ssh/id_rsa` and
`~/.ssh/id_rsa.pub`; if you don't have them yet, you can create them with the
`ssh-keygen` command. See the man page `ssh-keygen(1)` for details. The private
key should be protected with a passphrase, which you can enter when you
create the keys. Normally ssh will ask for the key's passphrase when you log
on to a machine, but this can be avoided by using ssh-agent, see the man
page `ssh-agent(1)`, and providing the passphrase to the agent with `ssh-add`.
Public-key authenticaton together with an ssh-agent will allow `cadofactor`
to use ssh to run commands on slave machines automatically.

Most of the recent Linux distributions will run an `ssh-agent` for you. But
if this is not the case with your distribution, or if you are running
`cado-nfs.py` inside a `screen` in order to logout from your desktop, you
will need to run the `ssh-agent` by hand. As a short recipe, you can type:
```
eval `ssh-agent`
ssh-add
```

You should also copy your public key, i.e., the contents of the file
`~/.ssh/id_rsa.pub`, into `$HOME/.ssh/authorized_keys` on the slave machines, to
allow logging in with public-key authentication.

Also, since localhost has an IP and key that varies, you should have
those 3 lines in your `$HOME/.ssh/config`:

```
Host    localhost
        StrictHostKeyChecking no
        UserKnownHostsFile /dev/null
```

If everything is correctly configured, when you type

ssh localhost

you should end up with a new shell on your machine, without having to
type any password/passphrase.


Restarting an interrupted factorization:
========================================

If you have started a factorization with the `cado-nfs.py` script, and it was
interrupted (for example because of a power failure) you can restart in
any of these two ways:

 * with the same `cado-nfs.py` command line if a work directory was
   explicitly provided on the command line:

   ```
   $ cado-nfs.py ... workdir=/path/to/workdir
   ```

 * with a single argument as in:

   ```
   $ cado-nfs.py     [[work directory]]/XXX.parameters_snapshot.YYY
   ```

   where `[[work directory]]` is the directory which has been chosen
   automatically, `XXX` is the "name" of the running factorisation, and `YYY`
   is the largest possible integer value for which such a file exists.

Factoring with SNFS:
====================

It is possible to take advantage of the special form of an integer and
use the Special Number Field Sieve. See
[`parameters/factor/parameters.F9`](parameters/factor/parameters.F9)
for that:

```
$ cado-nfs.py parameters/factor/parameters.F9 slaves.hostnames=localhost
```

Note in particular that you can force the special-q to be on the rational
side if this is more appropriate for your number, with
`tasks.sieve.sqside=0` on the `cado-nfs.py` command line or in the
parameter file (assuming side 0 is the rational side).

The default square root algorithm does not work in some very rare cases
that could possibly occur with SNFS polynomials (a degree 4 polynomial
with Galois group $Z/2 \times Z/2$ is the only reasonable example, next
case is for degree 8). The CRT approach is a workaround. See
[`sqrt/crtaglsqrt.c`](sqrt/crtaglsqrt.c) .

Big factorization (200 digits and more):
========================================

By default, to decrease memory usage, it is assumed that less than $2^32$
(~ four billion) relations or ideals are needed and that the ideals will
be less than $2^32$ (i.e., the `lpb0` and `lpb1` parameters are less or
equal to 32). In the case of factorizations of numbers of 200 digits and
more, these assumptions may not hold. In this case, you have to set some
variables in your `local.sh` script (see Configure section above for more
information on `local.sh` and section on big factorizations in
`local.sh.example`).

Factoring with two non-linear polynomials:
==========================================

You can find a bit of information on this topic in the development
version, in the GIT repository (see file
[`README.nonlinear`](README.nonlinear)).

Importing polynomials or relations:
===================================

If you have already computed a good polynomial and/or relations, you can
tell CADO-NFS to use them, see
[`scripts/cadofactor/README.md`](scripts/cadofactor/README.md).

Containers
==========

As a current work in progress, cado-nfs now ships pre-prepared containers
which should be enough to do some quick tests. For example, the following
command line should pull pre-compiled cado-nfs binaries, and run them in
a docker container (assuming you are using an `x86_64` CPU, haswell or
later).

```
docker run --rm registry.gitlab.inria.fr/cado-nfs/cado-nfs/factoring-full cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
```

Again, this is work in progress.


Using CADO-NFS under Windows:
=============================

Portability of CADO-NFS on Windows was not an initial goal of that project,
however we give here some hints that might help people wanting to use CADO-NFS
under Windows. We hope that the following information can be useful to
some extent, but the general message is that you're on your own.

* Cado-NFS uses the POSIX interface all over the place, which means that
  in one form of the other, you need to have the corresponding
  functionality. If you don't, you're out of luck. (e.g., cado-nfs cannot
  build as a "pure" visual studio project.)

* if you only need the siever to run on Windows, then you only need to compile
  the `las` program on Windows.

* [Cygwin](http://www.cygwin.com/) provides a Unix-like environment,
  where compilation should be easy.  However the binary requires a
  `cygwin.dll` file. We have been told of problems with shared libraries,
  which the following seems to address:
  ```
  PATH="installed/lib/cado-nfs-x.y.z:$PATH" ./cado-nfs.py [...]
  ```

* if you want a binary without any dependency, you might try
  [MinGW](http://www.mingw.org/). The INSTALL file from GNU MPFR contains
  detailed instructions on how to compile MPFR under Windows. Those
  instructions should work for CADO-NFS too.  See
  [`dev_docs/howto-MinGW.txt`](dev_docs/howto-MinGW.txt).

* you might try to use MPIR (<http://mpir.org/>) instead of GMP. MPIR
  is a fork of GMP, which claims to be more portable under Windows.
  Alternatively, Windows ports of GMP shouldn't be too hard to find.

* you might succeed in compiling the cado-nfs binaries with a
  cross-compiler for Windows (which does not waive the runtime
  requirements for `cado-nfs.py`, notably on unix-like utilities). See
  [`dev_docs/README.cross`](dev_docs/README.cross) in the source code
  repository for information on how to cross-compile.

Examples of basic usage:
========================

* Run a full factorization on the local machine, using all available
  cores:

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
```

* Run a full factorization on the local machine, using all available
  cores, but with a slightly modified set of default parameters.

  The difficulty here is that when `cado-nfs.py` uses a parameter file
  supplied on the command line, it does not automatically insert into the
  parameter set the options that are necessary for running jobs.
  Therefore, we need to add these options:

```
./cado-nfs.py --parameters ./parameters/factor/params.c60 90377629292003121684002147101760858109247336549001090677693 slaves.hostnames=localhost slaves.nrclients=2
```

* Run a full factorization on the local machine, using 8 threads for the
  server (this includes the linear algebra), and 4 jobs of 2 threads
  each for the polynomial selection and the sieving:

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 --slaves 4 --client-threads 2 --server-threads 8
```

* Run a factorization in the given directory (must be an absolute path),
  interrupt it (with Ctrl-C, or whatever unexpected event), and resume
  the computation:

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 workdir=/tmp/myfacto
[Ctrl-C]
./cado-nfs.py /tmp/myfacto/c60.parameters_snapshot.0
```

* Run a server on `machine1`, and a slave on `machine2`, disabling ssl:

```
machine1$ ./cado-nfs.py --server 90377629292003121684002147101760858109247336549001090677693 server.port=4242 server.ssl=no server.whitelist=machine2
machine2$ ./cado-nfs-client.py --server=http://machine1:4242 --bindir=...
```

Note: if you are on an insecure network, you'll have to activate ssl, and
then pass the appropriate sha1 certificate to the client (the server
prints the appropriate command-line to be copy-pasted on `machine2`).

* Run a factorization on `machine1`, and have it start automatically a
  slave on `machine2` via SSH:

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
slaves.hostnames=machine1,machine2
```

Note that, in that case, you have to specify `machine1` as well in the list
of `hostnames` if you want it to contribute to the polynomial selection and
the sieving.


Stopping the factorization during a specific step:
==================================================
It is possible to stop the factorization:

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 tasks.xxxx.run=false
```

This command works with: xxxx = polyselect, sieve, filter, linalg


Known problems:
===============

(some of these problems refer to versions or operating systems that are
no longer supported by cado-nfs anyway)

* when running the square root step in multi-thread mode (tasks.sqrt.threads=2
  or larger) with GMP <= 6.0, you might encounter an issue due to a "buglet"
  in GMP (<https://gmplib.org/list-archives/gmp-bugs/2015-March/003607.html>).
  Workaround: use tasks.sqrt.threads=1 or GMP >= 6.1.0.
* GCC 4.1.2 is known to miscompile CADO-NFS (see
  <https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/14490>),
  GCC 4.2.0, 4.2.1 and 4.2.2 are also affected.
* under NetBSD 5.1 amd64, Pthreads in the linear algebra step seem not to
  work, please use `-t 1` option in `cado-nfs.py` or `tasks.linalg.threads=1x1`.
* under AIX, if GMP is compiled in 64-bit mode, you should set the
  environment variable `OBJECT_MODE`, for example:
  export `OBJECT_MODE=64`
* if you encounter the error "Exceeded maximum number of failed workunits,
  maxfailed=100", you can restart the factorization with
  `tasks.maxfailed=200`.   But it would be wise, first, to try to
  understand _why_ workunits are failing. This should not appear. It
  might be that all your workunits are timing out, because your `adrange`
  and `qrange` parameters are too large.  Or it's a bug in cado-nfs, and
  then it should definitely be reported.


Contact, links:
===============

The website of the project is hosted at:
   <https://cado-nfs.inria.fr/>

You can get the latest development version with:
```
git clone https://gitlab.inria.fr/cado-nfs/cado-nfs.git
```
or
```
git clone git@gitlab.inria.fr:cado-nfs/cado-nfs.git
```
(use the latter if you have an account on Inria gitlab, and commit access
to cado-nfs)

There is now a unique mailing-list associated to Cado-nfs
[cado-nfs@inria.fr](https://sympa.inria.fr/sympa/info/cado-nfs). Please
do not use the old `cado-nfs-discuss` mailing list, the infrastructure
that hosts this mailing list has been removed in september 2021. All
mailing list content has been moved to the new mailing list, but links are
broken.

If you find a bug, if you have a problem compiling cado-nfs, if you want to
factor a large number and seek for advice for tuning the parameters, then
the cado-nfs list is the right place to ask.

On the <https://gitlab.inria.fr/cado-nfs/cado-nfs> web page you can also
find the cado-nfs bug tracker (a.k.a project issues). The bug tracker is
an important piece of the cado-nfs development cycle.  Submitting bugs
and merge requests there is welcome (you need an Inria gitlab account),
although if you are unsure, it might be better to speak up on the
cado-nfs mailing list first.
