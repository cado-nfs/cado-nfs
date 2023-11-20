The `cado-nfs.py` Python script is generally run with

```
./cado-nfs.py parameterfile
```

All the parameters for the factorization are read from the parameter file,
but it is possible to specify such parameters on the command line, after
parameterfile.

For example, running

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 tasks.workdir=/tmp/c59 server.whitelist=0.0.0.0/0 --server

./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 tasks.workdir=/tmp/c59 tasks.execpath=$HOME/build/cado-nfs/normal server.whitelist=0.0.0.0/0 --server
```

starts the `cado-nfs.py` script, which also starts the server. It does not
start any clients with this command line, so those would have to be
started manually. This is caused by the `--server flag`, which enforces the
"do not start clients" feature. The version above with `tasks.execpath`
gives you control over the place where the binaries will be found.

The `./cado-nfs.py` script, when not instructed to start clients
explicitly (that is, with `slaves.hostnames` not set), decides
nevertheless to do so when it is using a parameter file that was picked
automatically from the set of default parameter files that are shipped
with cado-nfs. If provided with an explicit parameter file, the default
is to obey what is in there, and therefore not start clients unless
`slaves.hostnames` is set.

To start additional clients, do (following the command line given in the
output of the server):

```
./cado-nfs-client.py --server=https://quiche.loria.fr:8001 --certsha1=[Certificate SHA1] --bindir=...
```

where the `--server` and `--certsha1` parameters should be given the URL and
certificate SHA1 value of the server, respectively, as printed by the
`cado-nfs.py` script. You can start an arbitrary number of client scripts,
on any machines that can connect to the server via HTTP.

The (optional) `--bindir` parameter tells where to find the binaries to use
instead of having them downloaded from the server. It is possible to pass
the root of the build tree as `--bindir`, for example with

```
--bindir $(eval `make show` ; echo $build_tree)
```

If you want to let the server automatically start clients, you need to supply a list of hostnames on which to start clients, using constructs similar to the following two examples.

```
./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 tasks.workdir=/tmp/c59 slaves.nrclients=2 slaves.hostnames=localhost

./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693 tasks.workdir=/tmp/c59 slaves.nrclients=2 slaves.hostnames=localhost tasks.execpath=$HOME/build/cado-nfs/normal slaves.scriptpath=$HOME/git/cado-nfs/
```

to let it start two clients on localhost. The second version gives you
more control power, but should hardly be needed.  The `slaves.scriptpath`
parameter must be the path to the directory on the client machine which
contains `cado-nfs-client.py` and `workunit.py`. The latter can also be
fetched from `scriptpath+"/scripts/cadofactor"`, whence giving the path to
the source tree as `scriptpath` should be fine.


For complex set-ups, it is preferable to write a parameter file. Some
examples are in
[`scripts/cadofactor/parameters`](scripts/cadofactor/parameters),
[`scripts/cadofactor/parameters.oar`](scripts/cadofactor/parameters.oar),
[`scripts/cadofactor/parameters.rsa512.oar`](scripts/cadofactor/parameters.rsa512.oar),
The [`parameters/factor/params.c90`](parameters/factor/params.c90) file
contains extensive comments describing the most common parameters.

The `.oar` parameter files are meant for `cado-nfs.py` scripts that run
*inside* an OAR submission on clusters that use OAR as the job scheduler,
such as Grid5000, as they read the list of slave hostnames from the OAR
node file.

For example,

```
oarsub -I
./cado-nfs.py parameters.oar
```

factors the usual c59 test number using the nodes reserved via OAR, in
this case one node.  The
[`scripts/cadofactor/parameters.oar`](scripts/cadofactor/parameters.oar),
file contains the line `slaves.catrel.nrclients = 8` which tells the
script to launch 8 clients on each unique host name (=node); the
parameter `threads=2` causes all the programs to use two threads,
resulting in 16 threads per node.

If a factorization is interrupted, it is usually possible to resume it where
it left off, simply by running `cado-nfs.py` again with the same command
line parameters, or the latest parameter file snapshot from the work
directory (found at `[[work directory]]/XXXX.parameters_snapshot.YYY`).
These command line parameters are printed to the screen and
written to the `.log` file in the factorization's working directory on each
run.

If a factorization cannot be resumed for whatever reason, it is often
possible to salvage some of the output files which can then be imported into
a new factorization run of the same number, but with a new working
directory.

Running with several groups of machines
=======================================

You can define several groups of machines slaves.xxx, slaves.yyy, ...
Here is an example:

```
slaves.big.hostnames = fondue,berthoud
slaves.big.nrclients = 4
slaves.big.scriptpath = /users/caramel/zimmerma/svn/cado-nfs
slaves.big.basepath = /tmp/nfs
slaves.small.hostnames = sel,poivre
slaves.small.nrclients = 2
slaves.small.scriptpath = /users/caramel/zimmerma/svn/cado-nfs
slaves.small.basepath = /tmp/nfs
```

This will run 4 clients (with by default 2 threads each) on the "big" machines
fondue and berthoud, and 2 clients (with by default 2 threads each) on the
"small" machines sel and poivre.

Importing a polynomial file
===========================

If you want to import a polynomial file (either an SNFS polynomial, or
if polynomial selection was already completed in an earlier run), use:

```
tasks.polyselect.import=xxx.poly
```

where `xxx.poly` is the name of the polynomial file (in CADO-NFS format).
You have to make sure that imported polynomials are irreducible and in
particular have content 1.

This imports the polynomial and prevents any additional polynomial selection
from running, i.e., the imported polynomial is used unconditionally.


The polynomial selection run by `cado-nfs.py` is performed in two phases:
the first phase searches for polynomials with good size property and
keeps in a priority queue the `nrkeep` most promising ones, then the
second phase root-optimizes these to find the overall best polynomial by
Murphy E value.

You can import files that were previously produced by phase 1 (resp.
phase 2) into phase 1 (resp. phase 2) again; the imported files will be
processed as if phase 1 (resp. phase 2) had generated them by itself.
Phase 1 will add any polynomials from the imported files into its
priority queue and forward the kept ones to phase 2, while phase 2 will
choose the best polynomial by Murphy E from the imported ones and the
ones it computes by itself.


If you want to import size-optimized polynomials into phase 1 of polynomial
selection and then continue polynomial search, use

```
tasks.polyselect.import_sopt=file
```

or

```
tasks.polyselect.import_sopt=@file
```

With `@file`, the file `file` should contain a list of files to import,
one file name per line. All polynomials from the given file(s) are
imported and added to the phase 1 priority queue, keeping the nrkeep best
ones. Polynomial selection phase 1 then continues normally. For example,
if an earlier run completed ad values up to `ad=200000` and you want to
resume the search from there, you could use

```
tasks.polyselect.import_sopt=@list_of_existing_files tasks.polyselect.admin=200000
```

which imports the existing files, then resumes searching at `ad=200000`
up to the `admax` value given in the parameter file.


If you want to import root-optimized polynomials into phase 2 of
polynomial selection and then continue polynomial search, use

```
tasks.polyselect.import_ropt=file
```

or

```
tasks.polyselect.import_ropt=@file
```

This reads the polynomial(s) in the given file(s), then root-optimizes
any polynomials found in phase 1; the best polynomial (rated by the
Murphy E value) will be used for the sieving.

Warning: if a polynomial file does not specify a Murphy E value and is
imported into phase 2, its Murphy E value is set to 0 by default. Since
any polynomials found by the polynomial search have positive Murphy E,
the imported one will always "lose". To import, e.g., an SNFS polynomial
without Murphy E, use `tasks.polyselect.import`, and not `import_ropt`.

Sieving with composite special-q's
==================================

By default the sieving only uses prime special-q's. To allow composite
special-q's you have to set `tasks.sieve.allow_compsq` to true, for example:

```
tasks.sieve.allow_compsq = true
tasks.sieve.qfac_min = 50
tasks.sieve.qfac_max = 100000
```

This will allow composite special-q's with smallest prime factor at least
50, and largest prime factor at most 100000.

If `tasks.sieve.qfac_max` is not set, it is considered as infinite, which
means that prime special-q's will always be used (in addition to
composite s-q's).

If sieving over [q0,q1], with `tasks.sieve.qfac_max < q0`, then all
special-q's will be composite.

Sieving with large special-q's
==============================

If you get an error like the following:

```
   The special q (27 bits) is larger than the large prime bound on side 1 (26 bits).
   You can disable this check with the -allow-largesq argument,
   It is for instance useful for the descent.
```

then you can solve it as follows:

```
tasks.sieve.allow_largesq = true
```

Importing relations
===================

Important: when importing relations, it is assumed one also imports the
corresponding polynomial (a few relations from each file will be checked
against that polynomial). If no polynomial is imported, anything can happen.

If you want to import already computed relations, use:

```
tasks.sieve.import=foo
```

where `foo` is the name of the relation file to be imported (in CADO-NFS
format). Alternatively you can do:

```
tasks.sieve.import=@foo
```

where `foo` is a file containing several names of relation files, one file
name per line.

Adjust the `tasks.qmin` parameter accordingly to prevent
already-sieved special-q ranges from being sieved again. To do that:

* go into the `cxx.upload/` directory
* list the `cxx.xxx-yyy.*.gz` files
* find the largest `yyy` value, and use it as new `qmin`

Be aware that importing relation files that were not produced by
CADO-NFS's siever (for instance, handcrafted relations) is not supported:
the script will assume various things, like the presence of certain
commented lines. The list of these particularities may vary over time,
so we don't even try to list them here.

When importing relations, some statistics are gathered by parsing the
relation files. This is a costly operation, especially when relations
from a large number of files are gathered. An expedient for this problem
is to create, for each file `cxx.xxx-yyy.zzzzzzz.gz`, a file
`cxx.xxx-yyy.zzzzzzz.gz.stats.txt`, which contain the relevant
statistics. The format of this file should match the format of the
statistics files produced by the cado-nfs siever with the `-stats-stderr`
option. Examples of such statistics files are found as the files
`*.stderr0` in the `$name.upload/` directory for a cado-nfs factorization.
The only crucial bit information is the line which gives the number of
relations in the file, which has to match the following format:

```
# Total $nrels reports [%.2fs/r, %.2fr/sq]
```

Here is an example of how to produce such files:

```
#!/usr/bin/env bash
find $workdir/*.upload/ -name '*.gz' | while read f ; do
   zcat $f | tail -n 1 > $f.stats.txt
done
```

File locking when using sqlite3 as a database backend
=====================================================

The SqLite database used by `cado-nfs.py` by default makes extensive use of
file locking.  This requires that the underlying file system properly
implements POSIX file locking. One example of a file system that claims
to, but does not, is glusterfs, which leads to random SqLite errors. See
[this thread](http://gluster.org/pipermail/gluster-users/2011-September/031720.html)

If you encounter `Database corrupted` or similar errors, try if storing the
working directory on a local (not network-mounted) file system resolves the
problem.

If you have no choice but to use a filesystem which does not handle
locking correctly, you may specify the database location to an external
place (see next)

Storing the database elsewhere, or using a different database back-end
======================================================================

The database which that the computation state can be backed by a mysql
server. This is optional, and requires the `python3-mysql.connector`
package (on Debian linux -- actual package name depends on your OS
distribution).
  
To activate this feature, pass the database URI as an extra parameter to
`cado-nfs.py`, e.g.  `database=db:mysql://USER:PASS@host:port/foobar`,
where `port`, `PASS`, and `USER` are optional. The name `foobar` is how
you choose to name your database inside the database server.

The same syntax can be used to force the use of an externally located
sqlite3 database, e.g. with `database=db:sqlite3:///home/john/rsa1024.db`

Controlling the filtering stage
===============================

By default CADO-NFS requires an initial excess (after the first singleton
removal step) of 0% more relations than ideals. This is controlled by:

```
tasks.filter.required_excess=0.0
```

If you want a larger excess, say 20%, just add on the cado-nfs.py command
line (or in the parameter file):

```
tasks.filter.required_excess=0.2
```

It also requires that the initial excess is at least some given value. If you
want to modify it:

```
tasks.filter.purge.keep=1000
```

By default CADO-NFS discards all the excess (larger than `keep`) during the
purge step. If you want to keep some excess for the merge step:

$ ./cado-nfs.py ... `tasks.filter.purge.keep=1000`

If your factorization already started the linear algebra step, and you
want to do more sieving, you can restart it with a larger `rels_wanted` than
the current number of relations. For example if you have say 1000000 relations
(grep "is now" in the log file), just add in the cado-nfs.py line:

```
tasks.sieve.rels_wanted=1000001
```

If on the contrary you want to start filtering as soon as you have a positive
excess (or you know you have enough relations in case of imported relations), use:

```
tasks.sieve.rels_wanted=1
```

If not enough relations are collected after the filtering step, the sieve
is executed again with an additional number of wanted relations.  The
parameter `tasks.filter.add_ratio` controls the number of additional
required relations as a ratio of the number of unique relations already
collected:

```
tasks.filter.add_ratio=0.05
```

specifies that you want 5% more relations before next filtering step.
The default value is 0.01 (i.e. 1% more relations).

Printing and manipulating the database
======================================

The `wudb.py` script can be used to print info on the workunits. For example,

```
./wudb.py -dbfile /tmp/work/testrun.db -dump -assigned
```

prints all currently assigned work units.

To cancel a workunit:

```
./wudb.py -dbfile /tmp/work/testrun.db -cancel -wuid wuname
```

where `wuname` appears after `Workunit` in the `-dump -assigned` call.

This script can also be used to manipulate the database, for example:

```
./wudb.py -dbfile /tmp/work/testrun.db -setdict sqrt next_dep int 0
```

sets the `next_dep` variable in the `SqrtTask`'s state to the integer value 0.

Using several threads for the square root step
==============================================

By default the given number of threads (`-t` option of `cado-nfs.py`, or
`tasks.threads`) is used in the square root step. To use say
16 threads in parallel, add the following to the `cado-nfs.py` command line:

```
tasks.sqrt.threads = 16
```

Re-running the linear algebra step
==================================

If the linear algebra step was already started, and you want to run it
again with different parameters or with a different set of relations (for
example with more sieving), or just because you'd like to run it again
from scratch, then you need to add the following to the cado-nfs.py
command line, in order to force the previous linear algebra subdirectory
to be removed:

```
tasks.linalg.force_wipeout = True
```

Debugging
=========

For debugging the `cado-nfs.py` script, you might use:

```
$ ./cado-nfs.py ... --filelog DEBUG
```

and this will print more information in the `cxxx.log` file.

