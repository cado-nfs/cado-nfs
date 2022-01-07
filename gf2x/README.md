
The `gf2x` software library
===========================

See the [releases](../../releases) pages for the different `gf2x`
releases.

`gf2x` is a C/C++ software package containing
routines for fast arithmetic in GF(2)[x] (multiplication, squaring,
GCD) and searching for irreducible/primitive trinomials.

Authors: Richard Brent, Pierrick Gaudry, Emmanuel Thomé, Paul Zimmermann.

License (for the library) is either:
 - If the archive contains a file named `toom-gpl.c` (not a trivial
   placeholder), the GNU General Public License, either version 3 of the
   License, or (at your option) any later version.
 - If the archive contains a file named `toom-gpl.c` which is a trivial
   placeholder, the GNU Lesser General Public License, either version 2.1
   of the License, or (at your option) any later version.

License for the [`apps/`](apps/) subdir is the GNU General Public
License, either version 2 of the License, or (at your option) any later
version.


Dependencies
============

`gf2x` has no external dependencies.

Some of the demos in the [`apps/`](apps/) subdirectory require the gmp
and NTL libraries. (use `./configure ; make` from that directory. You may
want to use `PKG_CONFIG_PATH` to have the autotools stuff there find
`gf2x` properly).


WARNING: The text below is somewhat outdated, and was originally written
for `gf2x-1.0` at best, and was marginally updated since. It is still
"reasonably accurate", but must be followed with caution.

Summary:
========

This README covers:

 - Package contents
 - Caution for gcc users
 - Instructions to install the package
 - Caution regarding installation
 - Hooking `gf2x` into ntl
 - Using the library
 - Contacting developers


Package contents:
=================

It contains the following files:

Miscellaneous doc files:
- [`README`](README)
- [`BUGS`](BUGS)
- [`src/TODO`](src/TODO)
- [`src/README`](src/README)
- [`already_tuned/tuned/README`](already_tuned/tuned/README)
- [`AUTHORS`](AUTHORS)
- [`ChangeLog`](ChangeLog)

Actual code:
- [`gf2x.h`](gf2x.h)            ; main api
- [`gf2x-impl.h`](gf2x-impl.h)  ; internal api
- [`gf2x.c`](gf2x.c)            ; top-level source for multiplication code
- [`gf2x-small.h`](gf2x-small.h); small-sized inlined multiplication routines
- [`toom.c`](toom.c)            ; main file for Karatsuba and Toom-Cook multiplication
- [`toom-gpl.c`](toom-gpl.c)    ; same, for GPL-tainted distribution
- [`fft.c`](fft.c)              ; multiplication using Fast Fourier Transform

Code adapted for the selected hardware
- [`already_tuned/`](already_tuned/)              ; pre-configured codes for selected architectures.
- [`already_tuned/*/gf2x-thresholds.h`](already_tuned/) ; pre-tuned thresholds files
- [`already_tuned/generic64/`](already_tuned/generic64/)    ; code that works on any 64-bit platform
- [`already_tuned/generic64/gf2x-thresholds.h`](already_tuned/generic64/gf2x-thresholds.h) ; placeholder thresholds
- [`already_tuned/generic32/`](already_tuned/generic32/)    ; code that works on any 32-bit platform
- [`already_tuned/generic32/gf2x-thresholds.h`](already_tuned/generic32/gf2x-thresholds.h) ; placeholder thresholds
- [`already_tuned/x86_64/`](already_tuned/x86_64/)       ; code that works on amd64 and intel core2
- [`already_tuned/x86_sse2/`](already_tuned/x86_sse2/)     ; code that works on x86 platforms supporting sse2
- [`already_tuned/generic/`](already_tuned/generic/)      ; code that works everywhere. Does _not_ include mul1
- [`gf2x/`](gf2x/)                  ; place where symlinks to the files above go

For testing:
- [`tests/check-mul.c`](tests/check-mul.c)       ; simple check program
- [`tests/do-check-mul.sh`](tests/do-check-mul.sh)   ; shell script driving check-mul

For tuning:
- [`lowlevel/mul*.c`](lowlevel/) ; various candidate code samples for basic routines
- [`src/tuneup.c`](src/tuneup.c) ; tuning program for basecase multiplication
- [`src/tunetoom.c`](src/tunetoom.c) ; tuning program for Karatsuba/Toom-Cook multiplication
- [`src/tunefft.c`](src/tunefft.c) ; tuning program for FFT multiplication
src/tune-lowlevel.pl
- [`src/gen_bb_mul_code.c`](src/gen_bb_mul_code.c)   ; program to generate many alternatives for mul1
- [`src/replace.h`](src/replace.h) ; helper code
- [`src/replace.c`](src/replace.c) ; helper code
- [`src/tuning-common.h`](src/tuning-common.h) ; helper code
- [`src/tuning-common.c`](src/tuning-common.c) ; helper code
- [`src/timing.h`](src/timing.h) ; helper code
- [`src/timing.c`](src/timing.c) ; helper code
- [`src/modify-thresholds.c`](src/modify-thresholds.c) ; helper code

Applications that use gf2x and NTL. These applications are covered by the
GPL.
- [`apps/halfgcd.hpp`](apps/halfgcd.hpp) ; subquadratic gcd over GF(2)[x]
- [`apps/halfgcd.cpp`](apps/halfgcd.cpp) ; subquadratic gcd over GF(2)[x]
- [`apps/factor.cpp`](apps/factor.cpp)  ; finds smallest irreducible factor of trinomial over GF(2)
- [`apps/check*.sh`](apps/)   ; some tests using factor

Caution for users of old gcc versions
=====================================

gcc versions 4.3.0 and 4.3.1 have a bug which affects gf2x in a an
unpredictable way. It is recommended to upgrade to at least 4.3.2, or
configure with `--disable-sse2`

Instructions to install the package:
====================================

Cautious users follow steps 1 to 5 below. Urged users follow only 1 and 5.

1) Type:
    ```
      ./configure && make
    ```

   See also the notes for specific systems below, regarding ABI selection
   or instruction set selection.

   Several flags and arguments can be passed to `./configure` (or set
   as environment variables) ; see `./configure --help`

2) Highly recommended ; run
    ```
      make check
    ```
   To guard against possible bugs (either from `gf2x` or from the compiler).

3) Optional, but recommended: tune low-level routines, Karatsuba/Toom-Cook and
   FFT multiplication:
    ```
      make tune-lowlevel
      make tune-toom
      make tune-fft
    ```
   Note that there are some minor issues related to tuning on amd64 -- see BUGS.

   Tuning unfortunately takes a long while ; it is possible to decrease
   the time spent in tuning using additional parameters to the `make
   tune-*` commands, more specifically: ``` make tune-toom
   TOOM_TUNING_LIMIT=<max size in words> make tune-fft
   FFT_TUNING_LIMIT=<max size in words> ``` The default values for
   `TOOM_TUNING_LIMIT` and `FFT_TUNING_LIMIT` are `2048` and `8000000`,
   respectively. `TOOM_TUNING_LIMIT` must be in the range `[30,2048]`,
   while `FFT_TUNING_LIMIT` is only limited by the available memory.

   More detailed information on tuning can be found in the
   [`src/tunetoom.c`](src/tunetoom.c) and [`src/tunefft.c`](src/tunefft.c)
   files, as well as [`src/README`](src/README).

   All the tuning targets automatically rebuild the library in order to
   incorporate the tuned parameters. So an additional `make` is
   unnecessary (although innocuous, since the library should be up to
   date).

4) Highly recommended ; run
    ```
      make check
    ```
   to see whether everything went ok.

5) Either run:
    ```
      make install
    ```
   Or, if you prefer, keep your build tree, and have your programs
   include `<buildtree>/gf2x.h`, and link against the library
   `<buildtree>/.libs/libgf2x.a` (static) or
   `<buildtree>/.libs/libgf2x.so` (dynamic). Usual pitfalls apply
   regarding dynamic linking: if you are familiar with none of the
   `-Wl,-rpath,<path>` or `LD_LIBRARY_PATH` ways of making this work,
   your easiest bet is to stick to the static library approach.

Notes for specific systems:
===========================

ABI selection
-------------

Your hardware platform may support several ABIs (Application Binary
Interfaces), corresponding to the type `unsigned long` being either
32-bit or 64-bit wide. The `gf2x` package accomodates for this under the
assumption that ABI selection is covered by the selection of the
appropriate compiler options. In order to compile for an ABI different
from the default one, you have to pass additional parameters to the
configure script:
```
  ./configure ABI=<bit width of unsigned long> CFLAGS=<corresponding CFLAGS> && make
```
For example on a Mac OS X computer with an Intel Core 2 processor using
gcc, one may use: `./configure ABI=64 CFLAGS="-O2 -m64" && make`

(equivalently, one may also use `ABI=64 CC='gcc -m64'`)

Note that ABI selection is sometimes tricky. Be sure that you select the
proper combination of ABI= and CFLAGS= parameters. You must also make
sure that those correspond to the CFLAGS that were used by any other
binary object you're linking `gf2x` with. That applies, in particular, to
the GMP library if you intend to compile the files in [`apps/`](apps/)

Instruction set selection
-------------------------

In order to select a specific instruction set, you should pass a `-march`
flag to the compiler (this implicitly means that we are using a compiler
that understands `-march`). For example something like `-march=x86-64` to
compile for generic x86-64 hardware.

When you do that, gf2x will **NOT** try to enable instruction sets that
are outside the instruction set that you selected. Otherwise, gf2x's
default behaviour is to go out of its way to enable all possible
instruction sets.

Note also that in some rare circumstances (some virtual machine setups),
it may be relevant to set `-march=native` explicitly. This will trust the
compiler's feature detection as the ultimate authority, and not try to
guess what runs and what does not run.

Other
-----

Under AIX the following might be required to build the `gf2x` binaries:

```
  $ export OBJECT_MODE=64
```

Caution regarding installation:
===============================

The header files installed in `$includedir/gf2x/` are architecure-dependent.

Hooking gf2x into ntl:
======================

`ntl-5.5` offers the possibility to replace the multiplication of
polynomials with the `gf2x` code. For this, you have to install the
`gf2x` library in your favorite `/path/to/gf2x/` (e.g. using `gf2x`'s
`./configure` --prefix=/path/to/gf2x/ && make && make install). Then
configure NTL with
```
./configure NTL_GF2X_LIB=on GF2X_PREFIX=/path/to/gf2x/
```

If for some reason `gf2x` was installed with custom directories, you may
also consider the `GF2X_INCDIR` and `GF2X_LIBDIR` options ; some setups will
typically require `GF2X_LIBDIR=/path/to/gf2x/lib64`, for example.

Using the library:
==================

`gf2x` exports one public header called [`gf2x.h`](gf2x.h) . This header
exports one function `gf2x_mul`, whose use is documented in
[`gf2x.h`](gf2x.h). Unlike with versions of `gf2x` up to 1.1, the
function `gf2x_mul` is now reentrant as of `gf2x-1.2`, and manages its
own allocation needs. The companion function `gf2x_mul_r` is here if you
want to control the allocation needs -- within a call to `gf2x_mul_r`,
below the FFT threshold, no allocation is performed. Above, there is
some, and this is an acknowledged bug.

Some "half-public"' headers are in the [`gf2x/`](gf2x/) subdirectory.
Some of these headers (`gf2x/gf2x-thresholds.h` and `gf2x/gf2x_mul*`) are
architecture-dependent. In order to use the inlined multiplication
routines for small sizes, you may use
[`gf2x/gf2x-small.h`](gf2x/gf2x-small.h).

By default, `gf2x` creates both static and dynamic libraries.

Contacting developers
=====================

The `gf2x` discussion mailing list can be reached at `gf2x@inria.fr` ;
the list info page is [here](https://sympa.inria.fr/sympa/info/gf2x)

Note that traffic on that mailing list is very small.

It is also possible to post issues to the gitlab site, but please [read
this
first](https://sympa.inria.fr/sympa/arc/cado-nfs/2020-10/msg00006.html).

(Very) Old versions
===================

The gzipped tar file for `gf2x` version `0.3.1` is
[here](http://maths-people.anu.edu.au/~brent/ftp/trinom/gf2x-0.3.1.tar.gz).

The gzipped tar file for `gf2x` version `0.2` is
[here](http://maths-people.anu.edu.au/~brent/ftp/trinom/gf2x-0.2.tar.gz).

The gzipped tar file for `gf2x` version `0.1` is
[here](http://maths-people.anu.edu.au/~brent/ftp/trinom/gf2x-0.1.tar.gz).

References
==========

[Richard P. Brent](http://maths-people.anu.edu.au/~brent) and
[Paul Zimmermann](https://members.loria.fr/PZimmermann/)
[A multi-level blocking distinct degree factorization
algorithm](http://maths-people.anu.edu.au/~brent/pub/pub230.html),
presented at the *Eighth International Conference on Finite Fields and
Applications (Fq8)*, Melbourne, 9-13 July 2007. Published in
*Contemporary Mathematics*, Vol. 461, 2008, 47-58.
Also appeared as [INRIA Tech Report RR-6331](http://hal.inria.fr/inria-00181029/) (October 2007, 16 pp.) and as [arXiv:0710.4410](http://arxiv.org/abs/0710.4410).

[Richard P. Brent](http://maths-people.anu.edu.au/~brent) and
[Pierrick Gaudry](https://members.loria.fr/PGaudry/)
[Emmanuel Thomé](https://members.loria.fr/EThome/)
[Paul Zimmermann](https://members.loria.fr/PZimmermann/)
[Faster multiplication in
GF(2)[x]](http://maths-people.anu.edu.au/~brent/pub/pub232.html)
*Proc. ANTS-VIII* (Banff, May 2008),
*Lecture Notes in Computer Science*, Vol. 5011,
Springer-Verlag, 2008, 153-166.
Also [INRIA Tech. Report RR-6359](http://hal.inria.fr/inria-00188261/) (Nov. 2007, 19 pp.)

[Richard P. Brent](http://maths-people.anu.edu.au/~brent) and
[Paul Zimmermann](https://members.loria.fr/PZimmermann/)
[Ten new primitive binary trinomials](http://maths-people.anu.edu.au/~brent/pub/pub233.html), *Mathematics of Computation*, Volume 78, Number 266, April 2009, Pages 1197–1199.
[article on Inria HAL](https://hal.inria.fr/inria-00337525)
