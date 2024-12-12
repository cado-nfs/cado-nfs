
This set of scripts implement pretty-printing for internal cado-nfs types
within gdb

Not all the previously existing gdb pretty printing code has been
converted. The lacking bits are:
 - the `polymat` layer. It's still in the git history.
 - floats. Commented out in floats.py
 - fft transforms. Commented out in base.py

In order to gain this functionality, cado-nfs should be built with shared
libraries (and in debug mode, of course). The required scripts should be
auto-loaded by gdb. (it's also possible to do without shared libs, but
it's more annoying)
