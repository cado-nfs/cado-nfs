This is the *old* modular arithmetic layer of cado-nfs

It's mostly based on macro renaming and typedefs. Newer code should use
the files in arithxx instead.

We have five distinct implementations.

  - [`mod_mpz`](mod_mpz.h)         -> based on mpz
  - [`mod_ul`](mod_ul)             -> based on unsigned longs
  - [`modredc_ul`](modredc_ul)     -> based on unsigned longs + redc
  - [`modredc_15ul`](modredc_15ul) -> based on two unsigned longs, with
                                      products of less than three words 
  - [`modredc_2ul2`](modredc_2ul2) -> based on two unsigned longs.

The overall structure is as follows. In the descriptions below, X stands
for one of `_mpz`, `_ul`, `redc_ul`, `redc_15ul`, `redc_2ul2`

 - `modX.h` defines types and inline functions. All the `modX` files can
   live one along the others. Implementations of non-inline functions are
   in `modX.c`. A `modX.h` file must of course adhere to several
   constraints, one of them being the necessity to define some specific
   macros that are later used for renaming.

 - `modX_default.h` is a proxy header, which uses macros in
   [`mod_rename.h`](mod_rename.h) in order to expose a coherent set of
   functions. These are only macros!

   Note in particular that a function with C linkage that uses these
   types, when duplicated in a library, should have its name multiplexed
   in the same way. Even if it is a function with C++ linkage, operator
   overloading doesn't save us because the underlying types of our
   different implementations are likely to collide (e.g. `mod_ul` and
   `modredc_ul` are both based on 1-size arrays of `unsigned long`s).

 - We have some lower level plumbing code in [`ularith.h`](ularith.h). 

 - Some of the code can be shared between the implementations, and
   therefore we can write it generically. More precisely:
    - [`mod_common.c`](mod_common.c) contains code that is common to all
      implementations.
    - [`mod_ul_common.c`](mod_ul_common.c) contains code that works for
      both `mod_ul` and `modredc_ul`.
    - [`modredc_2ul_common.c`](modredc_2ul_common.c) contains code that
      works for 
      both `modredc_15ul` and `modredc_2ul2` (both have same width data
      types).

 - [`modul_poly.h`](modul_poly.h) and [`modul_poly.c`](modul_poly.c) are
   an implementation of a generic polynomial layer based on the types
   above.
