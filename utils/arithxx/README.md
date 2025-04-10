This modular arithmetic layer reimplements some stuff from
[../arith](../arith).

We have the following implementations.

 - [`mod64.hpp`](mod64.hpp) ; this is similar to the `mod_ul` layer from
   the old arithmetic layer, but here we're based on 64-bit `uint64_t`,
   not just `unsigned long`.

 - [`modredc64.hpp`](modredc64.hpp) ; same, with redc

 - [`modredc126.hpp`](modredc126.hpp) ; two 64-bit integers minus two bits.

 - [`mod_mpz_new.hpp`](mod_mpz_new.hpp) ; an mpz-based type.

Note that **we do not have code for 96-bit moduli at this point**.

Implementations are both inlines in the `.hpp` files, and code in the
`.cpp` files.

Some of the code is shared between the implementations, and therefore we
can write it generically. More precisely,
[`arithxxx_api_impl.hpp`](arithxxx_api_impl.hpp) and
[`arithxxx_api64_impl.hpp`](arithxxx_api64_impl.hpp) contain code that is
common to all (or at least several) implementations.

We have some lower level plumbing code in [`u64arith.h`](u64arith.h)
as well as [`modint.hpp`](modint.hpp).

Finally, [`mod_stdop.hpp`](mod_stdop.hpp) defines the `ResidueOp`
template class, which can be used to get nice operator overloads. Note
that this does not provide copy-less compound operations.


TODO:
 - We need a 1.5 word layer, so that we can eventually get rid of the C
   interface.
