# The arithxx layer

This modular arithmetic layer reimplements some stuff from
[../arith](../arith).

We have the following implementations.

 - [`mod64.hpp`](mod64.hpp) ; this is similar to the `mod_ul` layer from
   the old arithmetic layer, but here we're based on 64-bit `uint64_t`,
   not just `unsigned long`.

 - [`modredc64.hpp`](modredc64.hpp) ; same, with redc

 - [`modredc96.hpp`](modredc96.hpp) ; 1.5 words, with redc

 - [`modredc126.hpp`](modredc126.hpp) ; two 64-bit integers minus two bits, with redc

 - [`mod_mpz_new.hpp`](mod_mpz_new.hpp) ; an mpz-based type.

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


# compatibility map between the arith and arithxx interfaces

see [`compatibility_arith_arithxx.md`](compatibility_arith_arithxx.md)

# Internal structure of arithxx

The `mod_mpz_new` differs a lot from the rest, and it's best if you don't
look into it just yet. And in fact, this layer is so cumbersome that it
might be reimplemented differently at some point.

The primary citizens (layers) are `mod64`, `modredc64`, `modredc96`, and
`modredc126`.

All layers define three types: the `Integer` type, the `Residue` type,
and the `Modulus` type.

The `Integer` type is a fixed size integer, stored in a constant-size
array of 64-bit words, that is wide enough to store the modulus.

The `Residue` type is really just an `Integer` in disguise. We use it in
order to provide distinct overloads. Note that by design, in order to
deal with a Residue (e.g. if you want to construct such an object), you
must have a Modulus around.

The `Modulus` type is where all the stuff happens.

The `Modulus` type is implemented via the CRTP (Curiously recurring
template pattern) technique. A `Modulus` actually inherits from a
hierarchy of parent classes, parameterized by the layer class, and which
are all allowed to downcast to the modulus class. This is made so that
each of them can handle a specific subset of the functionalities,
possibly overriding default values.

In order, the possible base ("traits") classes are:
 - [`api<layer>`](arithxx_api.hpp). This implements the main modulus
   constructor, and in particular the `m` data field. Some set/get
   functions are implemented, as well as others that do things like
   divisions by constant integers. All of them have access to
   compile-time information about the modulus. *Default* implementations
   for these methods are in
   [`arithxx_api_impl.hpp`](arithxx_api_impl.hpp)
 - `api_bysize<layer`, which exists as two (partial) specializations
   [`arithxx_api_bysize<*, Integer64>`](arithxx_api64.hpp) and
   [`arithxx_api_bysize<*, Integer128>`](arithxx_api128.hpp). There's
   very little code for these classes. In a sense, most of what could go
   here could also, and perhaps better, go to the `Integer` classes, the
   `api` traits class, or maybe one of the `redc` traits classes if
   applicable. For the moment there is still a specific implementation of
   inversion for two
   machine words in [`arithxx_api128_impl.hpp`](arithxx_api128_impl.hpp),
   and that's it.
 - [`redc<layer>`](arithxx_redc.hpp). This is of course very important to
   redc implementations. This defines the specific data fields of moduli
   classes that use redc, in particular a `Residue` representation of one
   is stored. Some set/get functions are there, and they overload some of
   the ones that exist in `api<layer>`.
 - [`redc64<layer>`](arithxx_redc64.hpp) and
   [`redc128<layer>`](arithxx_redc128.hpp). These contain actual
   (possibly assembly, for instance) code that implement the redc
   routine, for a given layer. The key functions are `redc1` and
   `unredc1`, and overloads such as `frommontgomery` and `tomontgomery`
   are defined from these.

A final `Modulus` class may inherit from these, **starting from the last
one that is relevant**. E.g. it's probably either `api_bysize` or
`redc64` or `redc128`.


Note that several traits classes define at least some of their functions
with external (not inline) linkage, and the implementations are in the
`_impl.hpp` files. This is a design choice, given that we do not have to
expose all the code to all potential users.

For this reason, a final `Modulus` class has to include all the
`_impl.hpp` from its `.cpp` compilation unit, **AND** explicitly
instantiate **all the upper layers**, and this is so **even for layers
that have some of their functions implicitly instantiated by the layers
that derive from them**. This last bit is honestly a bit annoying, but
it's an obvious consequence of lazy template instantiation.

It may be called a nitpick, but this last bit causes some issues with the
`mod_mpz_new` layer, for which some of our "generic" code actually
doesn't apply. Inheritance actually masks all the non-functioning
templates, and they don't have to be instantiated.  However, since we
*do* instantiate the templates for the above reasons (because some
"generic" code in these layers such as `gcd` needs to be emitted), and
yet we don't provide explicit specializations for `mod_mpz_new`, we
have a problem. For this reason, `mod_mpz_new.cpp` actually *deletes*
these functions (which is a way of specializing).
