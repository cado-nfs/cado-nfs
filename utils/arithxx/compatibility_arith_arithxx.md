# compatibility map between the arith and arithxx interfaces

This is as much about transition than about compatibility.

## the `mod_int*` interfaces

These concern the plain integers, not in any transformed form.

All of these are trivially handled by the `Integer_base` base class, or
by boilerplate ctors. Note that this also works with `cxx_mpz`'s. Names
are not exactly as below, though.
 - `init`, `clear`: -> ctor and dtor
 - `set`, `set_ul`, `set_uls`: -> ctors (copy, (`uint64_t`), (ptr+len))
 - `get_ul`, `get_uls`, `get_double`: -> explicit operators
 - `equal`, `equal_ul`, `cmp`, `cmp_ul`, `cmp_uint64`: -> comparisons
 - `fits_ul`: -> `fits_uint64_t`

These are fairly easy too.
 - `add`, `sub`: operators.
 - `mod`, `divexact`: `operator%` and `divexact.
 - `bits`: same name. This returns bit length - ctz, aka sizeinbase(2)


## the ```mod_*mod``` interfaces

A `Modulus` type holds an `Integer` type, and we of course have
ctor/dtor/conversion operations. Really nothing fancy here. Everything is
pretty much covered by the existing code.

 - `initmod_ul` `initmod_int`: -> ctor
 - `clearmod`: -> dtor
 - `getmod_ul` `getmod_int`: -> `getmod()`. There's no direct equivalent
   of `getmod_ul`, but it's probably not much of an issue.

## the residue interfaces.

### First, the various trivialities.
 - `init`, `init_noset0`, `clear`, `swap`: -> ctor/dtor
 - `set`, `set_ul`, `set_ul_reduced`, `set_int`, `set_int_reduced`: ->
   the `mod.set` and `mod.set_reduced` interfaces do this.
 - `get_ul`, `get_int`: -> `mod.get(r)` (takes a `Residue`, returns an
   `Integer`)
 - `set0`, `set1`: -> same name.
Note that arithxx includes several other means to create a residue, such
as `operator()`.

### Comparisons.

We don't have full-fledged comparisons of `Residue` types.
The existing interface only works via the `Modulus` member functions, and
only does the following comparisons. Under the hood, these could all be
done on the `Residue` types alone.
 - `equal`, `is0`, `is1`: -> `mod.{equal,is0,is1}`

### Addition/subtractions:
 - `add`, `add_ul`, `sub`, `sub_ul`, `neg`: -> `mod.{neg,add,sub}`
 - `add1`, `sub1`: -> `mod.{add1,sub1}`

### Mul/square.

These are of course speed critical, most (all?) interfaces
have these coded right in the header file.
 - `mul`, `sqr`: -> same name

### Divisions by constants.
 - `div2`, `div3`, `div5`, `div7`, `div11`, `div13`: -> same name

### Powers.
 - `pow_ul`, `2pow_ul`, `pow_mp`, `2pow_mp`: -> `pow`

### Primality tests:
 - `sprp`, `sprp2`, `isprime`: -> `is_strong_pseudoprime_base2` and
   `is_prime`. We also have `is_strong_lucas_pseudoprime`.

### Inversion / Euclidean algorithm
 - `inv`, `gcd`, `jacobi`: -> same name

### Iteration support

 - `next`, `finished`: the only use in our code base was not with the
   abstract interface but rather with modul_next and modul_finished
   directly.  The interface was clumsy isn't that of an iterator, and
   there's no test. -> **interface dropped**.

### batch inversion
 - `batchinv`: -> same name
 - `batch_Q_to_Fp` and `batch_Q_to_Fp_context`. These correspond to dead
   code in `las-fill-in-buckets.cpp`, meant to compute rational roots in
   batches. It got commented out at some point. The backend code is now
   tested again, but not put in production.

### I/O
 - `fprintf`, `printf`: -> all under the umbrella of libfmt overloads,
   now.

### Lucas sequences
 - `V`, `V_dadd`, `V_dbl`: TODO (Test)


## Interfaces that only exist in `arithxx`
 - `inv_odd` and `inv_powerof2`. These are in fact only specialized for
   the `arithxx_mod64` layer.
 - `intinv`. It's rather silly. Only specialized for `arithxx_modredc64`.



