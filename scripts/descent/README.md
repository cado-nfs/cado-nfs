Examples of parameters can be found in the parameters/dlp directory

Algorithms
  - UpperClass:
      This is the initialization of the descent. An extended gcd is
      done between the target and p to get two "rational
      reconstructions", and then a sieving procedure is done to find a
      linear combination of these that is smooth (idea taken from
      Joux-Lercier). This sieving is done with las, with two linear
      polynomials.
      At the end, we have target = num/den, where num and den are
      smooth, with a smoothness bound that is large than the one that
      was used in the sieving / linear algebra.
  - MiddleClass:
      For all the primes dividing num and den that are larger than the
      large prime bound, a "special-q descent" is performed, in order
      to rewrite them in terms of smaller and smaller elements, until
      everything is known. This is done with the `las_descent` program.
  - LowerClass:
      This step is just putting everything together.
      In practice, this means computing appropriate Schirokauer maps
      and propagating the known logarithms in relations in order to
      deduce the unknown ones. The main tool for that is the
      reconstructlog program. Some ugly modifications of input files
      are necessary.


TODO: do we want to have default values, as done here, for the `--init-*`
arguments ? I would say no.

TODO: of course this is currently kludgy, and does not properly wield
the power of the cadofactor python programs. I need to understand that
stuff better.

TODO: keep las awake, if needed.

TODO: make output less verbose.
