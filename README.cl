(still in development)

Using CADO-NFS for class group computation of imaginary quadratic fields.
-------------------------------------------------------------------------

**** Basic usage

$ ./cado-nfs.py -qs -cl -13784472188539234591

It will output the class number and some information on the structure of the
class group (its exponent, its invariant factors and the generators)
corresponding to the given discriminant.

**** Notes

Discriminant should be negative, congruent to 1, 2 or 3 modulo 4 or divisible by
4 with quotient congruent to 2 or 3 modulo 4.
It implies that the discriminant need not be fundamental. The computed class
group corresponds to the class group of the order (not necessarily maximal) with
the given discriminant.

But the discriminant cannot have a square factor below the large prime bound.
It is check by the cado-nfs.py script and the computation aborts in this
case.

$ ./cado-nfs.py -qs -cl -51690857612769343560844414967703313862777343615087425622841
...
Critical:Check Discriminant: The discriminant -51690857612769343560844414967703313862777343615087425622841 has a square factor 874537^2 belonging to the factor base
...


Rationale: cado-nfs factors ideals into prime ideals of the maximal order. If
for an odd prime p, p^2 divides the discriminant, what we really want is the
factorization as prime ideals of the order of conductor p. If p does not belong
to the factor base, it does not matter because every ideals that could
appear in a factorization are coprime to p and in this case the factorization in
the maximal order and in the order of conductor p are equivalent. But if ideals
above p appear in the factor base, the fact that the computed factorizations are
between ideals of the maximal order means that what we really compute at the
end is the class group of the discriminant divided by p^2.

**** Importing your own polynomial file

Not supported.

It is assumed in lots of places that the polynomial is x^2-D

**** Creating your own parameter files

If the parameter file for your target size is missing, you can create them
by interpolating/extrapolating between existing parameter files. You need to
create a file params.qs.dNNN in the parameters/cl directory.

**** Known issues

Linear algebra is performed using CRT with 64-bit modulus. It implies that
p1 and/or pz should be present in BWC_GFP_ARITHMETIC_BACKENDS (e.g.,
BWC_GFP_ARITHMETIC_BACKENDS="p1;pz"). If only pz is present, linear algebra will
be slower.
Widths up to 15 words are compiled by default.
