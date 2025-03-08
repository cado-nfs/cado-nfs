import functools

from sage.arith.functions import lcm
from sage.arith.misc import gcd, valuation
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.misc.flatten import flatten
from sage.misc.functional import sqrt
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.complex_double import CDF
from sage.rings.complex_mpfr import ComplexField
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Zp
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RealField
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.symbolic.ring import SR
from sage.functions.log import exp
from cypari2.pari_instance import Pari

# This class is only here to embed code that is attached to the cado-nfs
# diagram. Of course, almost everything we'll ever need is there in
# sagemath, but we still do miss a few things, or want to work around a
# few bugs.


class CadoNumberFieldWrapper(object):
    def __init__(self, K, trial_division_bound=10**7):
        # Many operations on ideals will internally call K.pari_nf(), which factors the discriminant.
        # This will stall when the discriminant is huge.
        # pari does let you give a trial division bound for factoring the discriminant, but sage doesn't properly expose it.
        # So we initialize K._pari_nf ourselves. The code below is K.pari_nf() except we provide a trial division bound.
        pari = Pari()
        nf = pari([K.pari_polynomial("y"), trial_division_bound])
        K._pari_nf = nf.nfinit()
        self.K = K

    def LogMap(self, *args, **kwargs):
        """
        Compute the log embeddings.

        The embeddings are returned so that the exponential of the sum of
        the log embeddings is the algebraic norm

        In essence, the code below is exactly the same as
        NumberField.log_embeddings from sagemath, except that we want to
        fix a sagemath bug, which is that the precision is ignored.

        https://github.com/sagemath/sage/blob/79c047c0a22a98bea4567d182c694fd4df1aea81/src/sage/rings/number_field/number_field.py#L9574

        further fixes:
         - clean up the code,
         - add better treatment of the precision (in the same way as we
           corrected K.places())
         - add better treatment of log(0)

        TESTS:

            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^3 - 2);
            sage: from cado_sage import CadoNumberFieldWrapper
            sage: N = CadoNumberFieldWrapper(F)
            sage: z = alpha^2 - 12 * alpha + 47
            sage: L = N.LogMap()
            sage: ratio = exp(sum(L(z))) / z.norm()
            sage: ratio > 0.99 and ratio < 1.01
            True

        """
        pl = self.places(*args, **kwargs)
        K = self.K
        r1 = len(K.defining_polynomial().real_roots())
        r2 = (len(K.defining_polynomial().complex_roots()) - r1) // 2
        R = getattr(pl[0].codomain(), 'real_field', pl[0].codomain())
        def closure_map(x):
            return vector([R(abs(sigma(x))).log() for sigma in pl[:r1]] +
                          [2 * R(abs(tau(x))).log() for tau in pl[r1:]])
        return Hom(K, VectorSpace(R, r1 + r2), Sets())(closure_map)


    # def LogDriftMap(self, prec=53):
    #     """
    #     This is almost like the log embeddings, except that we're
    #     actually returning the drift with respect to the norm.
    #
    #     At some point we thought it was useful.
    #     """
    #
    #     L = self.LogMap(prec)
    #
    #     def c(x):
    #         l = list(L(x))
    #         S = sum(l) / len(l)
    #         return vector([a - S for a in list(L(x))])

    #     return c

    # @functools.cache
    # def OrderCastMap(self, O=None):
    #     """
    #     return a function that, on x, returns a vector of rationals that
    #     gives the coordinates of x with respect to the specified order. x
    #     is integral if all coordinates are integers.
    #     """
    #     K = self.K
    #     if O is not None:
    #         assert O.number_field() == K
    #     else:
    #         O = K.maximal_order()
    #     B = [vector(v) for v in O.basis()]
    #     V = QQ**K.degree()
    #     W = V.submodule_with_basis(B)
    #     return lambda x: W.coordinate_vector(vector(x))

    def _schirokauer_maps_internal(self, ell, uniformizer_mode=False):
        K = self.K
        OK = K.maximal_order()
        f = K.defining_polynomial()
        lc = f.leading_coefficient()
        alpha_hat = K.gen() * lc
        f_hat = alpha_hat.minpoly()

        disc = f.discriminant()
        per_ideal = [(I, m, I.residue_class_degree(), K.uniformizer(I))
                     for I, m in OK.fractional_ideal(ell).factor()]

        # It takes something to make sure that we have an algebraic
        # integer with zero valuation at each of the prime ideals above
        # ell.
        # We'll make our life easier by asserting that lc is coprime to
        # ell. If it isn't, then life hates us. But honestly, we're
        # pretty much never going to encounter this situation anyway.
        #
        # while we're at it, let's also assume that it does not divide
        # the index of the equation order.
        assert gcd(lc, ell) == 1
        assert valuation(disc, ell) < 2

        # When we raise to the lcm({ell^f-1})-th power, we reach
        # something that is 1 modulo all _prime_ ideals above ell. If ell
        # ramifies, their product is not quite ell*OK, and this is really
        # an essential difference. There no further information mod ell
        # that is to be hoped.
        internal = []

        # before we use the ideal uniformizer as a way to rectify
        # possible valuations, we must pay attention to the fact that
        # _of course_, by construction, the uniformizer maps to
        # something of non-zero valuation down below, which is
        # precisely what we're trying to kill (the question of
        # whether the valuation is itself a multiple of ell is none
        # of our concern: it's dealt with in another way with the
        # ideal factorizations).
        prec = 2 + (1 if uniformizer_mode else 0)
        Zell = Zp(ell, prec)
        f_down = f_hat.change_ring(Zell)
        for g_down, m in f_down.factor():
            # match the corresponding ideal
            g_down_ell = g_down.change_ring(GF(ell))
            g_down_Z = g_down_ell.change_ring(ZZ)
            g_down_Z_radical = g_down_Z.factor()[0][0]
            match = [ (I, u)
                      for I, e, f, u in per_ideal
                      if g_down_Z_radical(alpha_hat) in I ]
            assert len(match) == 1
            I, u = match[0]

            if g_down.degree() == 1:
                # easy peasy.
                # alpha_hat maps to r, so alpha maps to r/lc
                alpha_hat_down = g_down.roots()[0][0]
                degree_down = 1

            elif r := g_down.roots(GF(ell)):
                # we're degree 1, but not in the "easy" case, which
                # means that we ramify. So alpha is going to r, but
                # in order to properly deal with the ramification
                # situation, we have to work a bit.  First convert
                # the factor of f_down to Eisenstein form, then convert
                # back...
                r = ZZ(r[0][0])
                x = g_down.parent().gen()
                uname = f"shifted_{K.variable_name()}_hat_ell"
                L = Zell.extension(g_down(x + r), names=(uname))
                # g_down(x+r)(u) = g_down(u+r) is zero, so alpha_hat maps
                # to u+r in this field
                alpha_hat_down = L.gen() + r
                degree_down = 1

            else:
                # we have a factor of degree larger than one with no root
                # mod ell. Since we mandate that the valuation of the
                # discriminant be at most one, we can't have ramification
                # here. Admittedly, it's slightly cheating.
                uname = f"{K.variable_name()}_hat_ell"
                L = Zell.extension(g_down, names=(uname))
                # alpha_hat maps to u in this field
                alpha_hat_down = L.gen()
                degree_down = g_down.degree()

            assert f_hat(alpha_hat_down) == 0
            alpha_down = alpha_hat_down / lc
            assert f(alpha_down) == 0

            hom = K.Hom(GF(ell)**degree_down, Sets())

            internal.append((hom, I, u, alpha_down, g_down))
        return internal

    def _schirokauer_map_helper(self, z, alpha_down, fix_valuation=False):
        ZP = ZZ['x']
        assert z != 0
        degree_down = alpha_down.parent().inertia_degree()
        ell = alpha_down.parent().residue_characteristic()
        z_down = z.polynomial()(alpha_down)
        assert z.minpoly()(z_down) == 0
        if fix_valuation:
            assert z_down.valuation() == 1
            z_down /= z_down.parent().uniformizer()
        else:
            assert z_down.valuation() == 0
        ze = (z_down**(ell**degree_down - 1)).expansion()
        assert ZP(ze[0]) == 1
        Pl = [ZP(ze[1])[i] for i in range(degree_down)]
        return (GF(ell)**degree_down)(Pl)

    def _schirokauer_map_second_helper(self, I, u, alpha_down, cui, z):
        v = z.valuation(I)
        numen = self._schirokauer_map_helper(z / u**v, alpha_down)
        denom = - v * cui
        return numen + denom

    def schirokauer_maps(self, ell):
        internal = self._schirokauer_maps_internal(ell)
        internal_uniformizer = self._schirokauer_maps_internal(ell, True)

        cu = [self._schirokauer_map_helper(u, alpha_down, True)
              for _, _, u, alpha_down, _ in internal_uniformizer]

        maps_to_fields = []

        for i in range(len(internal)):
            hom, I, u, alpha_down, g_down = internal[i]
            cui = cu[i]

            # caveat emptor: you must not return closures from a loop.
            # Use functools.partial instead.
            maps_to_fields.append(hom(functools.partial(
                self._schirokauer_map_second_helper,
                I, u, alpha_down, cui)))

        return maps_to_fields

    def schirokauer_map_simple(self, ell):
        K = self.K
        n = K.degree()
        hom = K.Hom(GF(ell)**n, Sets())
        OK = K.maximal_order()
        e_and_f_s = [(m, I.residue_class_degree())
                     for I, m in OK.fractional_ideal(ell).factor()]
        assert prod([e == 1 for e, f in e_and_f_s])
        Zell = Zp(ell, 2)
        ZP = ZZ['x']
        Kell = Zell['x'].quotient(K.defining_polynomial())
        expo = lcm([ell**f - 1 for e, f in e_and_f_s])

        def c(z):
            assert z.norm().valuation(ell) == 0
            ze = [ZP(list(l.expansion(start_val=0)))
                  for l in list(z.polynomial()(Kell.gen())**expo)]
            assert ze[0][0] == 1
            assert prod([ze[i][0] == (i == 0) for i in range(n)])
            Pl = [ZP(ze[i])[1] for i in range(n)]
            return (GF(ell)**n)(Pl)

        return hom(c)

    def real_representation_of_embeddings(self, z, prec=53):
        r"""
        given a number field element z return an n-uple of real numbers
        that has the same L2 norm as the vector $\left(z_{\theta_1},
        \ldots, z(\theta_{r_1}), z(\theta_{r_1+1}),
        z(\bar{\theta_{r_1+1}}), \ldots\right)$

        TESTS:
            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^3 - 2);
            sage: from cado_sage import CadoNumberFieldWrapper
            sage: N = CadoNumberFieldWrapper(F)
            sage: z = alpha^2 - 12 * alpha + 47
            sage: Rz = N.real_representation_of_embeddings(z)
            sage: pl = N.places()
            sage: r1, r2 = F.signature()
            sage: eR = [p(z) for p in pl[:r1]]
            sage: eC = flatten([(p(z), p(z).conjugate()) for p in pl[r1:]])
            sage: tensor_C = vector(eR + eC)
            sage: ratio = tensor_C.norm(2) / Rz.norm(2)
            sage: ratio > 0.99 and ratio < 1.01
            True
        """

        K = self.K
        pl = self.places(precision=prec, interval_based=True)
        r1, r2 = K.signature()
        embs = [ P(z) for P in pl[:r1] ]
        sqrt2 = sqrt(pl[0].codomain()(2))
        embs += flatten([list(sqrt2 * P(z)) for P in pl[r1:]])
        if any(x.diameter() > 1 for x in embs):
            raise FloatingPointError("insufficient precision")
        return vector(embs)

    def all_complex_roots(self, prec=53):
        K = self.K
        pl = self.places(all_complex=True, precision=prec)
        r1, r2 = K.signature()
        z = K.gen()
        embs = [ P(z) for P in pl[:r1] ]
        embs += flatten([(x,x.conjugate())
                         for x in [P(z) for P in pl[r1:]]])
        return vector(embs)

    def modules_of_embeddings_from_log_embeddings(self, logs):
        r"""
        given the log embeddings as returned by LogMap, returns the
        absolute values of all $n$ embeddings, which can a priori be
        compared to the vector returned by
        real_representation_of_embeddings (that is, we won't necessarily
        get a vector with only ones and minus ones because we're going to
        have cosines and sines as well, but in any case we should get a
        vector whose L2-norm is exactly $\sqrt{n}$ ---it's the rotation
        of the all ones vector by a unitary matrix.)

        TESTS:
            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^3 - 2);
            sage: n = F.degree()
            sage: from cado_sage import CadoNumberFieldWrapper
            sage: N = CadoNumberFieldWrapper(F)
            sage: z = alpha^2 - 12 * alpha + 47
            sage: Rz = N.real_representation_of_embeddings(z)
            sage: Lz = N.LogMap()(z)
            sage: M = N.modules_of_embeddings_from_log_embeddings(Lz)
            sage: ratio = vector([Rz[i] / M[i] for i in range(n)]).norm() / float(sqrt(n))
            sage: ratio > 0.99 and ratio < 1.01
            True

        """
        K = self.K
        r1, r2 = K.signature()

        rr = [exp(r) for r in logs[:r1]]

        # there's really a nasty thing with 2.
        cc = flatten([(exp(c/2),)*2 for c in logs[r1:]])

        return vector(rr + cc)

    def places(self, precision=53, all_complex=False, interval_based=False):
        """
        K.places() is practically unusable:
         - all_complex happens to be its _first_ argument, so if you call
           K.places(80), you really don't get what you would expect
         - if we _really_ want interval arithmetic, there's no reason why
           we would not to _also_ specify a precision value
         - RDF and CDF seem to go back to the libm and _not_ get correct
           rounding at all, which leads to all sorts of surprises (like
           roots being wrong by several dozens of ulps).
        """

        if interval_based:
            R = RealIntervalField(precision)
            C = ComplexIntervalField(precision)
        else:
            R = RealField(precision)
            C = ComplexField(precision)

        f = self.K.defining_polynomial()

        # first, find the intervals with roots, and see how much
        # precision we need to approximate the roots
        all_intervals = [x[0] for x in f.roots(C)]

        # first, set up the real places
        if interval_based:
            if all_complex:
                real_intervals = [z for z in all_intervals if 0 in z.imag()]
            else:
                real_intervals = [x for x,y in all_intervals if 0 in y]
        else:
            if all_complex:
                real_intervals = [z for z in all_intervals if z.imag().is_zero()]
            else:
                real_intervals = [x for x,y in all_intervals if y.is_zero()]

        real_places = [self.K.hom([i], check=False) for i in real_intervals]

        complex_places = [self.K.hom([i], check=False)
                          for i in all_intervals if i.imag() > 0]

        return real_places + complex_places


class CadoNumberTheory(object):
    def __init__(self, poly):
        self.poly = poly
        self._wrappers = [ CadoNumberFieldWrapper(K) for K in self.poly.K ]

    def wrappers(self):
        return self._wrappers

    def __getitem__(self, i):
        return self._wrappers[i]

    def maximal_orders(self):
        return [ K.maximal_order() for K in self.poly.K ]

    def J(self):
        return [ OK.fractional_ideal(1, OK.number_field().gen())**-1
                 for OK in self.maximal_orders() ]
