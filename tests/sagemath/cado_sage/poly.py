import os
import re
import sys
import math

# from sage.rings.integer_ring import ZZ
# from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.real_mpfr import RR
from sage.rings.number_field.number_field import NumberField
from .tools import get_verbose
from .tools import NOK
from .number_theory import CadoNumberTheory


class CadoPolyFile(object):
    """
    This class contains stuff to read a polynomial pair from a cado .poly
    file. Some effort goes into adapting to the various ways a poly file
    can be laid out, but these are certainly a few gaps.

    EXAMPLES:

    sage: from cado_sage import CadoPolyFile
    sage: data = '''
    ....: n: 123412341234618726123412341234611
    ....: skew: 418.982
    ....: c0: 1312885934640
    ....: c1: -181216146
    ....: c2: -19457299
    ....: c3: 280
    ....: Y0: -4970288929
    ....: Y1: 19570577
    ....: # MurphyE (Bf=1.311e+05,Bg=1.311e+05,area=6.554e+09) = 2.440e-04
    ....: # f(x) = 280*x^3-19457299*x^2-181216146*x+1312885934640
    ....: # g(x) = 19570577*x-4970288929
    ....: '''
    sage: p = CadoPolyFile()
    sage: p._CadoPolyFile__read(data.split('\n'))
    sage: print(p)
    NFS configuration for modulus 123412341234618726123412341234611
      K0: Number Field in alpha0 with defining polynomial
      19570577*x - 4970288929
      K1: Number Field in alpha1 with defining polynomial
      280*x^3 - 19457299*x^2 - 181216146*x + 1312885934640


    """
    def __init__(self,
                 filename=None,
                 wdir=None):

        self.filename = filename
        if wdir is None and filename is not None:
            self.wdir = os.path.dirname(self.filename)
        else:
            self.wdir = wdir

        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self.K = []
        self.f = []
        self.skewness = 0
        self.N = 0
        self.nt = None
        self.m = None

    def _check_mapping(self):
        # we only know how to do it when we have a polynomial of degree 1
        # somewhere.
        rat = [i for i, f in enumerate(self.f) if f.degree() == 1]
        assert len(rat)
        u, v = list(self.f[rat[0]])
        m, = -u/v
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.integer_ring import ZZ
        m = Integers(self.N)(m)
        assert all([f(m) == 0 for f in self.f])
        return ZZ(m)

    def read(self):
        try:
            if get_verbose():
                print(f"Reading {self.filename}")
            with open(self.filename, "r") as fm:
                self.__read(fm)

        except Exception as e:
            # We're really in bad shape if an exception occurs here.
            # We're not even trying to salvage the BwcMatrix object, as
            # the error is most probably obvious.
            print(f"Exception while reading {self.filename} {NOK}",
                  file=sys.stderr)
            raise e

        return self

    def __parseline(self, line):
        """

        EXAMPLES::

            sage: from cado_sage import CadoPolyFile
            sage: p = CadoPolyFile()
            sage: p._CadoPolyFile__parseline("n=-47")
            sage: p.N
            -47

            sage: p._CadoPolyFile__parseline("poly1=280*x^3-19*x^2-18*x+13")
            sage: p.f[1]
            280*x^3 - 19*x^2 - 18*x + 13

        """

        from sage.rings.integer_ring import ZZ
        ZP = ZZ['x']
        x = ZP.gen()

        # regexps for integers, real numbers, equal sign, and lists of
        # coefficients.
        nb_nocapture = r"-?\d+"
        nbf_nocapture = r"-?\d+(?:\.\d*)?"
        nb = f"({nb_nocapture})"
        nbf = f"({nbf_nocapture})"
        eq = r"(?::\s*|\s*=\s*)"
        coeffs = rf"({nb_nocapture}(?:\s*,\s*{nb_nocapture})*)"

        if not line or line.startswith("#"):
            return
        elif m := re.match(rf"^[nN]{eq}{nb}$", line):
            self.N = ZZ(m.group(1))
        elif m := re.match(rf"^[m]{eq}{nb}$", line):
            self.m = ZZ(m.group(1))
        elif m := re.match(rf"^skew{eq}{nbf}$", line):
            self.skewness = RR(m.group(1))
        elif m := re.match(rf"^c(\d+){eq}{nb}$", line):
            if not self.f:
                self.f = [ZP(0)] * 2
            self.f[1] += ZZ(m.group(2)) * x**int(m.group(1))
        elif m := re.match(rf"^Y(\d+){eq}{nb}$", line):
            if not self.f:
                self.f = [ZP(0)] * 2
            self.f[0] += ZZ(m.group(2)) * x**int(m.group(1))
        elif m := re.match(rf"^poly(\d+){eq}{coeffs}$", line):
            i = int(m.group(1))
            coeffs = [ZZ(c) for c in m.group(2).split(',')]
            while i >= len(self.f):
                self.f.append(ZP(0))
            self.f[i] = ZP(coeffs)
        elif m := re.match(rf"^poly(\d+){eq}(.*)$", line):
            i = int(m.group(1))
            while i >= len(self.f):
                self.f.append(ZP(0))
            self.f[i] = ZP(m.group(2))
        else:
            raise ValueError(f"Cannot parse line in poly file: {line}")

    def __create_fields(self):
        # avoid fully factoring the discriminant, just do trial division
        # up to 10**7
        fdisc = [f.discriminant().factor(limit=10**7) for f in self.f]
        round2_primes = [[p for p, e in dd if e > 1] for dd in fdisc]
        self.K = [NumberField(f,
                              maximize_at_primes=round2_primes[i],
                              names=f"alpha{i}")
                  for i, f in enumerate(self.f)]

    def __check_resultant(self):
        """
        Check that the resultant is a multiple of N. If we want to extend
        this to the multiple number field setting, there's a bit of a
        catch here because we really want _all_ polynomials to intersect
        at the same root above N (not just pairwise). E.g. x^2-3x+2,
        x^2-x-2, and x^2+N-1 should not pass this test, yet they
        currently do.
        """
        f = self.f[0]
        for i, g in enumerate(self.f[1:]):
            if f.resultant(g) % self.N != 0:
                raise ValueError("resultant of polynomials 0 and {i}" +
                                 " is not zero modulo N")
            if f(self.m) % self.N != 0:
                raise ValueError(f"polynomial {i} does not"
                                 f" vanish mod {self.N}: {f(self.m)}")

    def __read(self, contents):
        self.__clear_fields_for_read()

        for line in contents:
            self.__parseline(line)

        if len(self.f) < 2:
            raise ValueError("poly file must" +
                             " contain at least two polynomials")

        if self.m is None:
            self.m = self._check_mapping()

        self.__create_fields()
        self.__check_resultant()
        self.nt = CadoNumberTheory(self)

    def __str__(self):
        if not self.N:
            if self.filename:
                return f"NFS configuration, to be read from {self.filename}"
            else:
                return "NFS configuration (empty placeholder)"
        else:
            rep = f"NFS configuration for modulus {self.N}\n"
            for i, K in enumerate(self.K):
                rep += f"  K{i}: {K}\n"
        return rep

    def _normalized_polynomial_for_graphics(self):
        pi = RR.pi()

        f1 = self.f[1]
        x = f1.parent().gen()
        d = f1.degree()
        f1s = f1(self.skewness*x)/self.skewness**(d/2)

        # find max on a circle.
        S = 100
        B = 0
        for i in range(1, S):
            t = i * pi / S
            c = math.cos(t)
            s = math.sin(t)
            B = max(B, abs(f1s.homogenize()(c, s)))

        # return f together with its max over a circle

        return f1s.homogenize(), B

    def create_2d_valley_plot(self, output_filename,
                              style='heatmap',
                              heatmap_plot_points=500,
                              scale='auto'
                              ):
        """
        This creates the frequently seen "spider web" graphics for this
        polynomial. It's normalized with respect to the skewness that is
        stored in the polynomial, and to the max value that the
        polynomial takes on the unit circle.

        style can be 'heatmap' or 'contour'.

        For style=='heatmap', the 'heatmap_plot_points' gives the size of
        the sampling grid.
        """

        if not re.match(r'.*\.png', output_filename):
            raise RuntimeError("output_filename must be a .png file")

        F, B = self._normalized_polynomial_for_graphics()

        if scale == 'auto':
            F = F / B
        else:
            F = F / scale

        z1 = 0
        z0 = -10

        def G(u, v):
            return max(math.log(abs(F(u, v))), z0)

        from sage.symbolic.ring import SR
        x, y = SR.var('x y')

        if style == 'contour':
            ni = 10
            dz = (z1 - z0) / ni
            levels = [z0 + i * dz for i in range(ni+1)]
            from sage.plot.contour_plot import contour_plot
            figure = contour_plot(G,
                                  (x, -1, 1), (y, -1, 1),
                                  contours=levels,
                                  fill=False,
                                  cmap="jet", colorbar=False)
        elif style == 'heatmap':
            from sage.plot.density_plot import density_plot
            figure = density_plot(G,
                                  (x, -1, 1), (y, -1, 1),
                                  cmap="jet",
                                  plot_points=heatmap_plot_points)
        else:
            raise RuntimeError("The 'style' parameter must be"
                               " either 'contour' or 'heatmap'")

        pm = figure.matplotlib(axes=False, axes_pad=0)
        pm.subplots_adjust(left=0, right=1, top=1, bottom=0)
        import matplotlib.pyplot
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        pm.set_canvas(FigureCanvasAgg(pm))
        pm.set_size_inches(5, 5)
        ax = matplotlib.pyplot.Axes(pm, [0., 0., 5., 5.])
        ax.set_axis_off()
        pm.add_axes(ax)
        pm.savefig(output_filename, bbox_inches=0)

    def create_3d_valley_plot(self, output_filename,
                              scale='auto',
                              resolution=300):
        """
        This creates a threejs animated graphics object (about 20MB in
        size for 300 points) that can show the valleys of log(|F(a,b)|).
        """

        F, B = self._normalized_polynomial_for_graphics()

        if scale == 'auto':
            F = F / B
        else:
            F = F / scale

        z1 = 0
        z0 = -10

        def G(u, v):
            return max(math.log(abs(F(u, v))), z0)

        def G0(u, v):
            return math.log(abs(F(u, v)))

        from sage.plot.colors import colormaps

        def col2(x, y):
            z = G0(x, y)
            return float((z - z0) / (z1 - z0))

        from sage.symbolic.ring import SR
        x, y = SR.var('x y')

        from sage.plot.plot3d.plot3d import plot3d

        surf = plot3d(G,
                      (x, -1, 1), (y, -1, 1),
                      color=(col2, colormaps.jet),
                      plot_points=resolution)

        bottom = plot3d(z0-1,
                        (x, -1, 1), (y, -1, 1),
                        color=(col2, colormaps.jet),
                        plot_points=resolution,
                        opacity=0.5)

        Gr = (surf + bottom)

        if re.match(r'.*\.html', output_filename):
            Gr.save(output_filename,
                    online=True,
                    aspect_ratio=[1, 1, 2/(z1-z0+1)])
        else:
            Gr.save(output_filename,
                    online=False,
                    aspect_ratio=[1, 1, 2/(z1-z0+1)])

        # os.system(rf'sed -e "s/<script src=\".*three.min.js\">'
        #           rf'<\/script>/<script src=\"three.js\"><\/script>/"'
        #           rf' -i {base}-valleys.html')
