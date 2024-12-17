import time


class FailedDescent(Exception):
    def __init__(self, failed, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.failed = failed
    def __str__(self):
        return "Failed descents for: " + ", ".join(self.failed)


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args,
                                                                 **kwargs)
        return cls._instances[cls]


class HwlocFeatureFlag(metaclass=Singleton):
    def __init__(self):
        self.t = None

    def __call__(self):
        return self.t

    def set(self, t):
        self.t = t


def feature_set_hwloc(x):
    HwlocFeatureFlag().set(x)


def feature_get_hwloc():
    return HwlocFeatureFlag()()


def check_result(two, log2, z, logz, p, ell):
    assert (p - 1) % ell == 0
    h = ((p - 1) // ell)
    assert pow(z, log2 * h, p) == pow(2, logz * h, p)
    print("Final consistency check ok!")


def a_over_b_mod_p(a, b, p):
    if b % p == 0:
        return p
    ib = pow(b, p - 2, p)
    return (a * ib) % p


def is_a_over_b_equal_r_mod_pk(a, b, rk, p, pk):
    if b % p != 0:   # non-projective
        if rk >= pk:
            return False
        return (a - b * rk) % pk == 0
    else:  # projective
        if rk < pk:
            return False
        return (b - a * rk) % pk == 0


class object_holder():
    def __init__(self, v=None):
        self.v = v

    def set(self, v):
        self.v = v

    def unset(self):
        self.v = None


def prime_ideal_mixedprint(pr):
    p = pr[0]
    side = pr[1]
    if side == 0:
        machine = "%x 0 rat" % p
        human = "0,%d" % p
    else:
        r = pr[2]
        machine = "%x %d %x" % pr
        human = "%d,%d,%d" % (side, p, r)
    return machine, human


# http://stackoverflow.com/questions/107705/disable-output-buffering
# shebang takes only one arg...
# python3 doesn't grok sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
# setting PYTHONUNBUFFERED here is too late.
class drive_me_crazy(object):
    def __init__(self, stream, timestamp=False):
        self.stream = stream
        self.eol = 1
        self.timestamp = timestamp

    def write(self, data):
        if self.timestamp:
            p = 0
            while len(data) > p:
                d = data.find('\n', p)
                if d < 0:
                    break
                if self.eol:
                    self.stream.write(time.asctime() + " ")
                self.stream.write(data[p:d + 1])
                self.eol = True
                p = d + 1
            if len(data) > p:
                if self.eol:
                    self.stream.write(time.asctime() + " ")
                self.stream.write(data[p:])
                self.eol = False
        else:
            self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)
