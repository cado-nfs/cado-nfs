
import re
import subprocess

OK = "ok ✅"
NOK = "❌"
EXCL = "❗"
NOTHING_TO_DO = "💤"
HURRAH = "🥳"


# The functions below are used to parse binary data, which is found quite
# often in the Bwc files, for instance (but these functions can also be
# of more general interest)

def get_int(f, bytes=4, signed=False, repeat=None, may_fail=False):
    assert not (may_fail and repeat is not None)
    if repeat is not None:
        return [get_int(f, bytes=bytes, signed=signed) for i in range(repeat)]
    data = f.read(bytes)
    if not data and may_fail:
        return None
    assert data
    return int.from_bytes(data, 'little', signed=signed)


def u32(f, *args, **kwargs):
    return get_int(f, bytes=4, *args, **kwargs)


def u64(f, *args, **kwargs):
    return get_int(f, bytes=8, *args, **kwargs)


def s32(f, *args, **kwargs):
    return get_int(f, bytes=4, signed=True, *args, **kwargs)


def s64(f, *args, **kwargs):
    return get_int(f, bytes=8, signed=True, *args, **kwargs)


def cat_or_zcat(filename):
    """
    return an iterable over the file contents, which my entail
    decompression if the file happens to be compressed
    """
    if re.search(r"\.gz$", filename):
        return subprocess.Popen(["zcat", filename],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).stdout
    else:
        return open(filename)


# define a singleton class for a verbose flag that we use throughout the
# cado-nfs sage code.

class __verbose_flag(object):
    def __init__(self):
        self.v = False

    def __bool__(self):
        return self.v


cado_verbose = __verbose_flag()


def get_verbose():
    return bool(cado_verbose)


def set_verbose(flag=True):
    cado_verbose.v = flag
