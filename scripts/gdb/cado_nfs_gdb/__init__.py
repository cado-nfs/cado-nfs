import os
import gdb

print("## python code is loaded from the following location:")
print(f"## {os.path.dirname(__file__)}")

from . import base          # noqa: E402
from . import display       # noqa: E402
from . import integers      # noqa: E402
from . import polynomials   # noqa: E402
from . import matrices      # noqa: E402

if base.debug:
    import importlib
    importlib.reload(base)
    importlib.reload(display)
    importlib.reload(integers)
    importlib.reload(polynomials)
    importlib.reload(matrices)


def register_cado_nfs_printers(obj):
    # from .base import register_cado_nfs_printers
    gdb.printing.register_pretty_printer(obj, base.cado_nfs_printer)
