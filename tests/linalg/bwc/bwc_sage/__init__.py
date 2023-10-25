"""
This package is work in progress. As of late 2023, it is still far from
the functionality level of the original magma script

The work plan is:
    1. prime field case (nullspace=right), with arbitrary balancing
        ✅ This works for a few automated tests, at least as far lingen
    2. we want to use this to debug the double matrix setting
    3. extend to the binary case (nullspace=left). We know that there are
       a few monsters down that path, like m4ri data being hard to
       initialize totally freely from within sage
    4. address the remaining cases like interleaving and so on.
"""

from .BwcParameters import BwcParameters  # noqa: F401
from .BwcMatrix import BwcMatrix  # noqa: F401
from .BwcBalancing import BwcBalancing  # noqa: F401
from .BwcVector import BwcVector  # noqa: F401
from .BwcVectorSet import BwcVectorSet  # noqa: F401
from .BwcXVector import BwcXVector  # noqa: F401
from .BwcCheckData import BwcCheckData  # noqa: F401
from .BwcAFiles import BwcAFiles  # noqa: F401
from .BwcFFiles import BwcFFiles  # noqa: F401
from .BwcSVector import BwcSVector  # noqa: F401
from .BwcSVectorSet import BwcSVectorSet  # noqa: F401
