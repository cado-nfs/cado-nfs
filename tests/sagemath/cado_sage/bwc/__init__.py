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

from .parameters import BwcParameters  # noqa: F401
from .matrix import BwcMatrix  # noqa: F401
from .balancing import BwcBalancing  # noqa: F401
from .vector import BwcVector  # noqa: F401
from .vector_set import BwcVectorSet  # noqa: F401
from .x_vector import BwcXVector  # noqa: F401
from .check_data import BwcCheckData  # noqa: F401
from .a_files import BwcAFiles  # noqa: F401
from .f_files import BwcFFiles  # noqa: F401
from .right_hand_side import BwcRightHandSide  # noqa: F401
from .s_vector import BwcSVector  # noqa: F401
from .s_vector_set import BwcSVectorSet  # noqa: F401
