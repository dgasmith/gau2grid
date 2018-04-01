"""
Gau2grid base init
"""

from . import RSH
from . import c_generator as c_gen
from . import codegen
from . import np_generator as np_gen
from . import order
from . import python_reference as ref

# Pull in code from the c wrapper
from .c_wrapper import collocation, collocation_basis, c_compiled, cgg_path

# Pull in tests
from .extras import test

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
