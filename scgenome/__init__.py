""" Single-Cell Whole Genome Analysis in Python. """

from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets

from . import _version
__version__ = _version.get_versions()['version']

# has to be done at the end, after everything has been imported
import sys

sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['tl', 'pp', 'pl']})

del sys
