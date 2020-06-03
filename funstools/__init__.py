#
"""
Python toolkit for FUNS project

funstools is a toolkit that contains the useful functions and core tools needed
for performing FUNS project with Python. It is based on the 'numpy' and 'astropy'
packages, which have been validated by astronomers for a long time.
"""

from . import checkshape
from . import cube2map
from . import makemaps
from . import plots
from . import smooth
from . import userinput
from . import decompose

from .checkshape import *
from .cube2map import *
from .makemaps import *
from .plots import *
from .smooth import *
from .userinput import *
from .decompose import *
