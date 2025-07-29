import os
import inspect
from omfit_classes.utils_math import contourPaths
from omfit_classes.fluxSurface import boundaryShape
from omfit_classes.utils_base import is_int, compare_version, sanitize_version_number

src_import = """
# install minimal omift_classes used for ips-fastran
import re
import scipy
import matplotlib
from matplotlib import pyplot
import numpy as np
def _available_to_user_math(f): pass
"""
src = '\n'.join(
    [src_import,
    inspect.getsource(is_int),
    inspect.getsource(compare_version),
    inspect.getsource(sanitize_version_number),
    inspect.getsource(contourPaths),
    inspect.getsource(boundaryShape)]
)
src = src.replace("@", "#")

package = 'omfit_classes'
if not os.path.exists(package):
    print('make directory ', package)
    os.makedirs(package)

with open(os.path.join(package, 'fluxSurface.py'), 'w') as file:
    file.write(src)
