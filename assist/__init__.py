# -*- coding: utf-8 -*-

#Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

#Import shared C library
import os
pymodulespath = os.path.dirname(__file__)
from ctypes import *
clibassist = cdll.LoadLibrary(pymodulespath + '/../libassist' + suffix)

# Version
__version__ = c_char_p.in_dll(clibassist, "assist_version_str").value.decode('ascii')

# Build
__build__ = c_char_p.in_dll(clibassist, "assist_build_str").value.decode('ascii')
# Check for version

# Githash
__githash__ = c_char_p.in_dll(clibassist, "assist_githash_str").value.decode('ascii')

try:
    import pkg_resources
    moduleversion = pkg_resources.require("assist")[0].version
    libassistversion = __version__
    if moduleversion != libassistversion:
        print("WARNING: python module and libassist have different version numbers: ", moduleversion, libassistversion)
except:
    pass    # this check fails in python 3. Problem with setuptools

from .ephem import Ephem
from .extras import Extras


__all__ = ["__libpath__", "__version__", "__build__", "__githash__", "Extras", "Ephem"]
