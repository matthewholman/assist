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


def _libassist_str(name: str) -> str:
    """Extract a string from libassist, referencing it by symbol
    name. Raise a RuntimeError if the name is not present.

    """
    val = c_char_p.in_dll(clibassist, name).value
    if val is None:
        raise RuntimeError(f"unable to find symbol {name} in libassist")
    return val.decode("ascii")

__version__ = _libassist_str("assist_version_str")

__build__ = _libassist_str("assist_build_str")

__githash__ = _libassist_str("assist_githash_str")


def assist_error_messages(e: int) -> str:
    e_N = c_int.in_dll(clibassist, "assist_error_messages_N").value
    if e >= e_N:
        raise RuntimeError("An error occured while trying to process an ASSIST error message.")
    ccpp  = c_char_p * e_N
    message = ccpp.in_dll(clibassist, "assist_error_messages")[e]
    if message is None:
        raise RuntimeError("assist_error_messages is missing from libassist. Please report this issue on GitHub.")
    return message.decode("ascii")

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
