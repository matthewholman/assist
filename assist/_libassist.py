"""
Internal C library loader and shared helpers.

This module exists to avoid circular imports:
- `assist/__init__.py` re-exports public objects like `Ephem`.
- `assist/ephem.py` needs access to `clibassist` and `assist_error_messages`.

If `assist/ephem.py` imports `assist` (the package top-level) at runtime,
it creates an import cycle. Importing from `assist._libassist` avoids that.
"""

from __future__ import annotations

import os
import sysconfig
from ctypes import c_char_p, c_int, cdll


def _ext_suffix() -> str:
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    return suffix if suffix is not None else ".so"


suffix = _ext_suffix()

# Import shared C library
pymodulespath = os.path.dirname(__file__)

# Allow tests / power users to force a specific libassist binary, so Python and C
# test suites can be run against the *exact same* compiled library.
#
# If a relative path is provided, interpret it relative to the repo root
# (one directory above the `assist/` package).
_override = os.environ.get("ASSIST_LIBASSIST_PATH")
if _override:
    if not os.path.isabs(_override):
        _repo_root = os.path.abspath(os.path.join(pymodulespath, os.pardir))
        _override = os.path.join(_repo_root, _override)
    clibassist = cdll.LoadLibrary(_override)
else:
    clibassist = cdll.LoadLibrary(pymodulespath + "/../libassist" + suffix)

__libpath__ = str(getattr(clibassist, "_name", ""))


def _libassist_str(name: str) -> str:
    """Extract a string from libassist by symbol name."""
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
    ccpp = c_char_p * e_N
    message = ccpp.in_dll(clibassist, "assist_error_messages")[e]
    if message is None:
        raise RuntimeError("assist_error_messages is missing from libassist. Please report this issue on GitHub.")
    return message.decode("ascii")

