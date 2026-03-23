# -*- coding: utf-8 -*-

# Public API re-exports.
#
# Keep all `ctypes` / libassist loading details in `assist._libassist` to avoid
# circular imports (e.g. `assist.ephem` needs `clibassist` but should not import
# the package top-level during import time).
from ._libassist import (
    __build__,
    __githash__,
    __libpath__,
    __version__,
    assist_error_messages,
    clibassist,
)

from .ephem import Ephem
from .tools import simulation_convert_to_rebound

# Do not change the following line. Will be updated automatically with update_version.py
moduleversion = '1.1.9'
libassistversion = __version__
if moduleversion != libassistversion:
    print("WARNING: python module and libassist have different version numbers: ", moduleversion, libassistversion)

# Avoid importing heavy optional dependencies (like numpy via `assist.extras`)
# at package import time. This keeps `import assist` lightweight and also
# makes relocation/linking tests easier to isolate.
def __getattr__(name: str):  # pragma: no cover
    if name == "Extras":
        from .extras import Extras
        return Extras
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "__libpath__",
    "__version__",
    "__build__",
    "__githash__",
    "assist_error_messages",
    "clibassist",
    "Extras",
    "Ephem",
    "simulation_convert_to_rebound",
]
