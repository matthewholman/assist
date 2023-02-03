from . import clibassist
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import assist
import warnings

class Ephem(Structure):
    """
    Main object used for all ASSIST Ephemeris operations.
    Not tied to a particular REBOUND simulation.
    """
    
    def __new__(cls, planets_path=None, asteroids_path=None):
        assist = super(Ephems,cls).__new__(cls)
        clibassist.assist_ephem_init(byref(assist), c_char_p(planets_path.encode("ascii")), c_char_p(asteroids_path.encode("ascii")))
        return assist

    def __init__(self, sim, filename=None):
        pass

    def __del__(self):
        clibassist.assist_ephem_free_pointers(byref(self))

    _fields_ =  [("_sim", POINTER(rebound.Simulation)),
            ]

