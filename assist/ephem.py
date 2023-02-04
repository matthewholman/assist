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
    
    def __init__(self, planets_path=None, asteroids_path=None):
        if planets_path is not None:
            planets_path = planets_path.encode("ascii")
        if asteroids_path is not None:
            asteroids_path = asteroids_path.encode("ascii")
        clibassist.assist_ephem_init.restype = c_int
        ret = clibassist.assist_ephem_init(byref(self), c_char_p(planets_path), c_char_p(asteroids_path))
        if ret != 0:
            if ret == 1:
                raise RuntimeError("JPL planet ephemeris file not found.")
            if ret == 2:
                raise RuntimeError("JPL asteroid ephemeris file not found.")

            raise RuntimeError("An error occured while creating the ephemeris structure.")

    def get_particle(self, i, t):
        clibassist.assist_get_particle.restype = rebound.Particle
        p = clibassist.assist_get_particle(byref(self), c_int(i), c_double(t))
        return p



    def __del__(self):
        clibassist.assist_ephem_free_pointers(byref(self))

    _fields_ =  [("jd_ref", c_double),
                 ("_pl", c_void_p),
                 ("_spl", c_void_p),
            ]

