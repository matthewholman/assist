import warnings
from ctypes import (CFUNCTYPE, POINTER, Structure, byref, c_char, c_char_p,
                    c_double, c_int, c_long, c_uint, c_uint32, c_ulong,
                    c_void_p, cast)
from typing import Optional, Union

import rebound

import assist

from . import assist_error_messages, clibassist

ASSIST_BODY_IDS = {
        0: "Sun",
        1: "Mercury",
        2: "Venus",
        3: "Earth",
        4: "Moon",
        5: "Mars",
        6: "Jupiter",
        7: "Saturn",
        8: "Uranus",
        9: "Neptune",
        10: "Pluto",
	    11: "Camilla",
        12: "Ceres",
        13: "Cybele",
        14: "Davida",
	    15: "Eunomia",
        16: "Euphrosyne",
        17: "Europa",
        18: "Hygiea",
        19: "Interamnia",
        20: "Iris",
        21: "Juno",
        22: "Pallas",
        23: "Psyche",
        24: "Sylvia",
        25: "Thisbe",
        26: "Vesta",
        }

class Ephem(Structure):
    """
    Main object used for all ASSIST Ephemeris operations.
    Not tied to a particular REBOUND simulation.
    """

    def __init__(self, planets_path: Optional[str] = None, asteroids_path: Optional[str] = None):
        if planets_path is not None:
            c_planets_path = c_char_p(str(planets_path).encode("ascii"))
        else:
            c_planets_path = c_char_p(None)

        if asteroids_path is not None:
            c_asteroids_path = c_char_p(str(asteroids_path).encode("ascii"))
        else:
            c_asteroids_path = c_char_p(None)

        clibassist.assist_ephem_init.restype = c_int
        ret = clibassist.assist_ephem_init(byref(self), c_planets_path, c_asteroids_path)
        if ret != 0:
            raise RuntimeError(assist_error_messages(ret))

    def get_particle(self, body: Union[int, str], t: float) -> rebound.Particle:
        if isinstance(body, str):
            body_str = body.lower()
            body = -1
            for k in ASSIST_BODY_IDS:
                if body_str == ASSIST_BODY_IDS[k].lower():
                    body = k
            if body < 0:
                raise ValueError("Cannot find body '" + body_str + "'. Needs to be one of: " + ", ".join([ASSIST_BODY_IDS[k] for k in ASSIST_BODY_IDS]) + ".")
        if not isinstance(body, int):
            raise ValueError("Expecting integer for body id.")

        clibassist.assist_get_particle_with_error.restype = rebound.Particle
        e = c_int(0)
        p = clibassist.assist_get_particle_with_error(
            byref(self), c_int(body), c_double(t), byref(e)
        )
        if e.value:
            raise RuntimeError(assist_error_messages(e.value))
        return p

    def __del__(self) -> None:
        clibassist.assist_ephem_free_pointers(byref(self))

    _fields_ = [
        ("jd_ref", c_double),
        ("spk_global", c_void_p),
        ("spk_planets", c_void_p),
        ("spk_asteroids", c_void_p),
    ]

