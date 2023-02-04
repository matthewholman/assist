from . import clibassist
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import warnings
from .ephem import Ephem

ASSIST_FORCES = {
    "SUN"                : 0x01,
    "PLANETS"            : 0x02,
    "ASTEROIDS"          : 0x04,
    "NON_GRAVITATIONAL"  : 0x08,
    "EARTH_HARMONICS"    : 0x10,
    "SUN_HARMONICS"      : 0x20,
    "GR_EIH"             : 0x40,
    "GR_SIMPLE"          : 0x80,
    "GR_POTENTIAL"       : 0x100,
} 

class Extras(Structure):
    """
    Main object used for all ASSIST operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct assist_extras.
    """
    
    def __init__(self, sim, ephem):
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected_ 
        clibassist.assist_init(byref(self), byref(sim), byref(ephem))
        self.extras_should_free_ephem = 0

    def __del__(self):
        clibassist.assist_free_pointers(byref(self))

    def detach(self, sim):
        sim._extras_ref = None # remove reference to assist so it can be garbage collected 
        clibassist.assist_detach(byref(sim), byref(self))

    def integrate_or_interpolate(self, t):
        clibassist.assist_integrate_or_interpolate(byref(self), c_double(t))
    
    @property
    def forces(self):
        l = []
        for k in ASSIST_FORCES:
            if self._forces & ASSIST_FORCES[k]:
                l.append(k)
        return l
    @forces.setter
    def forces(self, value):
        if not isinstance(value, list):
            raise AttributeError("Forces need to be a list.")
        for v in value:
            if not isinstance(v, str):
                raise AttributeError("Each force needs to be a string.")
            if v.upper() not in ASSIST_FORCES:
                raise AttributeError("Force '"+v+"' not recognized. Needs to be one of the following: "+", ".join(ASSIST_FORCES))
        v = 0
        for k in ASSIST_FORCES:
            if k in value:
                v = v | ASSIST_FORCES[k]
        self._forces = c_int(v)

    @property
    def particle_params(self):
        raise AttributeError("Cannot get particle_params. Only setting is supported.")

    @particle_params.setter
    def particle_params(self, value):
        value_p = value.ctypes.data_as(POINTER(c_double))
        self._particle_params = value_p



    _fields_ =  [("_sim", POINTER(rebound.Simulation)),
                 ("ephem", POINTER(Ephem)),
                 ("_ephem_cache", c_void_p),
                 ("extras_should_free_ephem", c_int),
                 ("geocentric", c_int),
                 ("last_state", POINTER(rebound.Particle)),
                 ("current_state", POINTER(rebound.Particle)),
                 ("_particle_params", POINTER(c_double)),
                 ("steps_done", c_int),
                 ("_forces", c_int),
        ]

# avoid circular imports
