from . import clibassist
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import warnings
from .ephem import Ephem

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
    
    _fields_ =  [("_sim", POINTER(rebound.Simulation)),
                 ("ephem", POINTER(Ephem)),
                 ("_ephem_cache", c_void_p),
                 ("extras_should_free_ephem", c_int),
                 ("geocentric", c_int),
                 # not complete....
        ]

# avoid circular imports
