from . import clibassist
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import assist
import warnings

class Extras(Structure):
    """
    Main object used for all ASSIST operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct assist_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    """
    
    def __new__(cls, sim, filename=None):
        assist = super(Extras,cls).__new__(cls)
        return assist

    def __init__(self, sim, filename=None):
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected_ 
        clibassist.assist_initialize(byref(sim), byref(self))
        self.process_messages()

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibassist.assist_free_pointers(byref(self))

    def detach(self, sim):
        sim._extras_ref = None # remove reference to assist so it can be garbage collected 
        clibassist.assist_detach(byref(sim), byref(self))
    
# Need to put fields after class definition because of self-referencing
Extras._fields_ =  [("_sim", POINTER(rebound.Simulation)),
        ]

