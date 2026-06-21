
from ._libassist import clibassist
from ctypes import POINTER, byref, c_int, c_double
import rebound
    
def simulation_convert_to_rebound(sim, ephem, merge_moon=1):
    clibassist.assist_simulation_convert_to_rebound.restype = POINTER(rebound.Simulation)
    sim2 = clibassist.assist_simulation_convert_to_rebound(byref(sim), byref(ephem), c_int(merge_moon))
    return sim2.contents

def assist_create_interpolated_simulation(sa, t):
    clibassist.assist_create_interpolated_simulation.restype = POINTER(rebound.Simulation)
    sim = clibassist.assist_create_interpolated_simulation(byref(sa), c_double(t))
    return sim.contents
