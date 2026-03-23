import rebound
import assist
import unittest
import math
import numpy as np
from _data_paths import data_path

class TestAssistSPK(unittest.TestCase):
    def test_apophis(self):
        ephem = assist.Ephem(data_path("de440.bsp"), data_path("sb441-n16.bsp"))
        t_initial = 2.4621385359989386E+06 - ephem.jd_ref
        apophis_initial_helio = rebound.Particle(x = -5.5946538550488512E-01,
                                                 y =  8.5647564757574512E-01,
                                                 z =  3.0415066217102493E-01,
                                                 vx= -1.3818324735921638E-02,
                                                 vy= -6.0088275597939191E-03,
                                                 vz= -2.5805044631309632E-03)
        sun_initial = ephem.get_particle("sun", t_initial)
        apophis_initial = apophis_initial_helio + sun_initial
        sim = rebound.Simulation()
        sim.add(apophis_initial)
        sim.t = t_initial
        sim.ri_ias15.min_dt = 0.001
        extras = assist.Extras(sim, ephem)
        extras.gr_eih_sources = 11
        extras.particle_params = np.array([4.999999873689E-13, -2.901085508711E-14, 0.0])
        t_final = 2.4625030372426095E+06  - ephem.jd_ref
        sim.integrate(t_final)
        # Reference final state (from JPL/NASA reference integration), converted to barycentric below.
        apophis_final_ref = rebound.Particle(x =  1.7028330901729331E-02,
                                             y =  1.2193934090901304E+00,
                                             z =  4.7823589236374386E-01,
                                             vx= -1.3536187639388663E-02,
                                             vy=  5.3200999989786943E-04,
                                             vz= -1.6648346717629861E-05)
        apophis_final_ref += ephem.get_particle("sun", t_final)
        delta = sim.particles[0] - apophis_final_ref
        delta_pos = np.sqrt(delta.x**2 + delta.y**2 + delta.z**2)
        delta_pos *= 149597870700
        self.assertLess(math.fabs(delta_pos), 300)

if __name__ == '__main__':
    unittest.main()



