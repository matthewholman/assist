import math
import unittest
from _data_paths import data_path

import numpy as np
import rebound

import assist

class TestAssistASCIIApophis(unittest.TestCase):
    def test_apophis_ascii_440(self):
        """
        Apophis integration using the ASCII-derived binary DE440 planets file
        (linux_p1550p2650.440). This mirrors the SPK-based Apophis test but
        exercises the ASCII (.440/.441) pathway with its own tolerance.
        """
        ephem = assist.Ephem(data_path("linux_p1550p2650.440"), data_path("sb441-n16.bsp"))

        t_initial = 2.4621385359989386e06 - ephem.jd_ref  # Julian Days relative to jd_ref

        apophis_initial_helio = rebound.Particle(
            x=-5.5946538550488512e-01,  # AU
            y=8.5647564757574512e-01,
            z=3.0415066217102493e-01,
            vx=-1.3818324735921638e-02,  # AU/day
            vy=-6.0088275597939191e-03,
            vz=-2.5805044631309632e-03,
        )

        sun_initial = ephem.get_particle("sun", t_initial)
        apophis_initial = apophis_initial_helio + sun_initial

        sim = rebound.Simulation()
        sim.add(apophis_initial)
        sim.t = t_initial
        sim.ri_ias15.min_dt = 0.001
        extras = assist.Extras(sim, ephem)
        extras.gr_eih_sources = 11

        # Non-gravitational parameters
        extras.particle_params = np.array([4.999999873689e-13, -2.901085508711e-14, 0.0])

        t_final = 2.4625030372426095e06 - ephem.jd_ref
        sim.integrate(t_final)

        # Reference final state (from JPL/NASA reference integration), converted to barycentric below.
        apophis_final_ref = rebound.Particle(
            x=1.7028330901729331e-02,
            y=1.2193934090901304e00,
            z=4.7823589236374386e-01,
            vx=-1.3536187639388663e-02,
            vy=5.3200999989786943e-04,
            vz=-1.6648346717629861e-05,
        )
        # Convert from heliocentric to barycentric
        apophis_final_ref += ephem.get_particle("sun", t_final)

        delta = sim.particles[0] - apophis_final_ref
        delta_pos = np.sqrt(delta.x ** 2 + delta.y ** 2 + delta.z ** 2)
        delta_pos *= 149597870700  # AU -> meters

        # The ASCII-derived binary (.440/.441) planets path should behave equivalently to the
        # SPK (.bsp) planets path for this integration once asteroid masses are loaded.
        self.assertLess(math.fabs(delta_pos), 300.0)


if __name__ == '__main__':
    unittest.main()







