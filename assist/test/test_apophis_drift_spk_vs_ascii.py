import math
import unittest
from _data_paths import data_path

import numpy as np
import rebound

import assist

class TestApophisDriftSPKvsASCII(unittest.TestCase):
    """
    Compare *relative drift* between the SPK (.bsp) and ASCII-derived binary
    (.440/.441) planets backends for the same Apophis setup as the existing
    Apophis tests.

    Motivation:
    - The existing tests only check a single final-epoch error vs a reference.
    - Here we sample the SPK-vs-.440 residual over time to make regressions
      (e.g. missing asteroid masses in the .440 path) obvious as a growing drift.
    """

    def test_apophis_drift_over_time_spk_vs_440(self):
        au2meter = 149597870700.0

        ephem_spk = assist.Ephem(data_path("de440.bsp"), data_path("sb441-n16.bsp"))
        ephem_ascii = assist.Ephem(data_path("linux_p1550p2650.440"), data_path("sb441-n16.bsp"))

        # Epochs in JD, mirroring the existing Apophis tests.
        t0_jd = 2.4621385359989386e06
        t1_jd = 2.4625030372426095e06

        t0 = t0_jd - ephem_spk.jd_ref
        t1 = t1_jd - ephem_spk.jd_ref

        # Start both simulations from the same barycentric initial condition,
        # using the SPK Sun shift (the Sun shift differences are ~nanometers here).
        apophis_initial_helio = rebound.Particle(
            x=-5.5946538550488512e-01,
            y=8.5647564757574512e-01,
            z=3.0415066217102493e-01,
            vx=-1.3818324735921638e-02,
            vy=-6.0088275597939191e-03,
            vz=-2.5805044631309632e-03,
        )
        sun0 = ephem_spk.get_particle("sun", t0)
        apophis_initial = apophis_initial_helio + sun0

        ng_params = np.array([4.999999873689e-13, -2.901085508711e-14, 0.0])

        def integrate_to(ephem: assist.Ephem, t: float) -> rebound.Particle:
            sim = rebound.Simulation()
            sim.add(apophis_initial)
            sim.t = t0
            sim.ri_ias15.min_dt = 0.001

            extras = assist.Extras(sim, ephem)
            extras.gr_eih_sources = 11
            extras.particle_params = ng_params

            # One-shot integration to `t` (we intentionally do not step through
            # intermediate sample times, because adaptive stepping + output-time
            # constraints can change the final state at the ~100 m scale here).
            sim.integrate(t)
            return sim.particles[0]

        # A small set of checkpoints that captures:
        # - early near-zero drift
        # - substantial accumulated drift by ~1 year
        # - sign changes in dy/dz (residual direction rotates over time)
        checkpoints = [
            ("t0", t0),
            ("t0+120d", t0 + 120.0),
            ("t0+240d", t0 + 240.0),
            ("t0+330d", t0 + 330.0),
            ("t1", t1),
        ]

        residuals = {}
        for label, t in checkpoints:
            p_spk = integrate_to(ephem_spk, t)
            p_ascii = integrate_to(ephem_ascii, t)

            dx_m = (p_spk.x - p_ascii.x) * au2meter
            dy_m = (p_spk.y - p_ascii.y) * au2meter
            dz_m = (p_spk.z - p_ascii.z) * au2meter
            r_m = math.sqrt(dx_m * dx_m + dy_m * dy_m + dz_m * dz_m)

            residuals[label] = (dx_m, dy_m, dz_m, r_m)

        # Basic sanity: using the exact same initial particle, drift at t0 should be ~0.
        _, _, _, r0 = residuals["t0"]
        self.assertLess(r0, 1e-6)

        # Why there is any drift at all:
        # Even when both backends expose the same constants and mass model, SPK and ASCII-derived
        # binary (.440/.441) are not bit-identical sources. They differ slightly at the coefficient
        # / interpolation level (and they can differ in how floating-point operations are ordered
        # in the code). With an adaptive integrator, those tiny differences can cause different
        # step sequences and roundoff accumulation, which shows up as a growing but still small
        # positional residual over time.
        #
        # This test is meant to catch *big* regressions (missing asteroid masses, wrong backend,
        # bad constants), not enforce bit-exact trajectories.
        _, _, _, r_final = residuals["t1"]
        self.assertLess(r_final, 5000.0)  # 5 km upper bound

        # Ensure the drift isn't accidentally *exactly* zero (e.g. wrong backend used).
        self.assertGreater(r_final, 1e-3)  # 1 mm

        # \"Cyclical\" behavior in the *direction* of the residual vector:
        # dy and dz should take on both positive and negative values over the window.
        #
        # Use an epsilon so near-zero doesn't trip a sign test.
        eps = 1e-3  # 1 mm
        dy_vals = [residuals[k][1] for k, _t in checkpoints if k != "t0"]
        dz_vals = [residuals[k][2] for k, _t in checkpoints if k != "t0"]
        self.assertTrue(any(v > eps for v in dy_vals))
        self.assertTrue(any(v < -eps for v in dy_vals))
        self.assertTrue(any(v > eps for v in dz_vals))
        self.assertTrue(any(v < -eps for v in dz_vals))


if __name__ == "__main__":
    unittest.main()


