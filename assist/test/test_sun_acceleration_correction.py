import csv
import ctypes
import json
import math
import unittest
from pathlib import Path

import assist
import rebound
from _data_paths import data_path


class TestSunAccelerationCorrection(unittest.TestCase):
    def test_holman_against_cached_horizons(self):
        fixture = Path(__file__).with_name("fixtures") / "holman-horizons.json"
        response = json.loads(fixture.read_text())["result"]
        table = response.split("$$SOE")[1].split("$$EOE")[0]
        states = list(csv.reader(table.strip().splitlines()))
        epoch = float(states[0][0])
        initial = dict(zip(("x", "y", "z", "vx", "vy", "vz"), map(float, states[0][2:8])))
        ephem = assist.Ephem(data_path("linux_p1550p2650.440"), data_path("sb441-n16.bsp"))

        # Independent integrations from the same geometric barycentric state.
        self.assertEqual(len(states), 5)
        for row, baseline, corrected_limit in zip(states[1:], (0.01132, 1.044, 244.15, 6540.94),
                                                   (0.003, 0.05, 12.0, 300.0)):
            for enabled in (False, True):
                with self.subTest(jd=row[0], correction=enabled):
                    sim = rebound.Simulation()
                    sim.t = epoch - ephem.jd_ref
                    extras = assist.Extras(sim, ephem)
                    self.assertNotIn("SUN_ACCELERATION_CORRECTION", extras.forces)
                    if enabled:
                        extras.forces += ["SUN_ACCELERATION_CORRECTION"]
                    sim.add(**initial)
                    sim.integrate(float(row[0]) - ephem.jd_ref)
                    error_m = math.dist(sim.particles[0].xyz, tuple(map(float, row[2:5]))) * 149597870700
                    if enabled:
                        self.assertLess(error_m, corrected_limit)
                    else:
                        self.assertAlmostEqual(error_m, baseline, delta=max(0.0005, baseline * 0.02))

    def test_acceleration_matches_ephemeris_difference(self):
        ephem = assist.Ephem(data_path("linux_p1550p2650.440"), data_path("sb441-n16.bsp"))
        for t in (7304.5, 8304.5):
            sim = rebound.Simulation()
            sim.t = t
            extras = assist.Extras(sim, ephem)
            extras.forces = ["SUN_ACCELERATION_CORRECTION"]
            sim.add(x=1)
            sim.add(y=2)
            sim.add_variation()
            sun, *bodies = (ephem.get_particle(i, t) for i in range(27))
            expected = [sun.ax, sun.ay, sun.az]
            for body in bodies:
                offset = [b - s for b, s in zip(body.xyz, sun.xyz)]
                scale = body.m / math.hypot(*offset)**3
                for axis, distance in enumerate(offset):
                    expected[axis] -= scale * distance
            assist.clibassist.assist_additional_forces(ctypes.byref(sim))
            for i, particle in enumerate(sim.particles):
                acceleration = (particle.ax, particle.ay, particle.az)
                if i < 2:
                    self.assertLess(math.dist(acceleration, expected), 1e-21)
                else:
                    self.assertEqual(acceleration, (0., 0., 0.))
