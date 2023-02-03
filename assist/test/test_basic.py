import rebound
import assist
import unittest

class TestRebx(unittest.TestCase):
    def test_rebound(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1, inc=0.4)
        sim.integrate(1.0)
        self.assertEqual(sim.t,1.0)
    
    def test_ephem(self):
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        self.assertEqual(ephem.jd_ref, 2451545.0)
        p = ephem.get_particle(0,0) # Sun
        self.assertEqual(p.x, -0.007137179161607906)
        p = ephem.get_particle(1,100) # planet
        self.assertEqual(p.x, 0.12906301685045435)
        p = ephem.get_particle(20,200) #asteroid
        self.assertEqual(p.x, -2.62956381075119)
        del ephem



if __name__ == '__main__':
    unittest.main()
