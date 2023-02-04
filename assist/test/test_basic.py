import rebound
import assist
import unittest
import math

class TestAssist(unittest.TestCase):
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
    
    def test_ephem_names(self):
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        p1 = ephem.get_particle(0,0) # Sun
        p2 = ephem.get_particle("Sun",0) # Also Sun
        self.assertEqual(p1.x, p2.x)
        
        with self.assertRaises(ValueError) as context:
            ephem.get_particle("Planet 9",0) # Does not exist

    def test_holman(self):
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        
        sim = rebound.Simulation()
        extras = assist.Extras(sim, ephem)

        sim.t = 8416.5
        
        # Initial conditions
        sim.add(x= -2.724183384883979E+00, 
                y= -3.523994546329214E-02, 
                z= 9.036596202793466E-02,
                vx= -1.374545432301129E-04,
                vy= -1.027075301472321E-02,
                vz= -4.195690627695180E-03)

        sim.integrate(8446.5)

        # Final result from Horizons
        pf = rebound.Particle(x= -2.710320457933958E+00, 
                y= -3.424507930535848E-01, 
                z= -3.582442972611413E-02,
                vx= 1.059255302926290E-03, 
                vy= -1.018748422976772E-02, 
                vz= -4.207712906489264E-03)
    
        # difference
        d = sim.particles[0] - pf 
        
        
        au2meter = 149597870700
        
        self.assertLess(math.fabs(d.x*au2meter), 0.1) # 10cm accurary 
        self.assertLess(math.fabs(d.y*au2meter), 0.1) # 10cm accurary 
        self.assertLess(math.fabs(d.z*au2meter), 0.1) # 10cm accurary 
        self.assertLess(math.fabs(d.vx*au2meter), 5e-3) # 5mm/day accurary 
        self.assertLess(math.fabs(d.vy*au2meter), 5e-3) # 5mm/day accurary 
        self.assertLess(math.fabs(d.vz*au2meter), 5e-3) # 5mm/day accurary 



if __name__ == '__main__':
    unittest.main()
