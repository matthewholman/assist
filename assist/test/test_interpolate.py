import rebound
import assist
import unittest
import math

class TestInterpolate(unittest.TestCase):
    def test_interpolate(self):
        
        t = 8416.5
        sim = rebound.Simulation()
        sim.t = t
        sim.add(x= -2.724183384883979E+00, 
                y= -3.523994546329214E-02, 
                z= 9.036596202793466E-02,
                vx= -1.374545432301129E-04,
                vy= -1.027075301472321E-02,
                vz= -4.195690627695180E-03)

        sim2 = sim.copy()
        
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        extras = assist.Extras(sim, ephem)
        extras2 = assist.Extras(sim2, ephem)

        for i in range(10):
            t += 40.0
            sim.integrate(t)
            extras2.integrate_or_interpolate(t)

            d = sim.particles[0] - sim2.particles[0]
            au2meter = 149597870700
        
            self.assertLess(math.fabs(d.x*au2meter), 0.01) # 1cm accurary 
            self.assertLess(math.fabs(d.y*au2meter), 0.01) # 1cm accurary 
            self.assertLess(math.fabs(d.z*au2meter), 0.01) # 1cm accurary 



if __name__ == '__main__':
    unittest.main()
