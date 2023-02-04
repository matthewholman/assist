import rebound
import assist
import unittest
import math

class TestForces(unittest.TestCase):
    def test_errors(self):
        sim = rebound.Simulation()
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        extras = assist.Extras(sim, ephem)
        
        with self.assertRaises(AttributeError) as context:
            extras.forces = "no array"
        
        with self.assertRaises(AttributeError) as context:
            extras.forces = ["Magic"]
        
        with self.assertRaises(AttributeError) as context:
            extras.forces = [1,2,3]


    def test_forces(self):
        
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

        forces = extras.forces
        Nforces = len(forces)
        forces.remove("GR_EIH")
       
        extras.forces = forces
        
        self.assertEqual(Nforces-1, len(extras.forces)) 


        sim.integrate(t+60.)
        sim2.integrate(t+60.)


        d = sim.particles[0] - sim2.particles[0]
        au2meter = 149597870700
    
        self.assertLess(math.fabs(d.x*au2meter-100.), 20) # Expect about a 100m difference with GR off



if __name__ == '__main__':
    unittest.main()
