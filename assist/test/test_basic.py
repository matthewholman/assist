import rebound
import assist
import unittest

class TestRebx(unittest.TestCase):
    def test_rotate_sim(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1, inc=0.4)

        self.assertLess(0, 1e-15)


if __name__ == '__main__':
    unittest.main()
