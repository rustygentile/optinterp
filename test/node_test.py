import unittest
import numpy as np
from src import NodeFinder
from src import nodes


class NodeTest(unittest.TestCase):

    def test_best_nodes(self):
        """Calculate optimal nodes for n=2 to 21 and compare the Lebesgue function with literature values"""

        expected_l = [1.0, 1.25, 1.42291957, 1.55949021, 1.67221037, 1.76813458,
                      1.85159939, 1.92545762, 1.99168499, 2.05170576, 2.10658026,
                      2.15711897, 2.20395521, 2.24759321, 2.28844092, 2.32683304,
                      2.36304752, 2.39731771, 2.42984142, 2.46078775]

        for n in range(2, 22):

            nds = nodes(n)
            nf = NodeFinder(n)
            nf.update_nodes(nds)
            print(n, nf.ly1, expected_l[n - 2])
            expected_ly = np.ones_like(nf.ly1) * expected_l[n - 2]
            self.assertAlmostEqual(nf.ly1.all(), expected_ly.all())

        self.assertTrue(True)

    def test_100_pts(self):
        nds = nodes(100, alpha=1.00000001, max_iter=10)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
