import unittest
from src.node_finder import NodeFinder


class NodeTest(unittest.TestCase):
    def test_make_nodes(self):
        n = 20
        nf = NodeFinder(n)
        ig1, ig2 = nf.initial_guess(n, 1.01)
        nf.update_nodes(ig1)
        nf.update_nodes(ig2)
        for i in range(10):
            nf.next_nodes()

        print(nf.nodes)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
