from unittest import TestCase, main
from scgenome.TNode import TNode


class TestTNode(TestCase):

    def test_referencing(self):
        """Test if we can make tree from TNode objects"""
        #(sample_inds, left_child, right_child, pi, d, ll, r):
        n1 = TNode([0, 1],   None, None, 1, 2,  3, 4)
        n2 = TNode([10, 20], None, None, 5, 6,  7, 8)
        n3 = TNode([30, 40], None, None, 9, 10, 11, 12)
        n1.left_child = n2
        n1.right_child = n3

        clusters = [n1, n2, n3]
        del clusters[1]
        self.assertEqual(TNode([0, 1], n2, n3, 1, 2,  3, 4), clusters[0])
