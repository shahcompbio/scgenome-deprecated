from unittest import TestCase, main
from scgenome.TNode import TNode


class TestTNode(TestCase):

    def test_referencing(self):
        """Test if we can make tree from TNode objects"""
        #(sample_inds, left_child, right_child, pi, d, ll, r):
        n1 = TNode([0, 1], None, None, -1)
        n2 = TNode([10, 20], None, None, -1)
        n3 = TNode([30, 40], None, None, -1)
        n1.left_child = n2
        n1.right_child = n3

        clusters = [n1, n2, n3]
        del clusters[1]
        self.assertEqual(TNode([0, 1], n2, n3, -1), clusters[0])

    def test_get_leaves(self):
        n   = TNode([0, 1], None, None, -1)
        n_l = TNode([10, 20], None, None, -1)
        n_r = TNode([30, 40], None, None, -1)
        n_l_l = TNode([50, 60], None, None, -1)
        n_r_l = TNode([70, 80], None, None, -1)
        n_r_r = TNode([80, 90], None, None, -1)

        n.left_child = n_l
        n.right_child = n_r
        n_l.left_child = n_l_l
        n_r.left_child = n_r_l
        n_r.right_child = n_r_r

        exp_n_l_leaves = [n_l_l]
        exp_n_r_leaves = [n_r_l, n_r_r]
        exp_n_leaves = exp_n_l_leaves + exp_n_r_leaves

        obs = n_l.get_leaves()
        self.assertListEqual(obs, exp_n_l_leaves)

        obs = n_r.get_leaves()
        self.assertListEqual(obs, exp_n_r_leaves)

        obs = n.get_leaves()
        self.assertListEqual(obs, exp_n_leaves)


