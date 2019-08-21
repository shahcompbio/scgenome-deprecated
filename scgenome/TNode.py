from math import gamma

from scgenome.jointcnmodels import calculate_marginal_ll_simple
from .constants import ALPHA, NO_CHILDREN


class TNode:

    # TODO assert types
    def __init__(self, sample_indeces, left_child, right_child, pi, d, ll, r):
        self.sample_indeces = sample_indeces
        self.left_child = left_child
        self.right_child = right_child
        self.pi = pi
        self.d = d
        self.ll = ll
        self.r = r

    def is_leaf(self):
        return (self.left_child is None and
                self.right_child is None and
                len(self.sample_indeces) == 1)

    def calculate_pi_d(self, alpha=ALPHA):
        if self.left_child is None or self.right_child is None:
            raise ValueError(NO_CHILDREN)

        n_k = len(self.left_child.sample_indeces) + len(
            self.right_child.sample_indeces)
        gnk = gamma(n_k)
        self.d = alpha * gnk + self.left_child.d * self.right_child.d
        self.pi = alpha * gnk / self.d

        return self.pi, self.d

    def calculate_rk(self):
        # TODO make abstract class/interface with only this and calculate ll
        # And this is DPMM implementation
        pass

    def get_ll(self, measurement, variances, tr_mat):
        self.ll = calculate_marginal_ll_simple(
            measurement[self.sample_indeces, :],
            variances[self.sample_indeces, :],
            tr_mat
        )
        # TODO feels weird in case where there is only 1 sample index
        return self.ll

    def get_r(self):
        top = self.pi * self.ll
        bottom = (self.pi * self.ll +
            (1 - self.pi) * self.left_child.ll * self.right_child.ll
        )
        self.r = top / bottom
        return self.r

#    def calculate_pi_d(self):
#        n_k = len(self.left_child.sample_indeces) +
#              len(self.right_child.sample_indeces)
#
#        if self.left_child is None and self.right_child is None:
#            # Node is a leaf
#            self.d = self.alpha
#            self.pi = 1
#        else:
#            left_pi, left_d = calculate_pi(left_cluster.left_child,
#                                           left_cluster.right_child, alpha)
#            right_pi, right_d = calculate_pi(right_cluster.right_child,
#                                             right_cluster.right_child, alpha)
#
#            gnk = gamma(n_k)
#            d = alpha*gnk + left_d*right_d
#            pi = alpha*gnk / d
#
#        return pi, d
