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

    def get_pi_d(self, alpha=ALPHA):
        if self.left_child is None or self.right_child is None:
            raise ValueError(NO_CHILDREN)

        n_k = len(self.left_child.sample_indeces) + len(
            self.right_child.sample_indeces)
        gnk = gamma(n_k)
        self.d = alpha * gnk + self.left_child.d * self.right_child.d
        self.pi = alpha * gnk / self.d

        return self.pi, self.d

    # TODO is weird in the case where there is only 1 sample index
    def get_ll(self, measurement, variances, tr_mat):
        self.ll = calculate_marginal_ll_simple(
            measurement[self.sample_indeces, :],
            variances[self.sample_indeces, :],
            tr_mat
        )
        return self.ll

    def get_r(self):
        top = self.pi * self.ll
        bottom = (self.pi * self.ll +
            (1 - self.pi) * self.left_child.ll * self.right_child.ll
        )
        self.r = top / bottom
        return self.r

    def __str__(self):
        return f"sample_indeces : {self.sample_indeces}, " \
            f"left_child : {self.left_child}, " \
            f"right_child : {self.right_child}, " \
            f"pi : {self.pi}, " \
            f"d : {self.d}, " \
            f"ll : {self.ll}, " \
            f"r : {self.r}"

    def __eq__(self, other):
        return (
            self.sample_indeces == other.sample_indeces and
            self.left_child == other.left_child and
            self.right_child == other.right_child and
            self.pi == other.pi and
            self.d == other.d and
            self.ll == other.ll and
            self.r == other.r
        )
