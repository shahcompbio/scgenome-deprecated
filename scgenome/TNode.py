from scipy.special import loggamma

from scgenome.jointcnmodels import calculate_marginal_ll_simple
from .constants import ALPHA, NO_CHILDREN
from scipy.special import logsumexp
import numpy as np


class TNode:

    # TODO assert types
    def __init__(self, sample_inds, left_child, right_child, cluster_ind,
                 pi=None, d=None, ll=None, tree_ll=None, log_r=None):

        # Indeces of sampels in this cluster
        self.sample_inds = sample_inds
        self.left_child = left_child
        self.right_child = right_child
        # Cluster index, maps to i/j in linkage matrix
        self.cluster_ind = cluster_ind

        # pi and d are actally log(pi), log(d)
        self.pi = pi
        self.d = d
        self.ll = ll
        self.tree_ll = tree_ll
        self.log_r = log_r

    def is_leaf(self):
        return (self.left_child is None and
                self.right_child is None and
                len(self.sample_inds) == 1)

    def update_vars(self, measurement, variances, transmodel, alpha=ALPHA):
        self.update_pi_d(alpha)
        self.update_ll(measurement, variances, transmodel)
        self.update_tree_ll()
        self.update_log_r()

    def update_pi_d(self, alpha=ALPHA):
        log_alpha = np.log(alpha)
        if self.is_leaf():
            self.pi = np.log(1)
            self.d = log_alpha
        else:
            n_k = (len(self.left_child.sample_inds) +
                   len(self.right_child.sample_inds))

            log_gnk = loggamma(n_k) # Need this or else gamma(n_k) too large
            self.d = logsumexp([log_alpha + log_gnk,
                                self.left_child.d + self.right_child.d])
            self.pi = log_alpha + log_gnk - self.d

    def update_tree_ll(self):
        if self.is_leaf():
            self.tree_ll = self.ll
        else:
            first = self.pi + self.ll
            # pi becomes 0 sometimes at later iterations of bhc. This
            # prevents pmo from becoming infinity
            pmo = 0 if self.pi == 0 else np.log(1 - np.exp(self.pi))
            second = pmo + self.left_child.tree_ll + self.right_child.tree_ll
            self.tree_ll = logsumexp([first, second])

    def update_ll(self, measurement, variances, transmodel):
        self.ll = calculate_marginal_ll_simple(
            measurement[self.sample_inds, :],
            variances[self.sample_inds, :], transmodel)

    def update_log_r(self):
        top = self.pi + self.ll
        bottom = self.tree_ll
        self.log_r = top - bottom

    def get_leaves(self, leaves=None):
        if leaves is None:
            leaves = []
        if self.left_child is None and self.right_child is None:
            leaves.append(self)
            return leaves
        elif self.right_child is None:
            return self.left_child.get_leaves(leaves)
        elif self.left_child is None:
            return self.right_child.get_leaves(leaves)
        else:
            leaves = leaves + self.left_child.get_leaves()
            leaves = leaves + self.right_child.get_leaves()
            return leaves

    def __str__(self):
        return f"sample_inds : {self.sample_inds}, " \
            f"left_child : {self.left_child.__repr__()} " \
            f"right_child : {self.left_child.__repr__()} " \
            f"pi : {self.pi}, " \
            f"d : {self.d}, " \
            f"ll : {self.ll}, " \
            f"log_r : {self.log_r}, " \
            f"cluster_ind : {self.cluster_ind}"

    def __eq__(self, other):
        return (
                self.sample_inds == other.sample_inds and
                self.left_child == other.left_child and
                self.right_child == other.right_child and
                self.pi == other.pi and
                self.d == other.d and
                self.ll == other.ll and
                self.log_r == other.log_r and
                self.cluster_ind == other.cluster_ind
        )

