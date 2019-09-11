from math import gamma

from scgenome.jointcnmodels import calculate_marginal_ll_simple
from .constants import ALPHA, NO_CHILDREN
from scipy.special import logsumexp
import numpy as np


class TNode:

    # TODO assert types
    def __init__(self, sample_inds, left_child, right_child, cluster_ind,
                 pi=None, d=None, ll=None, tree_ll=None, log_r=None):

        self.sample_inds = sample_inds
        self.left_child = left_child
        self.right_child = right_child
        self.cluster_ind = cluster_ind

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
        self.update_tree_ll(measurement, variances, transmodel)
        self.update_log_r()

    #def update_pi_d(self, alpha=ALPHA):
    #    if self.is_leaf():
    #        self.pi = 1
    #        self.d = alpha
    #    else:
    #        n_k = (len(self.left_child.sample_inds) +
    #               len(self.right_child.sample_inds))
    #        gnk = gamma(n_k)
    #        self.d = alpha * gnk + self.left_child.d * self.right_child.d
    #        self.pi = alpha * gnk / self.d
    def update_pi_d(self, alpha=ALPHA):
        log_alpha = np.log(alpha)
        if self.is_leaf():
            self.pi = np.log(1)
            self.d = log_alpha
        else:
            n_k = (len(self.left_child.sample_inds) +
                   len(self.right_child.sample_inds))

            log_gnk = np.log(gamma(n_k))
            self.d = logsumexp([log_alpha + log_gnk,
                                self.left_child.d + self.right_child.d])
            self.pi = log_alpha + log_gnk - self.d

    def update_tree_ll(self, measurement, variances, transmodel):
        # TODO this should be named update, not get
        # TODO refactor so these are called in init and theres base cases
        if len(self.sample_inds) == 1:
            self.tree_ll = calculate_marginal_ll_simple(
                measurement[self.sample_inds, :],
                variances[self.sample_inds, :], transmodel)
        else:
            #first = np.log(self.pi) + self.ll
            first = self.pi + self.ll
            #second = (np.log(1 - self.pi) + self.left_child.tree_ll +
            second = (np.log(1) - self.pi + self.left_child.tree_ll +
                      self.right_child.tree_ll)
            self.tree_ll = logsumexp([first, second])

    def update_ll(self, measurement, variances, transmodel):
        if self.is_leaf():
            self.ll = 1
        else:
            self.ll = calculate_marginal_ll_simple(
                measurement[self.sample_inds, :],
                variances[self.sample_inds, :], transmodel)

    def update_log_r(self):
        #top = np.log(self.pi) + self.ll
        top = self.pi + self.ll
        #bottom = logsumexp([np.log(self.pi) + self.ll,
        #                    np.log(1 - self.pi) +
        #                    self.left_child.tree_ll +
        #                    self.right_child.tree_ll])
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

