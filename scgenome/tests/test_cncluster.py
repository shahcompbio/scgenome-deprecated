import numpy as np
from unittest import TestCase, main
from scgenome.cncluster import bayesian_cluster


class TestCluster(TestCase):

    def test_bayesian_cluster(self):
        num_bin = 100
        mean = [1, 1, -1, -1]
        np.random.seed(1)
        data = [np.random.normal(x, size=num_bin)[np.newaxis, :] for x in mean]
        mat = np.concatenate(data, axis=0)

        cell_ids = ["cl1_cell1", "cl1_cell2", "cl2_cell1", "cl2_cell2"]

        linkage, root, cell_ids = (
            bayesian_cluster(cn_mat=mat, cell_ids=cell_ids))

        print(linkage)
        print(root)
        print(cell_ids)


