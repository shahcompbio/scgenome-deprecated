import scgenome.simulation as sim
import numpy as np
import pandas as pd
from IPython.display import display
from scgenome import cncluster
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from unittest import TestCase, main
from scgenome.TNode import TNode


class TestCluster(TestCase):

    def test_bayesian_cluster(self):
        num_sample = 4
        num_bin = 500
        chr_names = ["1", "2"]
        max_cn = 7

        cluster1 = sim.cn_mat_poisson(num_sample, num_bin, init_lambda=1.,
                                      jump_lambda=1, seed=1, max_cn=max_cn)
        cluster2 = sim.cn_mat_poisson(num_sample, num_bin, init_lambda=3.,
                                      jump_lambda=0.1, seed=3, max_cn=max_cn)

        clst1_cell_ids = [f"cl1_cell{i}" for i in range(cluster1.shape[0])]
        clst2_cell_ids = [f"cl2_cell{i}" for i in range(cluster1.shape[0])]

        cn_mat = np.concatenate([cluster1, cluster2])
        cell_ids = clst1_cell_ids + clst2_cell_ids

        df_cn_mat = cn_mat_as_df(cn_mat, chr_names)
        cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=cell_ids)
        cn_data["cluster_id"] = cn_data["cell_id"].str.split("_",expand=True).iloc[:,0]
        cn_data["copy2"] = cn_data["copy"] + np.absolute(np.random.normal(size=cn_data.shape[0], scale=0.3))
        cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                           "cluster_id", "copy"]

        tlinkage, root, cl_cell_ids = (
            cncluster.rec_bayesian_cluster(cn_data, n_states=max_cn,
            value_ids=["copy"]))

        print(tlinkage)

