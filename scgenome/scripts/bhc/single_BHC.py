"""
Runs bayesian hierarchical clustering on a simulation of poisson random walk
"""
import scgenome.simulation as sim
import numpy as np
from scgenome import cncluster
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df

NUM_SAMPLE = 20
NUM_BIN = 6000
CHR_NAMES = ["1", "2"]
max_cn = 7

cn_mat = sim.cn_mat_poisson(NUM_SAMPLE, NUM_BIN, init_lambda=1., jump_lambda=1,
                            seed=None, max_cn=max_cn)
cell_ids = [f"cl1_cell{i}" for i in range(cn_mat.shape[0])]

CHR_NAMES = ["1", "2"]
df_cn_mat = cn_mat_as_df(cn_mat, CHR_NAMES)
cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=cell_ids)
cn_data["cluster_id"] = (
    cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
cn_data["copy2"] = (cn_data["copy"] +
                    np.absolute(np.random.normal(size=cn_data.shape[0],
                                                 scale=0.3)))
cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                   "cluster_id", "copy"]

tlinkage, root, cl_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=max_cn, value_ids=["copy"]))
