from scgenome import cncluster
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df
from unittest import TestCase, main
from .constants_test import *


class TestCluster(TestCase):

    def test_bayesian_cluster(self):
        df_cn_mat = cn_mat_as_df(E1_CN_MAT, E1_CHR_NAMES)
        cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=E1_CELL_IDS)
        cn_data["cluster_id"] = (
            cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
        cn_data["copy2"] = (cn_data["copy"] +
            np.absolute(np.random.normal(size=cn_data.shape[0], scale=0.3)))
        cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                           "cluster_id", "copy"]

        linkage, root, cl_cell_ids = (
            cncluster.bayesian_cluster(cn_data, n_states=11,
                                       value_ids=["copy", "state"]))

        print(linkage)
        print(root)
        print(cl_cell_ids)


if __name__ == "__main__":
    main()
