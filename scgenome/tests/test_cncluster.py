from scgenome import cncluster, jointcnmodels
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df
from unittest import TestCase, main
from .constants_test import *
import pandas as pd


class TestCluster(TestCase):

    def test_bayesian_cluster(self):
        df_cn_mat = cn_mat_as_df(E2_CN_MAT, E2_CHR_NAMES)
        cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=E2_CELL_IDS)
        cn_data["cluster_id"] = (
            cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
        np.random.seed(2)
        cn_data["copy2"] = cn_data["copy"] + abs(np.random.normal(scale=0.3, size=cn_data.shape[0]))
        cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                           "cluster_id", "copy"]
        cn_data["reads"] = 1

        alpha=10
        linkage, root, cl_cell_ids, matrix_data, measurement, variances = (
            cncluster.bayesian_cluster(cn_data, n_states=11, debug=True,
                                       value_ids=["copy", "state"],
                                       alpha=alpha))

        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(linkage)
            print("linkage\n")
            print(variances)
            print("variances\n")



if __name__ == "__main__":
    main()
