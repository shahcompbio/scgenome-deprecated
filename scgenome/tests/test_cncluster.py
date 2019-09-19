from scgenome import cncluster, jointcnmodels
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df
from unittest import TestCase, main
from .constants_test import *
import pandas as pd


class TestCluster(TestCase):

    def test_bayesian_cluster(self):
        df_cn_mat = cn_mat_as_df(E1_CN_MAT, E1_CHR_NAMES)
        cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=E1_CELL_IDS)
        cn_data["cluster_id"] = (
            cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
        np.random.seed(2)
        cn_data["copy2"] = cn_data["copy"] + abs(np.random.normal(scale=0.3, size=cn_data.shape[0]))
        cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                           "cluster_id", "copy"]

        linkage, root, cl_cell_ids, matrix_data, measurement, variances = (
            cncluster.bayesian_cluster(cn_data, n_states=11, debug=True,
                                       value_ids=["copy", "state"]))

        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(linkage)
        print(root)
        print(cl_cell_ids)
        print("expected----------")
        prob_same = 0.8
        transmodel = {"kind": "twoparam", "e0": prob_same,
                      "e1": 1 - prob_same}
        ll = [
            jointcnmodels.calculate_marginal_ll_simple(
                E1_CN_MAT[[i], :],
                variances[[i], :],
                transmodel
            )

            for i in range(E1_CN_MAT.shape[0])
        ]
        print(f"ll {ll}")
        print(f"variances\n{variances}")
        print(f"cn_data:\n{cn_data}")

    def test_bayesian_cluster3(self):
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
        print(root)
        print(cl_cell_ids)
        print("expected----------")
        prob_same = 0.8
        transmodel = {"kind": "twoparam", "e0": prob_same,
                      "e1": 1 - prob_same}
        ll = [
            jointcnmodels.calculate_marginal_ll_simple(
                E2_CN_MAT[[i], :],
                variances[[i], :],
                transmodel
            )

            for i in range(E1_CN_MAT.shape[0])
        ]
        np.testing.assert_almost_equal(variances, E2_VAR)
        print(f"ll {ll}")
        print(f"variances\n{variances}")
        print(f"cn_data:\n{cn_data}")
        print("*********************************-")
        ll_01 = jointcnmodels.calculate_marginal_ll_simple(E2_CN_MAT[[0,1], :], variances[[0,1], :], transmodel)
        pi_01 =  1 / (1+alpha**2)
        tll_01 = 1
        print(f"ll_01: {jointcnmodels.calculate_marginal_ll_simple(E2_CN_MAT[[0,1], :], variances[[0,1], :], transmodel)}")
        print(f"ll_01: {jointcnmodels.calculate_marginal_ll_simple(E2_CN_MAT[[0,1], :], variances[[0,1], :], transmodel)}")



if __name__ == "__main__":
    main()
