from scgenome import cncluster, jointcnmodels
from scgenome.utils import cn_mat_to_cn_data, cn_mat_as_df
from unittest import TestCase, main
from .constants_test import *
import pandas as pd
from scipy.special import logsumexp


class TestCluster(TestCase):

    def setUp(self):
        print(E2_VAR)
        print("E2_VAR\n")
        print(E2_COPY)
        print("E2_COPY\n")

        e2_pw_ll = np.zeros((E2_COPY.shape[0], E2_COPY.shape[0]))
        for i in range(E2_COPY.shape[0]):
            for j in range(E2_COPY.shape[0]):
                data = E2_COPY[[i] + [j], :]
                var = E2_VAR[[i] + [j], :]
                ll = jointcnmodels.calculate_marginal_ll_simple(
                    data, var, E2_TRANS
                )
                e2_pw_ll[i, j] = ll
        print(e2_pw_ll)
        print("e2_pw_ll\n\n")

        e2_d1 = 2.302585092994046
        e2_pi1 = 0
        e2_d2 = 4.700480365792417
        e2_pi2 = -2.3978952727983707

        ll = [
            jointcnmodels.calculate_marginal_ll_simple(
                E2_COPY[[i], :],
                E2_VAR[[i], :],
                E2_TRANS
            )

            for i in range(E1_CN_MAT.shape[0])
        ]
        print(ll)
        print("single ll\n")

        e2_pw_tree_ll = np.zeros((E2_COPY.shape[0], E2_COPY.shape[0]))
        for i in range(E2_COPY.shape[0]):
            for j in range(E2_COPY.shape[0]):
                e2_pw_tree_ll[i,j] = logsumexp(
                    [e2_pi2 + e2_pw_ll[i, j],
                     np.log(1-np.exp(e2_pi2)) + ll[i] + ll[j]]
                )
        print(e2_pw_tree_ll)
        print("e2_pw_tree_ll\n\n")

        e2_r1 = e2_pi2 + e2_pw_ll - e2_pw_tree_ll
        print(e2_r1)
        print("e2_r1\n\n")


    def test_bayesian_cluster(self):
        df_cn_mat = cn_mat_as_df(E2_CN_MAT, E2_CHR_NAMES)
        cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=E2_CELL_IDS)
        cn_data["cluster_id"] = (
            cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
        np.random.seed(E2_SEED)
        cn_data["copy2"] = (cn_data["copy"] +
            abs(np.random.normal(scale=0.3, size=cn_data.shape[0]))
        )
        cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                           "cluster_id", "copy"]
        cn_data["reads"] = 1

        linkage, root, cl_cell_ids, matrix_data, measurement, variances = (
            cncluster.bayesian_cluster(cn_data, n_states=11, debug=True,
                                       value_ids=["copy", "state"],
                                       clustering_id="copy",
                                       alpha=E2_ALPHA,
                                       prob_cn_change=E2_TRANS["e0"])
        )

        np.testing.assert_almost_equal(variances, E2_VAR)

        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(cn_data)
            print("cn_data\n")
            print(linkage)
            print("linkage\n")


if __name__ == "__main__":
    main()
