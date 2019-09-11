import pandas as pd
from scgenome import cncluster
import os
import time

OUT_DIR = "/Users/massoudmaher/data/test_do_bhc/"
CN_DATA_FP = "/Users/massoudmaher/data/sc1935_cn_data_qc.csv"
N_CELLS = 100
# bhc params
N_STATES = 12
ALPHA = 0.3
PROB_CN_CHANGE = 0.8

cn_data = pd.read_csv(CN_DATA_FP)
cn_data = cn_data.iloc[:, 1:]

if N_CELLS is not None:
    keep_cells = cn_data["cell_id"].unique()

n_cells = len(set(cn_data["cell_id"].unique()))
print(f"Doing BHC on {n_cells} cells")

start = time.time()
bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=N_STATES, alpha=ALPHA,
                               prob_cn_change=PROB_CN_CHANGE)
)
print(f"{time.time()-start}s for BHC on {n_cells} cells")

bhc_linkage.to_csv(os.path.join(OUT_DIR, "linkage.csv"))
bhc_cell_ids.to_csv(os.path.join(OUT_DIR, "cell_ids.csv"))
matrix_data.to_csv(os.path.join(OUT_DIR, "matrix_data.csv"))
measurement.to_csv(os.path.join(OUT_DIR, "measurement.csv"))
variances.to_csv(os.path.join(OUT_DIR, "variances.csv"))



