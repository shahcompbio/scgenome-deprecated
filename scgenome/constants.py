import numpy as np
# Default params
ALPHA = 0.3
CN_DATA_ID = "copy"
CELL_ID = "cell_id"
MAX_CN = 8
INIT_CN = 2
VALUE_IDS = ["reads", "state", "copy"]
INDEX_IDS = ['chr', 'start', 'cell_id']
COPY_ID = "copy"
SIM_META = {
    "samples_per_cluster": np.int64,
    "num_bin": np.int64,
    "max_cn": np.int64,
    "alpha": np.float64,
    "init_lambdas": object,
    "jump_lambdas": object,
    "cn_data": object,
    "cn_mat": object,
    "plinkage": object,
    "plot_data": object,
    "clustering": object,
    "prop_correct": np.float64,
    "cell_id": str
}
LINKAGE_COLS = ["i", "j", "r_merge", "naive_dist", "log_like",
                "i_count", "j_count"]
LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
CACHE_DIR = "Scratch/tantalus_data"

# Error messages
NO_CHILDREN = "Node has no children to comptue pi, d with"
TRANS_NOT_SYMMETRIC = "Transition matrix must be symmetric"
TRANS_NOT_SUM = "Rows, columns of transition matrix must sum to 1"
TRANS_NOT_SQUARE = "Transition matrix must be square"
NBIN_NCHR = "Number of bins must be divisible by number of chromosomes"
CHR_CLUSTERING = "Number of chromosome clusters exceeds number" \
                 " of transition matrices"

