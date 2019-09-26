import numpy as np
import pandas as pd

####
# E1_V1 = np.array([1, 2, 3, 4])[np.newaxis, :]
# E1_V2 = E1_V1 + 1
# E1_V3 = np.array([5, 20, 1, 8])[np.newaxis, :]
# E1_V4 = E1_V3 + 2
# E1_CN_MAT = np.concatenate([E1_V1, E1_V2, E1_V3, E1_V4])
E1_V1 = np.array([1, 1, 9, 9])[np.newaxis, :]
E1_V2 = E1_V1 + 1
E1_V3 = np.array([5, 5, 9, 9])[np.newaxis, :]
E1_V4 = E1_V3 + 2
E1_CN_MAT = np.concatenate([E1_V1, E1_V2, E1_V3, E1_V4])

E1_CELL_IDS = ["c1", "c2", "c3", "c4"]
E1_CHR_NAMES = ["1", "2"]

#######
E2_V1 = np.array([1, 1, 9, 9, 9, 9])[np.newaxis, :]
E2_V2 = E2_V1 + 1
E2_V3 = np.array([5, 5, 9, 9, 9, 9])[np.newaxis, :]
E2_V4 = E2_V3 + 2
E2_CN_MAT = np.concatenate([E2_V1, E2_V2, E2_V3, E2_V4])

E2_CELL_IDS = ["c1", "c2", "c3", "c4"]
E2_CHR_NAMES = ["1", "2"]
E2_VAR = [[0.05, 0.00584792, 0.05, 0.05, 0.05, 0.05,
           0.05, 0.05, 0.05, 0.02704375, 0.05],
          [0.05, 0.05, 0.02480254, 0.05, 0.05, 0.05,
           0.05, 0.05, 0.05, 0.05, 0.05155068],
          [0.05, 0.05, 0.05, 0.05, 0.05, 0.05213731,
           0.05, 0.05, 0.05, 0.0201155, 0.05, ],
          [0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
           0.05, 0.02456518, 0.05, 0.05, 0.05]]
E2_VAR = np.array(E2_VAR)

E2_ALPHA = 10

#######
E3_CN_DATA = pd.DataFrame(
    data=[
        ["1", 0, 1, 1, -1, 1, 1, "c1", "s", "l", "o1", 1, 0],
        ["1", 2, 3, 1, -1, 1, 1, "c1", "s", "l", "o1", 1, 0],
        ["1", 4, 5, 1, -1, 1, 1, "c1", "s", "l", "o1", 1, 0],
        ["1", 6, 7, 1, -1, 1, 1, "c1", "s", "l", "o1", 1, 0],
        ["2", 0, 1, 1, -1, 1, 2, "c1", "s", "l", "o1", 1, 0],
        ["2", 2, 3, 1, -1, 1, 2, "c1", "s", "l", "o1", 1, 0],
        ["2", 4, 5, 1, -1, 1, 2, "c1", "s", "l", "o1", 1, 0],
        ["2", 6, 7, 1, -1, 1, 2, "c1", "s", "l", "o1", 1, 0],

        ["1", 0, 1, 2, -2, 2, 2, "c2", "s", "l", "o1", 3, 0],
        ["1", 2, 3, 2, -2, 2, 2, "c2", "s", "l", "o1", 3, 0],
        ["1", 4, 5, 2, -2, 2, 2, "c2", "s", "l", "o1", 3, 0],
        ["1", 6, 7, 2, -2, 2, 2, "c2", "s", "l", "o1", 3, 0],
        ["2", 0, 1, 2, -2, 2, 3, "c2", "s", "l", "o1", 3, 0],
        ["2", 2, 3, 2, -2, 2, 3, "c2", "s", "l", "o1", 3, 0],
        ["2", 4, 5, 2, -2, 2, 3, "c2", "s", "l", "o1", 3, 0],
        ["2", 6, 7, 2, -2, 2, 3, "c2", "s", "l", "o1", 3, 0],

        ["1", 0, 1, 3, -3, 3, 3, "c3", "s", "l", "o2", 2, 1],
        ["1", 2, 3, 3, -3, 3, 3, "c3", "s", "l", "o2", 2, 1],
        ["1", 4, 5, 3, -3, 3, 3, "c3", "s", "l", "o2", 2, 1],
        ["1", 6, 7, 3, -3, 3, 3, "c3", "s", "l", "o2", 2, 1],
        ["2", 0, 1, 3, -3, 3, 4, "c3", "s", "l", "o2", 2, 1],
        ["2", 2, 3, 3, -3, 3, 4, "c3", "s", "l", "o2", 2, 1],
        ["2", 4, 5, 3, -3, 3, 4, "c3", "s", "l", "o2", 2, 1],
        ["2", 6, 7, 3, -3, 3, 4, "c3", "s", "l", "o2", 2, 1],

        ["1", 0, 1, 4, -4, 4, 4, "c4", "s", "l", "o2", 2, 1],
        ["1", 2, 3, 4, -4, 4, 4, "c4", "s", "l", "o2", 2, 1],
        ["1", 4, 5, 4, -4, 4, 4, "c4", "s", "l", "o2", 2, 1],
        ["1", 6, 7, 4, -4, 4, 4, "c4", "s", "l", "o2", 2, 1],
        ["2", 0, 1, 4, -4, 4, 5, "c4", "s", "l", "o2", 2, 1],
        ["2", 2, 3, 4, -4, 4, 5, "c4", "s", "l", "o2", 2, 1],
        ["2", 4, 5, 4, -4, 4, 5, "c4", "s", "l", "o2", 2, 1],
        ["2", 6, 7, 4, -4, 4, 5, "c4", "s", "l", "o2", 2, 1],

        ["1", 0, 1, 5, -5, 5, 5, "c5", "s", "l", "o4", 4, 2],
        ["1", 2, 3, 5, -5, 5, 5, "c5", "s", "l", "o4", 4, 2],
        ["1", 4, 5, 5, -5, 5, 5, "c5", "s", "l", "o4", 4, 2],
        ["1", 6, 7, 5, -5, 5, 5, "c5", "s", "l", "o4", 4, 2],
        ["2", 0, 1, 5, -5, 5, 6, "c5", "s", "l", "o4", 4, 2],
        ["2", 2, 3, 5, -5, 5, 6, "c5", "s", "l", "o4", 4, 2],
        ["2", 4, 5, 5, -5, 5, 6, "c5", "s", "l", "o4", 4, 2],
        ["2", 6, 7, 5, -5, 5, 6, "c5", "s", "l", "o4", 4, 2],

        ["1", 0, 1, 6, -6, 6, 6, "c6", "s", "l", "o4", 4, 2],
        ["1", 2, 3, 6, -6, 6, 6, "c6", "s", "l", "o4", 4, 2],
        ["1", 4, 5, 6, -6, 6, 6, "c6", "s", "l", "o4", 4, 2],
        ["1", 6, 7, 6, -6, 6, 6, "c6", "s", "l", "o4", 4, 2],
        ["2", 0, 1, 6, -6, 6, 7, "c6", "s", "l", "o4", 4, 2],
        ["2", 2, 3, 6, -6, 6, 7, "c6", "s", "l", "o4", 4, 2],
        ["2", 4, 5, 6, -6, 6, 7, "c6", "s", "l", "o4", 4, 2],
        ["2", 6, 7, 6, -6, 6, 7, "c6", "s", "l", "o4", 4, 2],

    ],
    columns=['chr', 'start', 'end', 'reads', 'gc', 'copy', 'state', 'cell_id',
             'sample_id', 'library_id', 'origin_id', 'origin_id_int',
             'bhc_cluster_id']
)