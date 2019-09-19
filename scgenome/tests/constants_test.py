import numpy as np

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
