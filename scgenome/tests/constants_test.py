import numpy as np

E1_V1 = np.array([1, 2, 3, 4])[np.newaxis, :]
E1_V2 = E1_V1 + 1
E1_V3 = np.array([2, 4, 6, 8])[np.newaxis, :]
E1_V4 = E1_V3 + 2
E1_CN_MAT = np.concatenate([E1_V1, E1_V2, E1_V3, E1_V4])

E1_CELL_IDS = ["c1", "c2", "c3", "c4"]
E1_CHR_NAMES = ["1", "2"]
