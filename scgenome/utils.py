import pandas as pd


chrom_names = [str(v) for v in range(1, 23)] + ['X']

chrom_idxs = pd.Series(chrom_names)
chrom_idxs.name = 'chr'
chrom_idxs.index.name = 'chr_index'
chrom_idxs = chrom_idxs.reset_index()


