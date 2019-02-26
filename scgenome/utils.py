import pandas as pd

from . import refgenome


chrom_names = refgenome.info.chromosomes

chrom_idxs = pd.Series(chrom_names)
chrom_idxs.name = 'chr'
chrom_idxs.index.name = 'chr_index'
chrom_idxs = chrom_idxs.reset_index()


