import numpy as np
import pandas as pd
import pkg_resources


def read_chromosome_lengths(genome_fasta_index):
    fai = pd.read_csv(genome_fasta_index, sep='\t', header=None, names=['chrom', 'length', 'V3', 'V4', 'V5'])
    fai = fai.set_index('chrom')['length']
    return fai.to_dict()


class RefGenomeInfo(object):
    def __init__(self, version):
        if version == 'hg19':
            self.chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/GRCh37-lite.fa.fai')

        elif version == 'grch38':
            self.chromosomes = [f'chr{a}' for a in range(1, 23)] + ['chrX', 'chrY']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/GRCh38.fa.fai')

        elif version == 'mm10':
            self.chromosomes = [str(a) for a in range(1, 20)] + ['X', 'Y']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/mm10.fa.fai')

        else:
            raise ValueError()

        self.chromosome_length = pd.Series(read_chromosome_lengths(self.genome_fasta_index)).reindex(self.chromosomes).astype(int)
        self.chromosome_length.index.name = 'chr'

        self.chromosome_end = np.cumsum(self.chromosome_length)
        self.chromosome_start = self.chromosome_end.shift(1)
        self.chromosome_start[0] = 0
        self.chromosome_start = self.chromosome_start.astype(int)
        self.chromosome_mid = (self.chromosome_start + self.chromosome_end) / 2.

        self.chromosome_info = pd.DataFrame({
            'chromosome_length': self.chromosome_length,
            'chromosome_end': self.chromosome_end,
            'chromosome_start': self.chromosome_start,
            'chromosome_mid': self.chromosome_mid,
        }).reset_index()

        self.chrom_idxs = pd.Series(self.chromosomes)
        self.chrom_idxs.name = 'chr'
        self.chrom_idxs.index.name = 'chr_index'
        self.chrom_idxs = self.chrom_idxs.reset_index()


info = None


def set_genome_version(version):
    global info
    info = RefGenomeInfo(version)

set_genome_version('hg19')

