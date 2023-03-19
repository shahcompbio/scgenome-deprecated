import numpy as np
import pandas as pd
import pkg_resources


def read_chromosome_lengths(genome_fasta_index):
    chromosome_lengths = pd.read_csv(
        genome_fasta_index, sep='\t',
        header=None, names=['chr', 'chromosome_length', 'V3', 'V4', 'V5'])
    chromosome_lengths = chromosome_lengths[['chr', 'chromosome_length']]
    return chromosome_lengths


def read_cytobands(cyto_filename, remove_chr_prefix=False):
    cytobands = pd.read_csv(
        cyto_filename, sep='\t',
        names=['chr', 'start', 'end', 'cyto_band_name', 'cyto_band_giemsa_stain'])
    if remove_chr_prefix:
        cytobands['chr'] = cytobands['chr'].str.replace('^chr', '', regex=True)
    return cytobands


class RefGenomeInfo(object):
    def __init__(self, version):
        if version == 'hg19':
            self.chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
            self.plot_chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/hg19.fa.fai')
            self.cyto_filename = pkg_resources.resource_filename('scgenome', 'data/hg19_cytoBand.txt.gz')
            self.cytobands = read_cytobands(self.cyto_filename, remove_chr_prefix=True)

        elif version == 'grch38':
            self.chromosomes = [f'chr{a}' for a in range(1, 23)] + ['chrX', 'chrY']
            self.plot_chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/grch38.fa.fai')
            self.cyto_filename = pkg_resources.resource_filename('scgenome', 'data/grch38_cytoBand.txt.gz')
            self.cytobands = read_cytobands(self.cyto_filename, remove_chr_prefix=False)

        elif version == 'mm10':
            self.chromosomes = [str(a) for a in range(1, 20)] + ['X', 'Y']
            self.plot_chromosomes = [str(a) for a in range(1, 20)] + ['X', 'Y']
            self.genome_fasta_index = pkg_resources.resource_filename('scgenome', 'data/mm10.fa.fai')
            self.cyto_filename = pkg_resources.resource_filename('scgenome', 'data/mm10_cytoBand.txt.gz')
            self.cytobands = read_cytobands(self.cyto_filename, remove_chr_prefix=True)

        else:
            raise ValueError()

        self.chromosome_info = read_chromosome_lengths(self.genome_fasta_index)

        # Subset and order according to list of chromosomes
        self.chromosome_info = self.chromosome_info.set_index('chr').loc[self.chromosomes].reset_index()
        self.chromosome_info['chr_index'] = range(self.chromosome_info.shape[0])

        # Add plotting names of chromosomes
        self.chromosome_info['chr_plot'] = self.plot_chromosomes

        # Add start end and mid
        self.chromosome_info['chromosome_end'] = np.cumsum(self.chromosome_info['chromosome_length'])
        self.chromosome_info['chromosome_start'] = self.chromosome_info['chromosome_end'].shift(1)
        self.chromosome_info.loc[self.chromosome_info.index[0], 'chromosome_start'] = 0
        self.chromosome_info['chromosome_start'] = self.chromosome_info['chromosome_start'].astype(int)
        self.chromosome_info['chromosome_mid'] = (self.chromosome_info['chromosome_start'] + self.chromosome_info['chromosome_end']) // 2


info = None


def set_genome_version(version):
    global info
    info = RefGenomeInfo(version)

set_genome_version('hg19')

