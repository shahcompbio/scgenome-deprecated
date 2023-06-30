import numpy as np
import pandas as pd
import pkg_resources


class RefGenome(object):
    def __init__(self, genome_version='hg19'):
        assert genome_version in ['hg19', 'grch38', 'mm10']
        self.genome_version = genome_version

    def chromosomes(self):
        chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
        if self.genome_version == 'grch38':
            chromosomes = [f'chr{v}' for v in chromosomes]
        return chromosomes

    def plot_chromosomes(self):
        return [str(a) for a in range(1, 23)] + ['X', 'Y']

    def cytobands(self):
        filepath = pkg_resources.resource_filename(
            'scgenome', f'data/{self.genome_version}_cytoBand.txt.gz'

        )
        cytobands = pd.read_csv(
            filepath, sep='\t',
            names=['chr', 'start', 'end', 'cyto_band_name', 'cyto_band_giemsa_stain']
        )

        if not self.genome_version == 'grch38':
            cytobands['chr'] = cytobands['chr'].str.replace('^chr', '', regex=True)
        return cytobands

    def chromosome_info(self):
        filepath = pkg_resources.resource_filename('scgenome', f'data/{self.genome_version}.fa.fai')

        chromosome_info = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=['chr', 'chromosome_length', 'V3', 'V4', 'V5']
        )
        chromosome_info = chromosome_info[['chr', 'chromosome_length']]

        # Subset and order according to list of chromosomes
        chromosome_info = chromosome_info.set_index('chr').loc[self.chromosomes()].reset_index()
        chromosome_info['chr_index'] = range(chromosome_info.shape[0])

        # Add plotting names of chromosomes
        chromosome_info['chr_plot'] = self.plot_chromosomes()

        # Add start end and mid
        chromosome_info['chromosome_end'] = np.cumsum(chromosome_info['chromosome_length'])
        chromosome_info['chromosome_start'] = chromosome_info['chromosome_end'].shift(1)
        chromosome_info.loc[chromosome_info.index[0], 'chromosome_start'] = 0
        chromosome_info['chromosome_start'] = chromosome_info['chromosome_start'].astype(int)
        chromosome_info['chromosome_mid'] = (chromosome_info['chromosome_start'] + chromosome_info[
            'chromosome_end']) // 2

        return chromosome_info


def chromosome_info(genome_version='hg19'):
    return RefGenome(genome_version=genome_version).chromosome_info()


def cytobands(genome_version='hg19'):
    return RefGenome(genome_version=genome_version).cytobands()


def chromosomes(genome_version='hg19'):
    return RefGenome(genome_version=genome_version).chromosomes()

def plot_chromosomes(genome_version='hg19'):
    return RefGenome(genome_version=genome_version).plot_chromosomes()
