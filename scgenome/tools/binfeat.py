import pyfaidx
import pyranges as pr
import numpy as np

import scgenome.refgenome


def create_bins(binsize):
    """ Create a regular binning of the genome

    Parameters
    ----------
    binsize : int
        length of bins

    Returns
    -------
    pyrange.PyRanges
        regular bins tiled across the genome
    """    
    
    chromsizes = scgenome.refgenome.info.chromosome_info[['chr', 'chromosome_length']].rename(
        columns={'chr': 'Chromosome', 'chromosome_length': 'End'}).assign(Start=0)[['Chromosome', 'Start', 'End']]

    chromsizes = pr.PyRanges(chromsizes)

    bins = pr.gf.tile_genome(chromsizes, binsize)

    return bins


def _chromosome_count_gc(data, **kwargs):
    genome_fasta = kwargs['genome_fasta']
    chromosome = kwargs['chromosome']
    column_name = kwargs['column_name']

    genome = pyfaidx.Fasta(genome_fasta)
    chrom_sequence = np.fromiter(str(genome[chromosome]).upper(), dtype='<U1')

    gc_indicator = ((chrom_sequence == 'G') | (chrom_sequence == 'C')) * 1
    gc_cumsum = gc_indicator.cumsum()

    data[column_name] = gc_cumsum[data['End'].values - 1] - gc_cumsum[data['Start'].values]

    return data


def count_gc(bins, genome_fasta, column_name='gc', proportion=False):
    """ Count gc in each bin

    Parameters
    ----------
    bins : pyranges.PyRanges
        ranges for which to count gc
    genome_fasta : str
        reference genome fasta
    column_name : str, optional
        column to add to `bins`, by default 'gc'
    proportion : bool, optional
        proportion of length, by default False

    Returns
    -------
    pyranges.PyRanges
        output ranges with additional column for gc count
    """    

    data = bins.apply(_chromosome_count_gc, genome_fasta=genome_fasta, column_name=column_name)

    if proportion:
        data = data.assign(column_name, lambda df: df[column_name] / (df['End'] - df['Start']))

    return data

