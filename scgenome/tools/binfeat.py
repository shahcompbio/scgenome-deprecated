import pyfaidx
import pyBigWig
import numba
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


@numba.jit(nopython=True)
def mean_ranges(values, sequence, starts, ends):
    for i in range(starts.shape[0]):
        values[i] = np.nanmean(sequence[starts[i]:ends[i]])


def _chromosome_mean_bigwig(data, **kwargs):
    bigwig_file = kwargs['bigwig_file']
    chromosome = kwargs['chromosome']
    column_name = kwargs['column_name']
    chr_prefix = kwargs['chr_prefix']

    chromosome = chr_prefix + chromosome
    
    bw = pyBigWig.open(bigwig_file, 'r')

    if chromosome not in bw.chroms():
        data[column_name] = None
        return data

    chromosome_length = bw.chroms()[chromosome]
    sequence = np.array(bw.values(chromosome, start=0, end=chromosome_length))

    x = np.zeros(data.shape[0])
    mean_ranges(x, sequence, data['Start'].values, data['End'].values)

    data[column_name] = x

    return data


def mean_from_bigwig(bins, bigwig_file, column_name, chr_prefix=''):
    """ Count gc in each bin

    Parameters
    ----------
    bins : pyranges.PyRanges
        ranges for which to count gc
    bigwig_file : str
        bigwig filename
    column_name : str
        column to add to `bins`
    chr_prefix : str
        prefix for chromosome names, default ''

    Returns
    -------
    pyranges.PyRanges
        output ranges with additional column for mean bigwig value per bin
    """    

    data = bins.apply(
        _chromosome_mean_bigwig,
        bigwig_file=bigwig_file,
        column_name=column_name,
        chr_prefix=chr_prefix,
    )

    return data

