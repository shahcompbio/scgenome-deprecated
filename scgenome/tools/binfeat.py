import pyfaidx
import pyBigWig
import numba
import numpy as np
import pandas as pd

import scgenome.refgenome
import scgenome.tools.ranges


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


def add_cyto_giemsa_stain(bins):
    """ Add bin specific giesma stain values

    Parameters
    ----------
    bins : pandas.DataFrame
        dataframe with columns 'chr', 'start', 'end' for which to add giesma stain values

    Returns
    -------
    pandas.DataFrame
        output dataframe with columns 'chr', 'start', 'end' and additional column for giesma stain values
    """    

    bins_pr = scgenome.tools.ranges.dataframe_to_pyranges(bins.rename_axis('_index').reset_index())
    cyto_pr = scgenome.tools.ranges.dataframe_to_pyranges(scgenome.refgenome.info.cytobands)

    intersect_1 = bins_pr.intersect(cyto_pr)
    intersect_2 = cyto_pr.intersect(bins_pr)

    intersect = pd.merge(
        scgenome.tools.ranges.pyranges_to_dataframe(intersect_1),
        scgenome.tools.ranges.pyranges_to_dataframe(intersect_2))

    intersect['_width'] = intersect['end'] - intersect['start'] + 1

    selected = intersect.sort_values('_width').drop_duplicates(['_index'], keep='last').set_index('_index')

    cols = ['cyto_band_name', 'cyto_band_giemsa_stain']
    bins = bins.merge(selected[cols], left_index=True, right_index=True, how='left')

    return bins


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

