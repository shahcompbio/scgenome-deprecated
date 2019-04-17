import os
import wget
import pandas as pd

from datamanagement.miscellaneous.hdf5helper import read_python2_hdf5_dataframe


def get_highest_snpeff_effect(snpeff_data):
    """ Select the highest ranked effect from snpeff data.
    """
    ordered_effect_impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    snpeff_data = snpeff_data.merge(
        pd.DataFrame({
            'effect_impact': ordered_effect_impacts,
            'effect_impact_rank': range(len(ordered_effect_impacts))}))
    snpeff_data = (
        snpeff_data.sort_values(['chrom', 'coord', 'ref', 'alt', 'effect_impact_rank'], ascending=True)
        .groupby(['chrom', 'coord', 'ref', 'alt'], sort=False)
        .first()
        .reset_index())
    return snpeff_data


def get_snv_results(dest, museq_filter=0.9, strelka_filter=20.):
    print('starting load')

    mappability = read_python2_hdf5_dataframe(dest, '/snv/mappability')
    mappability = mappability[mappability['mappability'] > 0.99]
    mappability['chrom'] = mappability['chrom'].astype(str)

    strelka_results = read_python2_hdf5_dataframe(dest, '/strelka/vcf').rename(columns={'score': 'strelka_score'})
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = read_python2_hdf5_dataframe(dest, '/museq/vcf').rename(columns={'score': 'museq_score'})
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = read_python2_hdf5_dataframe(dest,'/snv/cosmic_status')
    print('cosmic', cosmic.shape)
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = read_python2_hdf5_dataframe(dest,'/snv/snpeff')
    snpeff = get_highest_snpeff_effect(snpeff)
    print('snpeff', snpeff.shape)

    tnc = read_python2_hdf5_dataframe(dest,'/snv/tri_nucleotide_context')

    data = read_python2_hdf5_dataframe(dest, '/snv_allele_counts')
    print('total', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(mappability)
    print('post mappability', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(strelka_results, how='left')
    print('post strelka', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(museq_results, how='left')
    print('post museq', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(cosmic, how='left')
    print('post cosmic', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(snpeff, how='left')
    print('post snpeff', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(tnc, how='left')

    if museq_filter is not None:
        data = data[data['museq_score'] > museq_filter]
        print('post museq filter', data[['chrom', 'coord']].drop_duplicates().shape)

    if strelka_filter is not None:
        data = data[data['strelka_score'] > strelka_filter]
        print('post strelka filter', data[['chrom', 'coord']].drop_duplicates().shape)

    print('finishing load', data[['chrom', 'coord']].drop_duplicates().shape)

    return data


