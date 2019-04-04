import os
import wget
import pandas as pd

from datamanagement.miscellaneous.hdf5helper import read_python2_hdf5_dataframe


def get_snv_results(dest):
    print('starting load')

    mappability = read_python2_hdf5_dataframe(dest, '/snv/mappability')
    mappability = mappability[mappability['mappability'] > 0.99]
    mappability['chrom'] = mappability['chrom'].astype(str)

    strelka_results = read_python2_hdf5_dataframe(dest, '/strelka/vcf').rename(columns={'score': 'strelka_score'})
    print('strelka', strelka_results.shape)
    strelka_results = strelka_results[strelka_results['strelka_score'] > 20.]
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = read_python2_hdf5_dataframe(dest, '/museq/vcf').rename(columns={'score': 'museq_score'})
    print('museq', museq_results.shape)
    museq_results = museq_results[museq_results['museq_score'] > .9]
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = read_python2_hdf5_dataframe(dest,'/snv/cosmic_status')
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = read_python2_hdf5_dataframe(dest,'/snv/snpeff')[['chrom', 'coord', 'ref', 'alt', 'effect_impact']].drop_duplicates()
    snpeff['value'] = 1
    snpeff = snpeff.set_index(['chrom', 'coord', 'ref', 'alt', 'effect_impact'])['value'].unstack(fill_value=0)
    snpeff = snpeff.rename(columns=str).reset_index()

    tnc = read_python2_hdf5_dataframe(dest,'/snv/tri_nucleotide_context')

    data = read_python2_hdf5_dataframe(dest, '/snv_allele_counts')
    print('total', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(mappability)
    print('post mappability', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(strelka_results)
    print('post strelka', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(museq_results)
    print('post museq', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(cosmic, how='left')
    print('post cosmic', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(snpeff, how='left')
    print('post snpeff', data[['chrom', 'coord']].drop_duplicates().shape)
    data = data.merge(tnc, how='left')

    print('finishing load', data[['chrom', 'coord']].drop_duplicates().shape)

    return data


