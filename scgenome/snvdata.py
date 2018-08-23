import os
import wget
import pandas as pd


def get_snv_results(sas, url, dest):
    print 'starting load'

    if not os.path.exists(dest):
        filename = wget.download(url + sas)
        os.rename(filename, dest)

    mappability = pd.DataFrame()
    for chunk in pd.read_hdf(dest, '/snv/mappability', chunksize=int(1e5)):
        mappability = pd.concat([mappability, chunk[chunk['mappability'] > 0.99]], ignore_index=True)
    mappability['chrom'] = mappability['chrom'].astype(str)

    store = pd.HDFStore(dest, 'r')

    strelka_results = store['/strelka/vcf'].rename(columns={'score': 'strelka_score'})
    print 'strelka', strelka_results.shape
    strelka_results = strelka_results[strelka_results['strelka_score'] > 20.]
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = store['/museq/vcf'].rename(columns={'score': 'museq_score'})
    print 'museq', museq_results.shape
    museq_results = museq_results[museq_results['museq_score'] > .9]
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = store['/snv/cosmic_status']
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = store['/snv/snpeff'][['chrom', 'coord', 'ref', 'alt', 'effect_impact']].drop_duplicates()
    snpeff['value'] = 1
    snpeff = snpeff.set_index(['chrom', 'coord', 'ref', 'alt', 'effect_impact'])['value'].unstack(fill_value=0)
    snpeff = snpeff.rename(columns=str).reset_index()

    tnc = store['/snv/tri_nucleotide_context']

    data = store['/snv_allele_counts']
    data = data.merge(mappability)
    data = data.merge(strelka_results)
    data = data.merge(museq_results)
    data = data.merge(cosmic, how='left')
    data = data.merge(snpeff, how='left')
    data = data.merge(tnc, how='left')

    print 'finishing load', data[['chrom', 'coord']].drop_duplicates().shape

    return data


