import pandas as pd


def filter_cells(cn, scores, threshold=0.75):
    cn.set_index(['chr', 'start', 'end', 'width'], inplace=True)

    scores = scores[
        (scores['joint'] >= threshold) &
        (scores['cell_call'].isin(['C1', 'C2']))
    ]

    for rm_cond in ['gDNA', 'GM', 'NCC', 'NTC']:
        mask = ~scores['experimental_condition'].str.contains(rm_cond)
        scores = scores[mask]

    good_cells = scores.loc[
        scores['cell_id'].isin(cn.columns.tolist()), 'cell_id'
    ].tolist()
    good_cells.sort()

    return cn[good_cells].reset_index()


def filter_qdnaseq_bins(cn, blacklist):
    blacklist.rename(columns={'chromosome': 'chr'}, inplace=True)

    cn = cn.merge(blacklist, how='left')

    # residual cutoff is 4 * madDiff(residual)
    cn = cn[
        (
            (
                ~pd.isnull(cn['residual']) &
                (cn['residual'].abs() <= 0.0696924)
            ) |
            pd.isnull(cn['residual'])
        ) &
        (cn['blacklist'] == 0)
    ]

    rm_cols = ['bases', 'gc', 'mappability', 'blacklist', 'residual', 'use']
    for rm_col in rm_cols:
        del cn[rm_col]

    return cn


def remove_contiguous_duplicate_bins(cn):
    cn.sort_values(by=['chr', 'start', 'end'], inplace=True)
    cn.set_index(['start', 'end', 'width'], inplace=True)
    cn = cn[(cn.shift() != cn).any(axis=1)]
    cn = cn.reset_index().set_index(['chr', 'start', 'end', 'width'])
    cn.reset_index(inplace=True)
    return cn
