from argparse import ArgumentParser
import pandas as pd


def get_args():
    descr = 'Remove bad bins from copynumber matrix'
    p = ArgumentParser(descr)

    p.add_argument('cn', help='copynumber tsv file')
    p.add_argument('blacklist', help='QDNAseq 500kb blacklist tsv file')
    p.add_argument('filtered', help='output filtered tsv file')

    return p.parse_args()


def filter_blacklist_bins(cn, blacklist):
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


def main():
    argv = get_args()

    cn = pd.read_table(argv.cn, dtype={'chr': str})
    blacklist = pd.read_table(argv.blacklist, true_values=['TRUE'])

    filtered_cn = filter_blacklist_bins(cn, blacklist)
    filtered_cn.to_csv(argv.filtered, sep='\t', index=False)
    return


if __name__ == '__main__':
    main()
