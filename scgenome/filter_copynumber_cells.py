from argparse import ArgumentParser
import pandas as pd


def get_args():
    descr = 'Remove bad cells from copynumber matrix'
    p = ArgumentParser(descr)

    p.add_argument('cn', help='copynumber tsv file')
    p.add_argument('scores', help='cell classification scores csv file')
    p.add_argument('filtered', help='output filtered tsv file')

    return p.parse_args()


def filter_cells(cn, scores):
    scores = scores[
        (scores['joint'] >= 0.75) &
        (scores['cell_call'].isin(['C1', 'C2']))
    ]
    for rm_cond in ['gDNA', 'GM', 'NCC', 'NTC']:
        mask = ~scores['experimental_condition'].str.contains(rm_cond)
        scores = scores[mask]

    good_cells = scores.loc[
        scores['cell_id'].isin(cn.columns.tolist()), 'cell_id'
    ].tolist()
    good_cells.sort()

    return cn[good_cells]


def main():
    argv = get_args()

    cn = pd.read_table(
        argv.cn, index_col=['chr', 'start', 'end', 'width'], dtype={'chr': str}
    )
    scores = pd.read_csv(argv.scores)

    filtered_cn = filter_cells(cn, scores)
    filtered_cn.to_csv(argv.filtered, sep='\t')
    return


if __name__ == '__main__':
    main()
