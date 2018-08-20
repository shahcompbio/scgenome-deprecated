from argparse import ArgumentParser
import pandas as pd


def get_args():
    descr = 'Remove bad and duplicated bins from copynumber matrix'
    p = ArgumentParser(descr)

    p.add_argument('cn', help='copynumber tsv file')
    p.add_argument('filtered', help='output filtered tsv file')

    return p.parse_args()


def remove_contiguous_duplicates(cn):
    cn.sort_values(by=['chr', 'start', 'end'], inplace=True)
    cn.set_index(['start', 'end', 'width'], inplace=True)
    cn = cn[(cn.shift() != cn).any(axis=1)]
    cn = cn.reset_index().set_index(['chr', 'start', 'end', 'width'])
    cn.reset_index(inplace=True)
    return cn


def main():
    argv = get_args()

    cn = pd.read_table(argv.cn, dtype={'chr': str})

    filtered_cn = remove_contiguous_duplicates(cn)
    filtered_cn.to_csv(argv.filtered, sep='\t', index=False)
    return


if __name__ == '__main__':
    main()
