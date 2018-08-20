from argparse import ArgumentParser
import numpy as np
import pandas as pd
import umap

np.random.seed(2794834348)


def get_args():

    descr = 'Perform dimensionality reduction on cell copy number data'
    p = ArgumentParser(descr)

    p.add_argument('cn', help='cell copy number tsv file')
    p.add_argument('embedding', help='output cell embedding tsv file')

    return p.parse_args()


def reduce_dimensions(cn):
    # large n_neighbors, low min_dist
    embedder = umap.UMAP(n_neighbors=100, min_dist=0.001, metric='manhattan')
    embedding = embedder.fit_transform(cn.values.T)
    embedding = pd.DataFrame(
        embedding, columns=['x', 'y'],
        index=pd.Index(cn.columns, name='cell')
    )
    return embedding


def main():
    argv = get_args()

    cn = pd.read_table(
        argv.cn, index_col=['chr', 'start', 'end', 'width'], dtype={'chr': str}
    )
    embedding = reduce_dimensions(cn)
    embedding.to_csv(argv.embedding, sep='\t')


if __name__ == '__main__':
    main()
