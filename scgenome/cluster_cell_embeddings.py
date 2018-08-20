from argparse import ArgumentParser
import hdbscan
import numpy as np
import pandas as pd
import string

np.random.seed(2794834348)


def get_args():

    descr = 'Cluster cells after UMAP dimensionality reduction'
    p = ArgumentParser(descr)

    p.add_argument('embedding', help='cell embedding tsv file')
    p.add_argument('clusters', help='output cell clusters tsv file')
    p.add_argument(
        '--min-samples', type=int, default=10,
        help='HDBSCAN min_samples argument'
    )
    p.add_argument(
        '--min-cluster-size', type=int, default=20,
        help='HDBSCAN min_cluster_size argument'
    )

    return p.parse_args()


def cluster_cell_embeddings(embedding, min_samples, min_cluster_size):
    # low min_samples, quite large min_cluster_size
    clusterer = hdbscan.HDBSCAN(
        min_samples=min_samples, min_cluster_size=min_cluster_size,
        approx_min_span_tree=False
    )
    cluster_labels = clusterer.fit_predict(embedding)

    if cluster_labels.max() <= 26:
        cluster_labels = [
            string.ascii_letters[i] if i > -1 else i for i in cluster_labels
        ]
    else:
        cluster_labels = [str(x) for x in cluster_labels]

    cluster_labels = pd.DataFrame(
        cluster_labels, columns=['cluster'], index=embedding.index
    )
    return cluster_labels


if __name__ == '__main__':
    argv = get_args()

    embedding = pd.read_table(argv.embedding, index_col='cell')
    cluster_labels = cluster_cell_embeddings(
        embedding, argv.min_samples, argv.min_cluster_size
    )
    cluster_labels.to_csv(argv.clusters, sep='\t')
