import umap
import hdbscan
import seaborn
import logging
import itertools
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sklearn.cluster
import scipy.spatial
from adjustText import adjust_text


def umap_hdbscan_cluster(
        cn,
        n_components=2,
        n_neighbors=15,
        min_dist=0.1,
    ):
    """ Cluster using umap and hdbscan.

    Args:
        cn: data frame columns as cell ids, rows as segments

    Returns:
        data frame with columns:
            cluster_id
            cell_id
            umap1
            umap2

    """
    embedding = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=42,
        metric='euclidean',
    ).fit_transform(cn.fillna(0).values.T)

    clusters = hdbscan.HDBSCAN(
        min_samples=10,
        min_cluster_size=30,
    ).fit_predict(embedding)

    df = pd.DataFrame({
        'cell_id': cn.columns, 'cluster_id': clusters,
        'umap1': embedding[:, 0], 'umap2': embedding[:, 1]
    })
    df = df[['cell_id', 'cluster_id', 'umap1', 'umap2']]
    df = df.dropna()

    return df


def compute_bic(kmeans, X):
    """ Computes the BIC metric for a given k means clustering

    Args:
        kmeans: a fitted kmeans clustering object
        X: data for which to calculate bic
    
    Returns:
        float: bic
    
    Reference: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
    """
    centers = [kmeans.cluster_centers_]
    labels  = kmeans.labels_
    n_clusters = kmeans.n_clusters
    cluster_sizes = np.bincount(labels)
    N, d = X.shape

    # Compute variance for all clusters
    cl_var = (1.0 / (N - n_clusters) / d) * sum([sum(scipy.spatial.distance.cdist(X[np.where(labels == i)], [centers[0][i]], 
             'euclidean')**2) for i in range(n_clusters)])

    const_term = 0.5 * n_clusters * np.log(N) * (d+1)

    bic = np.sum([cluster_sizes[i] * np.log(cluster_sizes[i]) -
               cluster_sizes[i] * np.log(N) -
             ((cluster_sizes[i] * d) / 2) * np.log(2*np.pi*cl_var) -
             ((cluster_sizes[i] - 1) * d/ 2) for i in range(n_clusters)]) - const_term

    return bic


def kmeans_cluster(
        cn,
        min_k=2,
        max_k=100,
    ):
    """ Cluster using kmeans and bic.
    """

    X = cn.T.values
    ks = range(min_k, max_k + 1)

    logging.info(f'trying with max k={max_k}')

    kmeans = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.cluster.KMeans(n_clusters=k, init="k-means++").fit(X)
        bic = compute_bic(model, X)
        kmeans.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = kmeans[opt_k]

    embedding = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        metric='euclidean',
    ).fit_transform(cn.fillna(0).values.T)

    clusters = pd.DataFrame({
        'cell_id': cn.columns, 'cluster_id': model.labels_,
        'umap1': embedding[:, 0], 'umap2': embedding[:, 1]
    })

    return clusters


def get_cluster_palette(n_col):
    if n_col <= 10:
        palette = plt.get_cmap("tab10")
    elif n_col <= 21:
        palette = mpl.colors.ListedColormap([
            '#1d1d1d', '#ebce2b', '#702c8c', '#db6917', '#96cde6', '#ba1c30',
            '#c0bd7f', '#7f7e80', '#5fa641', '#d485b2', '#4277b6', '#df8461',
            '#463397', '#e1a11a', '#91218c', '#e8e948', '#7e1510', '#92ae31',
            '#6f340d', '#d32b1e', '#2b3514'
        ])
    else:
        palette = plt.get_cmap("hsv")
    return palette


def get_cluster_color_map(cluster_ids):
    num_colors = len(np.unique(cluster_ids[cluster_ids >= 0]))
    pal = get_cluster_palette(num_colors)

    color_map = {}

    cluster_ids = np.sort(np.unique(cluster_ids))

    for cluster_id in np.sort(np.unique(cluster_ids)):
        if cluster_id < 0:
            color_map[cluster_id] = (0.75, 0.75, 0.75, 1.0)

    cluster_ids = cluster_ids[cluster_ids >= 0]

    idx = 0.
    for cluster_id in itertools.chain(cluster_ids[::2], cluster_ids[1::2]):
        color_map[cluster_id] = pal(float(idx) / float(num_colors - 1))
        idx += 1

    return color_map


def get_cluster_colors(cluster_ids, color_map=None, return_map=False):
    if color_map is None:
        color_map = get_cluster_color_map(cluster_ids)

    color_mat = []
    for cluster_id in cluster_ids:
        color_mat.append(color_map[cluster_id])

    if return_map:
        return color_mat, color_map 

    return color_mat


def cluster_labels(cluster_ids):
    counts = cluster_ids.value_counts().astype(str)
    labels = counts.index.to_series().astype(str) + ' (' + counts + ')'
    if -1 in labels:
        labels.at[-1] = labels[-1].replace('-1', 'Filt.')
    return dict(list(zip(counts.index, labels)))


def plot_umap_clusters(ax, df):
    """ Scatter plot of umap clusters.

    Args:
        ax: matplotlib axis
        df: clusters dataframe with columns:
            cluster_id
            umap1
            umap2

    """
    labels = cluster_labels(df['cluster_id'])
    color_map = get_cluster_color_map(df['cluster_id'].values)

    if -1 in labels:
        df_noise = df[df['cluster_id'] < 0]
        ax.scatter(
            df_noise['umap1'].values,
            df_noise['umap2'].values,
            color=color_map[-1],
            s=2,
            label=labels[-1],
        )

    text_labels = []
    for cluster_id, cluster_df in df[df['cluster_id'] >= 0].groupby('cluster_id'):
        ax.scatter(
            cluster_df['umap1'].values,
            cluster_df['umap2'].values,
            color=color_map[cluster_id],
            s=2,
            label=labels[int(cluster_id)],
        )

    label_pos = df.groupby('cluster_id').mean()
    text_labels = [
        ax.text(label_pos.at[c, 'umap1'], label_pos.at[c, 'umap2'], c)
        for c in list(labels.keys()) if c >= 0
    ]
    adjust_text(
        text_labels, ax=ax,
        force_points=(0.1, 0.1)
    )

    ax.legend(
        frameon=False, markerscale=5,
        scatterpoints=1, bbox_to_anchor=(0.96, 0.85))
    ax.set_xlabel('Comp. 1')
    ax.set_ylabel('Comp. 2')
    seaborn.despine(ax=ax, offset=0, trim=True)


