import umap
import hdbscan
import seaborn
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scgenome.jointcnmodels import get_variances, get_tr_probs
from itertools import combinations
from .TNode import TNode
from .constants import ALPHA, MAX_CN, VALUE_IDS
from .utils import cn_data_to_mat_data_ids


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


def get_cluster_palette(n_col):
    if n_col <= 10:
        palette = plt.get_cmap("tab10")
    else:
        palette = mpl.colors.ListedColormap([
            '#1d1d1d', '#ebce2b', '#702c8c', '#db6917', '#96cde6', '#ba1c30',
            '#c0bd7f', '#7f7e80', '#5fa641', '#d485b2', '#4277b6', '#df8461',
            '#463397', '#e1a11a', '#91218c', '#e8e948', '#7e1510', '#92ae31',
            '#6f340d', '#d32b1e', '#2b3514'
        ])
    return palette


def get_cluster_color_map(cluster_ids):
    num_colors = len(np.unique(cluster_ids[cluster_ids >= 0]))
    pal = get_cluster_palette(num_colors)

    color_map = {}
    idx = 0.
    for cluster_id in np.sort(np.unique(cluster_ids)):
        if cluster_id < 0:
            color_map[cluster_id] = (0.75, 0.75, 0.75, 1.0)
        else:
            color_map[cluster_id] = pal((idx) / (num_colors - 1))
            idx += 1

    return color_map


def get_cluster_colors(cluster_ids):
    color_map = get_cluster_color_map(cluster_ids)

    color_mat = []
    for cluster_id in cluster_ids:
        color_mat.append(color_map[cluster_id])

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
            c=color_map[-1],
            s=2,
            label=labels[-1],
        )

    text_labels = []
    for cluster_id, cluster_df in df[df['cluster_id'] >= 0].groupby('cluster_id'):
        ax.scatter(
            cluster_df['umap1'].values,
            cluster_df['umap2'].values,
            c=color_map[cluster_id],
            s=2,
            label=labels[int(cluster_id)],
        )

    label_pos = df.groupby('cluster_id').mean()
    text_labels = [
        ax.text(label_pos.at[c, 'umap1'], label_pos.at[c, 'umap2'], c)
        for c in list(labels.keys()) if c >= 0
    ]
    adjust_text(
        text_labels, ax=ax, x=df['umap1'], y=df['umap2'],
        force_points=(0.1, 0.1)
    )

    ax.legend(
        frameon=False, markerscale=5,
        scatterpoints=1, bbox_to_anchor=(0.96, 0.85))
    ax.set_xlabel('Comp. 1')
    ax.set_ylabel('Comp. 2')
    seaborn.despine(ax=ax, offset=0, trim=True)


# TODO set alpha
# TODO save more info as needed eg. Rs of subtrees
# TODO maybe cache values
# TODO next_level, r are redundant
# TODO return more stuff
def bayesian_cluster(cn_data, cluster_col="bayes_cluster_id", n_states=MAX_CN,
                     alpha=ALPHA, value_ids=VALUE_IDS):
    matrix_data, measurement, cell_ids = (
        cn_data_to_mat_data_ids(cn_data, value_ids=value_ids))
    n_cells = measurement.shape[0]
    n_segments = measurement.shape[1]
    variances = get_variances(cn_data, matrix_data, n_states)
    tr_probs = get_tr_probs(n_segments, n_states)
    tr_mat = np.log(tr_probs)

    clusters = [TNode([i], None, None, i, alpha) for i in range(n_cells)]
    #def __init__(self, sample_inds, left_child, right_child, cluster_ind,
    #             alpha=ALPHA, pi=None, d=None, ll=None, log_r=None,
    #             marginal=None):
    linkage = pd.DataFrame(data=None,
                           columns=["i", "j", "r_merge", "i_count", "j_count"],
                           index=list(range(n_cells-1)))
    li = 0
    while len(clusters) > 1:
        r = np.empty((len(clusters), len(clusters)))
        r.fill(np.nan)
        next_level = [[None for i in range(len(clusters))]
                      for j in range(len(clusters))]
        for i, j in combinations(range(len(clusters)), 2):
            left_cluster = clusters[i]
            right_cluster = clusters[j]
            merge_cluster = TNode(
                clusters[i].sample_inds + clusters[j].sample_inds,
                left_cluster, right_cluster, None, alpha,
            )

            pi, d = merge_cluster.get_pi_d()
            ll = merge_cluster.get_ll(measurement, variances, tr_mat)
            r[i, j] = merge_cluster.get_log_r(measurement, variances, tr_mat)
            next_level[i][j] = merge_cluster

        max_r_flat_ind = np.nanargmax(r)
        i_max, j_max = np.unravel_index(max_r_flat_ind, r.shape)
        selected_cluster = next_level[i_max][j_max]

        selected_cluster.cluster_ind = n_cells + li
        left_ind = selected_cluster.left_child.cluster_ind
        right_ind = selected_cluster.right_child.cluster_ind
        linkage.iloc[li] = [left_ind, right_ind, r.flatten()[max_r_flat_ind],
                            len(clusters[i_max].sample_inds),
                            len(clusters[j_max].sample_inds)]

        li += 1
        clusters[i_max] = selected_cluster
        del clusters[j_max]

    linkage["merge_count"] = linkage["i_count"] + linkage["j_count"]
    return linkage, clusters[0], cell_ids
