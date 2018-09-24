import umap
import hdbscan
import seaborn
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text


def umap_hdbscan_cluster(cn):
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
        n_neighbors=30,
        min_dist=0.0,
        n_components=2,
        random_state=42,
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


def cluster_labels(cluster_ids):
    counts = cluster_ids.value_counts().astype(str)
    labels = counts.index.to_series().astype(str) + ' (' + counts + ')'
    labels.at[-1] = labels[-1].replace('-1', 'Filt.')
    return dict(zip(counts.index, labels))


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

    df_noise = df[df['cluster_id'] < 0]
    ax.scatter(
        df_noise['umap1'].values,
        df_noise['umap2'].values,
        c='0.85',
        s=2,
        label=labels[-1],
    )

    num_colors = len(df[df['cluster_id'] >= 0]['cluster_id'].unique())
    pal = get_cluster_palette(num_colors)
    text_labels = []
    idx = 0.
    for cluster_id, cluster_df in df[df['cluster_id'] >= 0].groupby('cluster_id'):
        c = pal((idx) / (num_colors - 1))
        ax.scatter(
            cluster_df['umap1'].values,
            cluster_df['umap2'].values,
            c=c,
            s=2,
            label=labels[int(cluster_id)],
        )
        idx += 1.

    label_pos = df.groupby('cluster_id').mean()
    text_labels = [
        ax.text(label_pos.at[c, 'umap1'], label_pos.at[c, 'umap2'], c)
        for c in labels.keys() if c >= 0
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
