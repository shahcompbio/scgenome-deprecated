import umap
import hdbscan
import seaborn
import pandas as pd
import matplotlib.pyplot as plt


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

    df = pd.Series(clusters, index=cn.columns, name='cluster_id').reset_index()
    df['umap1'] = embedding[:, 0]
    df['umap2'] = embedding[:, 1]
    df = df.dropna()

    return df


def plot_umap_clusters(ax, df):
    """ Scatter plot of umap clusters.

    Args:
        ax: matplotlib axis
        df: clusters dataframe with columns:
            cluster_id
            umap1
            umap2

    """

    df_noise = df[df['cluster_id'] < 0]
    ax.scatter(
        df_noise['umap1'].values,
        df_noise['umap2'].values,
        c='0.85',
        s=2,
        label='Filt.',
    )

    num_colors = len(df[df['cluster_id'] >= 0]['cluster_id'].unique())
    idx = 0.
    for cluster_id, cluster_df in df[df['cluster_id'] >= 0].groupby('cluster_id'):
        c = plt.get_cmap("tab10")((idx) / (num_colors - 1))
        ax.scatter(
            cluster_df['umap1'].values,
            cluster_df['umap2'].values,
            c=c,
            s=2,
            label=int(cluster_id),
        )
        idx += 1.

    ax.legend(
        frameon=False, markerscale=5,
        scatterpoints=1, bbox_to_anchor=(0.96, 0.85))
    ax.set_xlabel('Comp. 1')
    ax.set_ylabel('Comp. 2')
    seaborn.despine(ax=ax, offset=0, trim=True)


