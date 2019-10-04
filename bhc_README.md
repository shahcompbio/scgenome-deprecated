# Bayesian Hierarchical Clustering
The function that does BHC is `scgenome/cncluster::bayesian_cluster()`. The
majority of mathematical implementation is in `scgenome/TNode.py`.
## Useful files
Good example of how to use BHC code
### Scripts
In `scgenome/scripts/bhc`

`do_bhc.py`: Script for running BHC on a dataset

`do_clustering.py`: Does BHC, Naive and UMAP+HDBSCAN to cluster, plots results

`do_naive.py`: Script for running naive clustering on a dataset

`heatmap.py`: Script for plotting heatmap for the result of a clustering run.
(eg. result of `do_bhc.py`)

`do_qc.py`: Do quality control on an HMMCopy results dataset

`multi_clustering_vary_num_read.py`: Do BHC, Naive and UMAP/HDBSCAN on 
submixtures from datasets with only cells below a certain read count. 
Calculate accuracy, plot results

`multi_clustering_vary_prop.py`: Do BHC, Naive and UMAP/HDBSCAN on 
Do BHC, Naive and UMAP/HDBSCAN on mixtures of sample, varying the proportion of
samples in mixture and plot results

`select_prune.py`: Varies tree pruning threshold for a linkage matrix generated 
either by BHC or Naive clustering, plots the resulting number of clusters. 
Used for selecting pruning threshold.

`single_BHC.py`: Runs bayesian hierarchical clustering on a simulation of 
poisson random walk

### Notebooks
In `notebooks`

`test_cnplot.ipynb`: Toy example of plotting heatmap with dendrogram

`bayesian_spike_in-Copy1.ipynb`: Example of running BHC on mixture experiment

`bhc_runtime.ipynb`: Plots runtime of BHC

## Example workflow
Given an HMMCopy results dataset, we do the following:
`do_qc.py` --> QCed HMMCopy dataset -->

`do_clustering.py` --> BHC, Naive linkage matrices; UMAP+HDBSCAN clustering; 
heatmaps with clusterings and dendrograms for each -->

`select_prune.py` on bhc linkage --> plots of threshold v. num_clusters; 
heatmap with dendrogram using with better clustering since we picked a 
different threshold
## Code Quirks
### Cell id ordering
`bayesian_cluster()` and 
`plot_clustered_cell_cn_matrix_figure()` both turn the long parameter `cn_data`
into a muli-indexed pandas DataFrame. However, they use slightly different 
methods so the resulting arrays don't always have the same order. Thus,
`plot_clustered_cell_cn_matrix_figure()` takes an argument `cell_id_order` 
that, after making `cn_data` a multi-indexed DataFrame uses the `cell_id` 
column so the order of the multi-indexed DataFrame matches that specified
by `cell_id_order`. `bayesian_cluster()` returns a list of cell_ids that 
should be passed into `plot_clustered_cell_cn_matrix_figure()` whenever it 
is being used to plot a dendrogram.

### Mixture experiment
The phrase "mixture experiment" and "spike in" are used interchangably
throughout comments and code. Both refer to when we take a set of samples
eg. SC-1935, SC-1936 and SC-1937, mix them in known proportion 
(eg. 33%, 33%, 34%), and cluster the mixed dataset, and see the clustering
algorithm's ability to recall what cell came from what sample.
