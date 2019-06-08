# A repository of code for analyzing single cell genomes

## Installation

It is recommended that you install all prerequisites with pip in a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Analyses

### Cell cycle

The `compute_cell_cycle_state.py` script can be used to run the cell cycle classifier and upload results into tantalus.

### Extract cellenone features

The `extract_cellenone_features.py` script searches for unprocessed cellenone data and extracts a table of features (Diameter, Elongation, Circularity) from the cellenone tables.

### Clonal inference

The `infer_clones.py` script can be used to run the full clonal analysis on a library or set of libraries.

#### Usage

The `infer_clones.py` script runs in 3 stages: `retrieve_cn`, `cluster_cn`, and `pseudobulk_analysis.py`, to be run in that order.

The CLI requires you to select the stage, and also specify the results prefix and the local storage directory for caching files.  Results of the analysis will be stored in files named with the results prefix and tantalus data will be cached in the local storage directory.

The `retrieve_cn` stage requests metadata from tantalus and caches the data locally.  This stage requires one or more library ids and sample ids.  Copy number tables will be stored with the results prefix.

The `cluster_cn` stage runs the current copy number clustering, produces tables including the cluster labels, and distance from each cell to each cluster.  Additional plots will be output to files with the results prefix.

The `pseudobulk_analysis.py` stage runs the current pseudobulk analyses on a given raw pseudobulk run.  This includes:
- an analysis of the data as bulk with plots that include mutation signatures and genome wide SNV frequencies
- SNV phylogenetic analysis of clone specific SNVs
- Inference of allele specific copy number
Tables and plots are output with the results prefix.

## Scripts

Filter cells and bins:
```
filter_copynumber copynumber_matrix.tsv copynumber_matrix_filt.tsv \
    --cell-scores all_metrics_summary_classified.csv \
    --qdnaseq-blacklist QDNAseq_1.14.0_500kb_blacklist.tsv
```

Same, but also remove consecutive bins with the same copy number across cells:
```
filter_copynumber copynumber_matrix.tsv copynumber_matrix_filt.tsv \
    --cell-scores all_metrics_summary_classified.csv \
    --qdnaseq-blacklist QDNAseq_1.14.0_500kb_blacklist.tsv \
    --filter-contig-dup-bins
```

Cluster cells:
```
cluster_cells copynumber_matrix_filt.tsv \
    copynumber_cell_clusters.tsv \
    --plot copynumber_cell_clusters.pdf
```

![cell cluster scatterplot](https://user-images.githubusercontent.com/381464/45980923-56f2b300-c021-11e8-9b0e-9dcf4b53f9c7.png)

There is a sample Snakemake file included in the pipelines directory. You can run it like this:
```
snakemake --config \
    copynumber_matrix=cn_matrix.csv \
    classified=all_metrics_summary_classified.csv \
    qdnaseq_blacklist=QDNAseq_1.14.0_500kb_blacklist.tsv
```
where classified is a file with cell classifications, and qndaseq_blacklist is a tsv file generated from QDNAseq with the following columns:
* chromosome
* start
* end
* bases
* gc
* mappability
* blacklist
* residual
* use
