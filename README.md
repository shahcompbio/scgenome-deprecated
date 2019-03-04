# A repository of code for analyzing single cell genomes

## Installation

The sisyphus repo is a prerequisite for the scgenome/analyses scripts.  Install as follows:

```
git clone https://github.com/shahcompbio/sisyphus
cd sisyphus
python setup.py develop
pip install -r requirements.txt
```

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
