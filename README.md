A repository of code for analyzing single cell genomes

Filter cells using the cell copynumber classifier output:
```
filter_copynumber_cells copynumber_matrix.tsv all_metrics_summary_classified.csv copynumber_matrix_filt-cell.tsv
```

Filter copynumber genome bins using the QDNAseq blacklist:
```
filter_copynumber_bins copynumber_matrix_filt-cell.tsv QDNAseq_1.14.0_500kb_blacklist.tsv \
    copynumber_matrix_filt-cell_filt-bin.tsv
```

Filter contiguously duplicate copynumber genome bins:
```
filter_copynumber_contiguous_duplicate_bins copynumber_matrix_filt-cell_filt-bin.tsv \
    copynumber_matrix_filt-cell_filt-bin_filt-dup.tsv
```

Cluster cells:
```
cluster_cells copynumber_matrix_filt-cell_filt-bin_filt-dup.tsv \
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
