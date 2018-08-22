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

Generate cell embeddings:
```
reduce_copynumber_dims copynumber_matrix_filt-cell_filt-bin_filt-dup.tsv copynumber_cell_embeddings.tsv
```

Cluster cell embeddings:
```
cluster_cell_embeddings copynumber_cell_embeddings.tsv copynumber_cell_clusters.tsv
```

Plot cell embeddings:
```
make_cell_embedding_scatterplot copynumber_cell_embeddings.tsv copynumber_cell_clusters.tsv \
    copynumber_embedding.pdf
```

![cell embedding scatterplot](https://user-images.githubusercontent.com/381464/44435403-afd3b500-a564-11e8-9365-98c5b66ac202.png)

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
