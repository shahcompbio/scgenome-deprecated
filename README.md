A repository of code for analyzing single cell genomes

Filter cells using the cell copynumber classifier output:
```
filter_copynumber_cells copynumber_matrix.tsv all_metrics_summary_classified.csv copynumber_matrix_filt-cell.tsv
```

Filter copynumber genome bins using the QDNAseq blacklist:
```
filter_copynumber_bins copynumber_matrix_filt-cell.tsv QDNAseq_1.14.0_500kb_blacklist.tsv copynumber_matrix_filt-cell_filt-bin.tsv
```

Filter contiguously duplicate copynumber genome bins:
```
filter_copynumber_contiguous_duplicate_bins copynumber_matrix_filt-cell_filt-bin.tsv copynumber_matrix_filt-cell_filt-bin_filt-dup.tsv
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
make_cell_embedding_scatterplot copynumber_cell_embeddings.tsv copynumber_cell_clusters.tsv copynumber_embedding.pdf
```
