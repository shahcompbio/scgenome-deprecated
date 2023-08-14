
def test_imports():
    import scgenome
    import scgenome.pl
    import scgenome.pp
    import scgenome.tl
    import scgenome.cnplot
    import scgenome.cncluster
    import scgenome.refgenome
    import scgenome.plotting.heatmap
    import scgenome.preprocessing.transform
    import scgenome.preprocessing.load_cn
    import scgenome.utils
    import scgenome.tools.ranges
    import scgenome.tools.getters

    print(scgenome.pl.plot_cn_profile)
    print(scgenome.pl.plot_var_profile)
    print(scgenome.pl.plot_cell_cn_matrix)
    print(scgenome.pl.plot_cell_cn_matrix_fig)
    print(scgenome.pl.cn_legend)
    print(scgenome.pl.plot_gc_reads)
    print(scgenome.pl.plot_tree_cn)


    print(scgenome.pp.filter_cells)
    print(scgenome.pp.calculate_filter_metrics)
    print(scgenome.pp.create_cn_anndata)
    print(scgenome.pp.read_dlp_hmmcopy)
    print(scgenome.pp.convert_dlp_hmmcopy)
    print(scgenome.pp.convert_dlp_signals)
    print(scgenome.pp.read_bam_bin_counts)
    print(scgenome.pp.read_medicc2_cn)
    print(scgenome.pp.read_snv_genotyping)


    print(scgenome.tl.cluster_cells)
    print(scgenome.tl.detect_outliers)
    print(scgenome.tl.aggregate_clusters_hmmcopy)
    print(scgenome.tl.aggregate_clusters)
    print(scgenome.tl.compute_umap)
    print(scgenome.tl.pca_loadings)
    print(scgenome.tl.sort_cells)
    print(scgenome.tl.sort_clusters)
    print(scgenome.tl.count_gc)
    print(scgenome.tl.mean_from_bigwig)
    print(scgenome.tl.add_cyto_giemsa_stain)
    print(scgenome.tl.read_ensemble_genes_gtf)
    print(scgenome.tl.aggregate_genes)
    print(scgenome.tl.ad_concat_cells)
    print(scgenome.tl.prune_leaves)
    print(scgenome.tl.align_cn_tree)
    print(scgenome.tl.create_bins)
    print(scgenome.tl.rebin)
    print(scgenome.tl.rebin_regular)
    print(scgenome.tl.weighted_mean)
    print(scgenome.tl.bin_width_weighted_mean)
    print(scgenome.tl.get_obs_data)

