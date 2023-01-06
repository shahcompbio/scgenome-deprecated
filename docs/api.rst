.. module:: scgenome
.. automodule:: scgenome
   :noindex:

API
===================================

Import scgenome as::

   import scgenome

Preprocessing: `pp`
-------------------

.. module:: scgenome.pp
.. currentmodule:: scgenome

Data loading and pre-processing functionality.

Data loading
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   pp.read_dlp_hmmcopy
   pp.convert_dlp_hmmcopy
   pp.convert_dlp_signals
   pp.read_bam_bin_counts
   pp.read_snv_genotyping

Filtering
~~~~~~~~~

.. autosummary::
   :toctree: generated/

   pp.filter_cells
   pp.calculate_filter_metrics


Tools: `tl`
-----------

.. module:: scgenome.tl
.. currentmodule:: scgenome

Any transformation of the data matrix that is not *preprocessing*. In contrast to a *preprocessing* function, a *tool* usually adds an easily interpretable annotation to the data matrix, which can then be visualized with a corresponding plotting function.

Clustering
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.cluster_cells
   tl.cluster_cells_kmeans
   tl.aggregate_clusters_hmmcopy
   tl.aggregate_clusters
   tl.sort_cells

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.compute_umap
   tl.pca_loadings

Generating binned data
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.create_bins
   tl.count_gc
   tl.mean_from_bigwig

Gene regions
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.read_ensemble_genes_gtf
   tl.aggregate_genes

Phylogenetics
~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.prune_leaves
   tl.align_cn_tree

Anndata Manipulation
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.ad_concat_cells


Plotting: `pl`
--------------

.. module:: scgenome.pl
.. currentmodule:: scgenome

The plotting module :mod:`scgenome.pl` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. note::
   TODO: more plotting functions matching tools

Copy number profiles and heatmaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   pl.plot_cn_profile
   pl.plot_cell_cn_matrix
   pl.plot_cell_cn_matrix_fig
   pl.plot_gc_reads

Phylogenetics
~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   pl.plot_tree_cn

