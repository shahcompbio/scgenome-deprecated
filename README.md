# A repository of code for analyzing single cell genomes

## Installation

It is recommended that you install all prerequisites with pip in a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install numpy cython
pip install -r requirements.txt
python setup.py develop
```

Note that you will have to install numpy and cython prior to other requirements.

## Analyses

### Cell cycle

The `compute_cell_cycle_state.py` script can be used to run the cell cycle classifier and upload results into tantalus.

### Extract cellenone features

The `extract_cellenone_features.py` script searches for unprocessed cellenone data and extracts a table of features (Diameter, Elongation, Circularity) from the cellenone tables.

### Clonal inference

The `infer_clones.py` script can be used to run the full clonal analysis on a library or set of libraries.

#### Usage

The `infer_clones.py` script runs in 3 stages: `retrieve-cn`, `cluster-cn`, and `pseudobulk-analysis`, to be run in that order.

The CLI requires you to select the stage, and also specify the results prefix and the local storage directory for caching files.  Results of the analysis will be stored in files named with the results prefix and tantalus data will be cached in the local storage directory.

The `rretrieve-cn` stage requests metadata from tantalus and caches the data locally.  This stage requires one or more library ids and sample ids.  Copy number tables will be stored with the results prefix.

The `cluster-cn` stage runs the current copy number clustering, produces tables including the cluster labels, and distance from each cell to each cluster.  Additional plots will be output to files with the results prefix.

The `pseudobulk-analysis` stage runs the current pseudobulk analyses on a given raw pseudobulk run.  This includes:
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

## API

The API allows access to both HMMCopy and Pseudobulk data stored in blob and managed by tantalus.

### Prerequisites

#### Software

You should set up an environment with the requirements from `requirements.txt` including sisyphus.

#### Accounts

You must have a tantalus account and access to azure blob storage.  Tantalus and azure blob credentials should be in your environment.

### HMMCopy Data

The following example code snippet will provide access to HMMCopy data for the OV cell line data:

```
import dbclients
import scgenome.utils

from scgenome.loaders.qc import load_qc_data
from scgenome.db.qc import cache_qc_results


tantalus_api = dbclients.tantalus.TantalusApi()

hmmcopy_tickets = [
    'SC-1935',
    'SC-1936',
    'SC-1937',
]

sample_ids = [
    'SA1090',
    'SA921',
    'SA922',
]

local_cache_directory = '/Your/Local/Cache'

cn_data = []
segs_data = []
metrics_data = []
align_metrics_data = []

for jira_ticket in hmmcopy_tickets:
    analysis = tantalus_api.get(
        'analysis',
        analysis_type__name='hmmcopy',
        jira_ticket=jira_ticket)

    cache_qc_results(jira_ticket, local_cache_directory)
    ticket_directory = os.path.join(local_cache_directory, ticket_id)
    hmmcopy_data = load_qc_data(ticket_directory)

    cn_data.append(hmmcopy_data['hmmcopy_reads'])
    segs_data.append(hmmcopy_data['hmmcopy_segs'])
    metrics_data.append(hmmcopy_data['hmmcopy_metrics'])
    align_metrics_data.append(hmmcopy_data['align_metrics'])

cn_data = scgenome.utils.concat_with_categories(cn_data)
segs_data = scgenome.utils.concat_with_categories(segs_data)
metrics_data = scgenome.utils.concat_with_categories(metrics_data)
align_metrics_data = scgenome.utils.concat_with_categories(align_metrics_data)
```

### Pseudobulk Data

The following example code snippet will provide access to pseudobulk SNV, allele and breakpoint data for the OV cell line data:

```
import scgenome.snvdata
import scgenome.loaders.snv
import scgenome.loaders.allele
import scgenome.loaders.breakpoint

ticket_id = 'SC-1939'

results_prefix = './results'

local_cache_directory = '/Your/Local/Cache'

museq_score_threshold = None
strelka_score_threshold = None
snvs_num_cells_threshold = 2
snvs_sum_alt_threshold = 2

ticket_directory = os.path.join(local_cache_directory, ticket_id)

snv_results = scgenome.loaders.snv.load_snv_data(
    ticket_directory,
    museq_filter=museq_score_threshold,
    strelka_filter=strelka_score_threshold,
)

allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
    ticket_directory,
)

breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data(
    ticket_directory,
)

```

For additional filtering and annotation see `scgenome.analyses.infer_clones.retrieve_pseudobulk_data`.

