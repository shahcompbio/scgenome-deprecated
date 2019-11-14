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

The list of credentials below are also required in your environment at runtime:
```
TANTALUS_API_PASSWORD
TANTALUS_API_USERNAME
CLIENT_ID
TENANT_ID
SECRET_KEY
AZURE_KEYVAULT_ACCOUNT
```

For access to Tantalus, please contact [Diljot Grewal](mailto:dgrewal@bccrc.ca).

For information on the other environment variables, please see [here](https://shahcompbio.atlassian.net/wiki/spaces/SOFT/pages/26378254/Scgenome+and+Tantalus+credentials+for+blob+access).

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

## Retrieving QC data

QC data is referenced by the jira ticket of the analysis that created that data.  The primary data store is in the `singlecellresults` azure blob storage.  A secondary store is on juno in the directory `/work/shah/tantalus/`.  Filesystem layout is `{storage_prefix}/{jira_ticket}` with analyses as subdirectories under the jira ticket directory.

For instance, library `A96199A` was analyzed with ticket `SC-1711`.  The alignment and hmmcopy results are stored under the directories `/work/shah/tantalus/SC-1711/results/results/alignment/` and `/work/shah/tantalus/SC-1711/results/results/hmmcopy_autoploidy/` respectively.  In Azure, the results are stored with blob prefixes `SC-1711/results/results/alignment/` and `SC-1711/results/results/hmmcopy_autoploidy/` respectively in the results container of the singlecellresults storage account.

### Dataset search

Given a library id, it will be necessary to search tantalus for a jira ticket corresponding to a QC analysis of that library in order to access the QC results.  This can be done with the `search_hmmcopy_analysis` function of the `scgenome.db.search` module.  

> Note that the logging setup is important in order to be able to obtain useful messages regarding caveats about the identified datasets, including whether they include all the available sequence data.

```
import sys
import logging
import scgenome.db.search

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

library_id = 'A96157C'

analysis = scgenome.db.search.search_hmmcopy_analysis(
    library_id,
    aligner_name='BWA_MEM_0_7_6A',
)

print(analysis['jira_ticket'])

>>>
SC-3041
```

### Loading from azure with caching

Given a list of jira tickets, it is possible to download data from azure blob storage to a local cache.  Repeated loads will not re-download.  Sample ids can be specified to filter for only a required list of samples.

> The environment variable `TANTALUS_CACHE_DIR` is assumed to have been set to the path to the local cache directory.

> You must have a tantalus account and access to azure blob storage.  Tantalus and azure blob credentials should be in your environment.

```
import sys
import logging
import scgenome.db.qc

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

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

results_tables = scgenome.db.qc.get_qc_data(
    hmmcopy_tickets,
    local_cache_directory,
    sample_ids=sample_ids,
    do_caching=True,
)

cn_data, metrics_data = (
    results_tables['hmmcopy_reads'],
    results_tables['annotation_metrics'],
)

print(cn_data.head())

>>> 
  chr    start      end  reads        gc      copy  state                cell_id sample_id library_id
0   1        1   500000     13 -1.000000       NaN      6  SA922-A90554B-R34-C70     SA922    A90554B
1   1   500001  1000000    442 -1.000000       NaN      6  SA922-A90554B-R34-C70     SA922    A90554B
2   1  1000001  1500000    461  0.598332  6.672340      6  SA922-A90554B-R34-C70     SA922    A90554B
3   1  1500001  2000000    478  0.539498  5.211916      6  SA922-A90554B-R34-C70     SA922    A90554B
4   1  2000001  2500000    594  0.594508  8.384862      6  SA922-A90554B-R34-C70     SA922    A90554B
```

### Loading from a storage directory

Given a list of jira tickets, it is also possible to load data directly from an existing storage, such as juno.  Simply set the storage directory to the directory containing the data and remember to leave do_caching as `False` (default).

```
import sys
import logging
import scgenome.db.qc

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

storage_directory = '/work/shah/tantalus/'

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

results_tables = scgenome.db.qc.get_qc_data(
    hmmcopy_tickets,
    storage_directory,
    sample_ids=sample_ids,
)

cn_data, metrics_data = (
    results_tables['hmmcopy_reads'],
    results_tables['annotation_metrics'],
)

print(cn_data.head())

>>> 
  chr    start      end  reads        gc      copy  state                cell_id sample_id library_id
0   1        1   500000     13 -1.000000       NaN      6  SA922-A90554B-R34-C70     SA922    A90554B
1   1   500001  1000000    442 -1.000000       NaN      6  SA922-A90554B-R34-C70     SA922    A90554B
2   1  1000001  1500000    461  0.598332  6.672340      6  SA922-A90554B-R34-C70     SA922    A90554B
3   1  1500001  2000000    478  0.539498  5.211916      6  SA922-A90554B-R34-C70     SA922    A90554B
4   1  2000001  2500000    594  0.594508  8.384862      6  SA922-A90554B-R34-C70     SA922    A90554B
```

### Retrieving pseudobulk data

Pseudobulk data is referenced by the jira ticket of the analysis that created that data.  As with QC data, the primary data store is in the `singlecellresults` azure blob storage.  A secondary store is on juno in the directory `/work/shah/tantalus/`.  Filesystem layout is `{storage_prefix}/{jira_ticket}` with analyses as subdirectories under the jira ticket directory.

For instance, the ovarian cell 2295 lines were analyzed with ticket `SC-1939`.  The pseudobulk results including SNVs, breakpoints, and haplotype allele results are stored under the directories `/work/shah/tantalus/SC-1939/results`.  In Azure, the results are stored with blob prefix `SC-1939/results` in the results container of the singlecellresults storage account.

Pseudobulk data can either be loaded directly from an existing server storage (eg `/work/shah/tantalus/`) or downloaded from the primary data store in azure and cached locally.

#### Caching pseudobulk data

To download from the primary data store, use the `datamanagement.transfer_files.cache_dataset` function from sisyphus.  The following code is searching for all results produced by an analysis with the jira ticket SC-1939, and then caching those results to the provided cache directory.

> The environment variable `TANTALUS_CACHE_DIR` is assumed to have been set to the path to the local cache directory.

> You must have a tantalus account and access to azure blob storage.  Tantalus and azure blob credentials should be in your environment.

```
import sys
import logging

import dbclients.tantalus
import datamanagement.transfer_files

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

ticket_id = 'SC-1939'

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

tantalus_api = dbclients.tantalus.TantalusApi()

ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

for results in ticket_results:
    filepaths = datamanagement.transfer_files.cache_dataset(
        tantalus_api,
        results['id'],
        'resultsdataset',
        'singlecellresults',
        local_cache_directory,
    )
```

#### Loading SNV data

The following code will load SNV tables for `SC-1939` from a local cache.

> Note that to load from a server storage simply replace `local_cache_directory` with the storage directory (eg `/work/shah/tantalus/`).

> Note that loading SNV data can be memory intensive

```
import os
import sys
import logging

import scgenome.loaders.snv

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

ticket_id = 'SC-1939'

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

ticket_directory = os.path.join(local_cache_directory, ticket_id)

snv_results = scgenome.loaders.snv.load_snv_data(
    ticket_directory,
)

print(snv_results.keys())
print(snv_results['snv_data'].head())

>>>
dict_keys(['snv_data', 'snv_count_data'])
   chrom    coord ref alt  alt_counts_sum  ref_counts_sum  mappability  \
17     1   985349   G   A            32.0            27.0          1.0   
24     1  1079129   G   T            80.0             0.0          1.0   
66     1  2032634   T   C            65.0            62.0          1.0   
68     1  2063033   C   A           354.0             1.0          1.0   
74     1  2117392   G   A           430.0             2.0          1.0   

   is_cosmic gene_name effect effect_impact amino_acid_change  \
17      True       NaN    NaN           NaN               NaN   
24       NaN       NaN    NaN           NaN               NaN   
66       NaN       NaN    NaN           NaN               NaN   
68       NaN       NaN    NaN           NaN               NaN   
74       NaN       NaN    NaN           NaN               NaN   

   tri_nucleotide_context  max_strelka_score  max_museq_score  
17                    CGT                 93             0.99  
24                    TGC                316             1.00  
66                    NaN                124             0.95  
68                    ACA                479             0.99  
74                    CGG                433             0.99  
```

#### Loading breakpoint data

The following code will load breakpoint tables for `SC-1939` from a local cache.  Note that to load from a server storage simply replace `local_cache_directory` with the storage directory (eg `/work/shah/tantalus/`).

```
import os
import sys
import logging

import scgenome.loaders.breakpoint

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

ticket_id = 'SC-1939'

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

ticket_directory = os.path.join(local_cache_directory, ticket_id)

breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data(
    ticket_directory,
)

print(breakpoint_results.keys())
print(breakpoint_results['breakpoint_data'].head())

>>>
dict_keys(['breakpoint_data', 'breakpoint_count_data'])
   prediction_id chromosome_1 strand_1  position_1 chromosome_2 strand_2  \
0            456            1        -    17084897            1        -   
1            569            1        -   234914910            1        -   
2           1032            1        -   148902774            1        -   
3           1313            1        -   204051209            1        -   
4           2081            1        -    64520958            1        -   

   position_2  homology  num_split inserted  ...  dgv_ids  is_germline  \
0    17084866         8          3      nan  ...      NaN        False   
1   234914883        13          6      nan  ...      NaN        False   
2   148902756         6          2      nan  ...      NaN        False   
3   204051202         0          2      nan  ...      NaN        False   
4    64520949         3          2      nan  ...      NaN        False   

   is_dgv  num_patients  is_filtered  dist_filtered  balanced  \
0   False             1        False         2769.0     False   
1   False             1        False         2695.0     False   
2   False             1        False         3302.0     False   
3   False             1        False        15009.0     False   
4   False             1        False        28133.0     False   

   rearrangement_type library_id  sample_id  
0            foldback    A90554A      SA921  
1            foldback    A90554A      SA921  
2            foldback    A90554A      SA921  
3            foldback    A90554A      SA921  
4            foldback    A90554A      SA921  

[5 rows x 37 columns]
    prediction_id                cell_id  read_count library_id sample_id
1               0  SA921-A90554A-R10-C22           1    A90554A     SA921
34              1  SA921-A90554A-R03-C08           1    A90554A     SA921
35              1  SA921-A90554A-R03-C19           1    A90554A     SA921
36              1  SA921-A90554A-R03-C27           1    A90554A     SA921
37              1  SA921-A90554A-R03-C36           1    A90554A     SA921
```

#### Loading haplotype data

The following code will load haplotype allele tables for `SC-1939` from a local cache.  Note that to load from a server storage simply replace `local_cache_directory` with the storage directory (eg `/work/shah/tantalus/`).

> Note that loading SNV data can be memory intensive

```
import os
import sys
import logging

import scgenome.loaders.allele

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

ticket_id = 'SC-1939'

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

ticket_directory = os.path.join(local_cache_directory, ticket_id)

allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
    ticket_directory,
)

print(allele_results.keys())
print(allele_results['allele_counts'].head())

>>>
dict_keys(['allele_counts'])
   allele_id                cell_id chromosome      end  hap_label  readcount  \
0          0  SA921-A90554A-R12-C09          1  1000000         27          1   
1          0  SA921-A90554A-R12-C09          1  2500000        151          1   
2          1  SA921-A90554A-R12-C09          1  3845268        259          1   
3          0  SA921-A90554A-R12-C09          1  4500000        285          1   
4          0  SA921-A90554A-R12-C09          1  7000000        406          1   

     start  
0   521368  
1  2000000  
2  3500000  
3  4000000  
4  6500000  
```
### Testing on ceto/juno

Run at the root of your repo directory, in your virtual environment with all required credentials available.

For more information on using the Juno High Performance Computing cluster, please see [here](https://shahcompbio.atlassian.net/wiki/spaces/SOFT/pages/28475398/Juno).

>test_load_qc
```
bsub -Is -R "rusage[mem=50]select[type==CentOS7]" python scgenome/tests/test_load_qc.py test-cached-single-ticket SC-2140 --local_storage_name juno
```
>test_load_pseudobulk
```
bsub -Is -R "rusage[mem=50]select[type==CentOS7]" python scgenome/tests/test_load_pseudobulk.py test-cached-single-ticket SC-2373 --local_storage_name juno
```
