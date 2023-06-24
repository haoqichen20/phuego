

### phuEGO: a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets
---

A general description of phuego.

## Installation

```bash
pip install phuego
```

## Usage

When using phuego for the first time, a support dataset that include semantic similarity, randomized Omnipath networks and multiple genesets need to be downloaded. The support dataset are hosted on (https://zenodo.org/), and the package provided a function for downloading the dataset.



```python
from phuego import dataprep

# Do not omit the forward slash at the end.
support_data_folder = "path/to/support_data_folder/"

# Download support dataset. If the user wish the keep the zip files for easier 
# data transfer between different local user/folders, change remove_zip_file to 
# False.
dataprep(support_data_folder=support_data_folder, 
            need_fisher=True, 
            need_gic_sim=True, 
            need_networks=True,
            remove_zip_file=True)
```

After downloading, perform a test run with the phuego test dataset, to familirize with the input data and output results.

```python

```

To run phuego on a single dataset, use the following code by changing the settings. Network propagation (using random walk with restart algorithm, rwr) is the most time consuming step of phuego. Therefore, when this step is finished, phuego stored the results in res_folder. User can set | use_existing_rwr = True | to reuse the network propagation results, and explore different settings of kde_cutoff or fisher genesets.

Note that for each damping factor, phuego would run a separate network propagation process.

```python
from phuego import phuego

# User input: paths.
support_data_folder = "/nfs/research/petsalaki/users/hchen/ph_phuego_test/support_data/"
res_folder = "/nfs/research/petsalaki/users/hchen/ph_phuego_test/result_10/package/"
test_path = "/nfs/research/petsalaki/users/hchen/pypi_phuego/phuego/data/EGF_vs_untreated_short.txt"

# User input: algorithm setting.
ini_pos = 'False'
ini_neg = 'False'
damping = 0.85
fisher_geneset = ["C","F","D","P","R","K","RT","B"]
fisher_threshold = 0.05
use_existing_rwr = False
kde_cutoff = [0.85, 0.9]

# Run phuego.
number_of_nodes, number_of_genes = phuego(
        support_data_folder=support_data_folder,
        res_folder=res_folder,
        test_path=test_path,
        fisher_geneset=fisher_geneset,
        fisher_threshold=fisher_threshold,
        ini_pos=ini_pos,
        ini_neg=ini_neg,
        damping=damping,
        kde_cutoff=kde_cutoff,
        use_existing_rwr=use_existing_rwr,
        )
```



## Development

phuego can be integrated into any phosphoproteomics analysis pipeline. It will soon become available through NF-core/Nextflow for easier pipeline integration.

## Citation

Please cite phuego if you use it in your analysis.

```BibTeX
```

## Contributors

The algorithm and python functions of phuego is developed by Girolamo Giudice ([@girolamogiudice](https://github.com/girolamogiudice)) and Evangelia Petsalaki at [EMBL-EBI] (https://www.ebi.ac.uk/).

The Python package and application is developed by Haoqi Chen ([@haoqichen20] (https://github.com/haoqichen20) at [EMBL-EBI] (https://www.ebi.ac.uk/)).



