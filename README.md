

### phuEGO: a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets
---

phuego is a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets. It combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets.

To run phuego, the user can either integrate the python functions into their own script, or directly calls the CLI application from command line. Both methods work equivalently, and it's up to the user to decide which method suits their work.

## Installation

```bash
# Install phuego in your python environment with python > 3.9
pip install phuego
phuego --version
phuego --help
```

## Using the CLI

The CLI application is an easy way to run phuego on single or a whole batch of dataset directly from the command line. 

### 1. Downloading supporting dataset.
When using phuego for the **first time**, a support dataset that contains three files (https://zenodo.org/record/8094690) need to be downloaded. After that, the user could always refer to the support data folder for running phuego. To download using the CLI application, the user could run the following. Make sure the target disk has at least **20 gb** of free space. 

```bash
# Download all three dataset, compare md5 checksum, unzip and removed zip files.
phuego -sf "path/to/desired/folder/" --need_fisher True --need_gic_sim True --need_networks True --remove_zip_file True
# If one of the file does not successfully download, for example the network file, then rerun above with only the missing file.
phuego -sf "path/to/desired/folder/" --need_networks True --remove_zip_file True
```


### 2. Performing a test/mock run with the phuego test dataset.
To test whether 
### Running your own dataset. 

### Submitting batch job to a lsf server.






## Using the python package.


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

The Python package and CLI application is developed by Haoqi Chen ([@haoqichen20] (https://github.com/haoqichen20) at [EMBL-EBI] (https://www.ebi.ac.uk/)).



