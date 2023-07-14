<img src="https://github.com/haoqichen20/phuego/blob/master/docs/images/logo.png" width="400">

### phuEGO: a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets
---

phuego is a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets. It combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets.

To run phuego, the user can either directly calls the CLI application from command line, or integrate the python functions into their own script. Both methods work equivalently, and it's up to the user to decide which method suits their work.

## Installation

```bash
# Install phuego in your python environment with 3.13 > python > 3.9
pip install phuego
# Check version to ensure successful installation.
phuego --version
# Check all parameters of CLI application.
phuego --help
```

## Using the CLI

The CLI application is an easy way to run phuego on single or a whole batch of dataset directly from the command line. 

### 1. Downloading supporting dataset.
When using phuego **for the first time**, a support dataset that contains three files (https://zenodo.org/record/8094690) need to be downloaded. After that, the user could always refer to the support data folder for running phuego. To download using the CLI application, the user could run the following. Make sure the target disk has at least **20 gb** of free space. 

```bash
# Download all three dataset, compare md5 checksum, unzip and removed zip files.
phuego -sf "path/to/desired/folder/" --need_fisher True --need_gic_sim True --need_networks True --remove_zip_file True

# If one of the file does not successfully download, for example the network file, then rerun above with only the missing file.
phuego -sf "path/to/desired/folder/" --need_networks True --remove_zip_file True
```

### 2. Performing a test/mock run with the phuego test dataset.
phuego contains a test dataset, which is a list of differentially phosphorylated proteins between untreated cells and those treated with epidermal growth factor (EGF) (PMID: 19651622). 

To understand the input and output of phuego, the user could perform a test run or a mock run using the test dataset. A test run is a complete run with 1000 propagation on randomized networks, and typically takes > 1hr to finish. A mock run performs 10 propagations and thus isn't sufficient for reliable statistical testing, but it typically finish within a few minutes and can help understand the format of input and output. 

Both test will print the head of the test dataset to standard output (e.g., your screen). To view the full test dataset, see next session "Using the python package: understanding the input data."

```bash
# Performing a mock run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -mock True

# Performing a test run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -test True
```

### 3. Running your own protein list.
To run phuego using your own list of protein, the data should be formatted in the same way as the test dataset. For detailed explanation of parameters, see section "Further documentations".

**Important:** the network propagation step (rwr, random walk with restart) is the most time consuming step of phuego. For each protein list and each damping factor, a separate propagation would be run and the result will be stored in the provided result folder (pvalues.txt, rwr_scores.txt, start_seeds.txt). After this, the user can reuse the results (so make sure you don't delete them!) to test different KDE cutoff and perform gene set enrichment analysis with various geneset by setting **-ru True**.

```bash
# Performing a run for the first time with example parameters. 
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ru False -d 0.85 -k 0.85 -fg "B"

# Performing a run reusing network propagation result, testing two different KDE values, and two different gene sets.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ru True -k 0.8 -k 0.9 -fg "C" -fg "K"
```

### 4. (optional) Batch job submission (LSF cluster).
Each phuego run works with one protein list and one damping factor. If you have multiple protein lists (e.g., from a set of experiment), and/or would like to test multiple damping factors, you could use a .sh script to call phuego multiple times, and submit jobs in batch manner. Below we provide a .sh script for a LSF cluster as an example.

To do so, first create a test_datasets.txt file that store the path to all your protein list files:

```text
path/to/protein_list_1.txt
path/to/protein_list_2.txt
path/to/protein_list_3.txt
```

Then submit your jobs using the following .sh script:

```bash
#!/bin/bash

dataset_file="path/to/test_datasets.txt"  # Path to the dataset file

# Read dataset names from the file into an array
readarray -t datasets < "$dataset_file"

# Iterate over each dataset and submit a job
for (( i=0; i<${#datasets[@]}; i++ )); do
    dataset="${datasets[i]}"
    job_name="job_$((i+1))"  # Use numeric index as the job name
    result_dir="./result/$job_name/"

    # Create the result directory for the current job
    mkdir -p "$result_dir"
    
    # Submit the job using bsub with the desired job configuration and dataset
    bsub -q standard -n 4 -M 4096 -R "rusage[mem=4096]" -o log_ph.txt -e err_ph.txt -J "$job_name" phuego -sf "./support_data/" -rf "$result_dir" -ru False -tpath "$dataset" -d 0.85 -k 0.85 -k 0.9 
done
```



## Using the python package.

The functions in the phuego package can be imported and used in your own python scripts. This makes it easier for integrating phuego into your own workflow.

### 1. Downloading supporting dataset.

```python
from phuego import dataprep

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

### 2. Understanding the input data.
```python
from phuego import load_test_example

# The function return the absolute path of the test dataset, and the dataset itself as a dataframe.
test_path, test_df = load_test_example()
print(test_df)
```

### 3. Performing a test/mock run with the phuego test dataset.

Performing a mock run.

```python
from phuego import phuego_mock

# User input: paths.
support_data_folder = "path/to/support_data_folder/"
res_folder = "path/to/desired_result_folder/"

# Loading the test dataset.
test_path, test_df = load_test_example()

# When calling the functions, the user need to manually set all parameters.
fisher_geneset = ["C","F","D","P","R","K","RT","B"]
fisher_threshold = 0.05
fisher_background = "intact"
ini_pos = ["False"]
ini_neg = ["False"]
damping = 0.85
kde_cutoff = [0.85]
use_existing_rwr = False

# Run phuego mock.
print("Run phuego_mock with test dataset, whose first few lines are: \n",test_df.head())
phuego_mock(
    support_data_folder=support_data_folder,
    res_folder=res_folder,
    test_path=test_path,
    fisher_geneset=fisher_geneset,
    fisher_threshold=fisher_threshold,
    fisher_background=fisher_background,
    ini_pos=ini_pos,
    ini_neg=ini_neg,
    damping=damping,
    kde_cutoff=kde_cutoff,
    use_existing_rwr=use_existing_rwr,
    )
```

Performing a full test run.

```python
from phuego import phuego

# User input: paths.
support_data_folder = "path/to/support_data_folder/"
res_folder = "path/to/desired_result_folder/"

# Loading the test dataset.
test_path, test_df = load_test_example()

# When calling the functions, the user need to manually set all parameters.
fisher_geneset = ["B"]
fisher_threshold = 0.05
fisher_background = "intact"
ini_pos = ["False"]
ini_neg = ["False"]
damping = 0.85
kde_cutoff = [0.85]
use_existing_rwr = False

# Run phuego with test dataset..
print("Run phuego with test dataset, whose first few lines are: \n",test_df.head())
phuego(
    support_data_folder=support_data_folder,
    res_folder=res_folder,
    test_path=test_path,
    fisher_geneset=fisher_geneset,
    fisher_threshold=fisher_threshold,
    fisher_background=fisher_background,
    ini_pos=ini_pos,
    ini_neg=ini_neg,
    damping=damping,
    kde_cutoff=kde_cutoff,
    use_existing_rwr=use_existing_rwr,
    )
```

### 4. Running your own protein list.
To run phuego on your own protein list, simply provide the **test_path** to the above code, and remove the **load_test_example()** line. 

To reuse network propagation result and explore different KDE cutoff / genesets, set **use_existing_rwr = True** (also see above).

## Detailed documentations.
### Supporting dataset. 
#### gic_sim folder, gic_sim_std.txt
some text
#### networks folder
some text
#### proteome and intact
some text
#### uniprot_to_gene.tab
some text
#### pfam_domains.txt
some text

### phuego function variables
#### iteam A
#### iteam B

### Additional phuego CLI variables
#### iteam A
#### iteam B

### Input
#### iteam A
#### iteam B

### Output
#### iteam A
#### iteam B


## Development

phuego can be integrated into any phosphoproteomics analysis pipeline. It will soon become available through NF-core/Nextflow for easier pipeline integration.

## Citation

Please cite phuego if you use it in your analysis.

```BibTeX
```

## Contributors

The algorithm and initial scripts of phuego are developed by Girolamo Giudice ([@girolamogiudice](https://github.com/girolamogiudice)) and Evangelia Petsalaki at [@EMBL-EBI] (https://www.ebi.ac.uk/).

The Python package and CLI application are developed by Haoqi Chen ([@haoqichen20] (https://github.com/haoqichen20) at [@EMBL-EBI] (https://www.ebi.ac.uk/)).