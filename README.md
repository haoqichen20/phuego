<img src="https://github.com/haoqichen20/phuego/blob/master/docs/images/logo.png" width="400">

### phuEGO: a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets
---

phuego is a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets. It combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets.

To run phuego, the user can either directly calls the [CLI](#using-the-cli) application from command line, or integrate the [python functions](#using-the-python-package) into their own script. We recommend the user to use the CLI and submit jobs to a server where possible, to ensure robust execution of the tasks.

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

**Specific note for Windows users**: please provide the paths with forward slash '/' instead of backward slash '\'.

### 1. Downloading supporting dataset.
When using phuego **for the first time**, a support dataset that contains three zipped files (https://zenodo.org/record/8094690) need to be downloaded. After that, the user could always refer to the support data folder for running phuego. 

To download using the CLI application, the user could run the following. Make sure the target disk has at least **20 gb** of free space. 

```bash
# Download all three dataset, compare md5 checksum, unzip and removed zip files.
phuego -sf "path/to/desired/folder/" --need_fisher True --need_gic_sim True --need_networks True --remove_zip_file True

# If one of the file does not successfully download, for example the network file, then rerun above with only the missing file.
phuego -sf "path/to/desired/folder/" --need_networks True --remove_zip_file True
```

### 2. Performing a test/mock run with the phuego test dataset.
Phuego contains a test dataset, which is a list of differentially phosphorylated proteins between untreated cells and those treated with epidermal growth factor (EGF) (PMID: 19651622). To view the test dataset, see [here](#2-checking-the-input-data).

To understand the [input](#input) and [output](#output) of phuego, the user could either perform a mock run or a test run using the test dataset. A mock run performs 10 propagations on randomized networks, and thus isn't sufficient for reliable statistics, but typically finishs within a few minutes. A test run is a complete run with 1000 propagation on randomized networks, and typically takes > 1hr to finish.

```bash
# Performing a mock run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock

# Performing a mock run, with two kde_cutoff and two geneset database for geneset enrichment analysis.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock -k 0.85 -k 0.9 -fg "K" -fg "B"

# Performing a test run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_test
```

### 3. Running your own protein list, reusing rwr results.
To run phuego using your own list of protein, provide the input file path to phuego.

**Important:** the network propagation (rwr) step is the most time consuming step of phuego. For each combination of protein list and damping factor, a separate propagation would be run and the result will be stored in the provided result folder (pvalues.txt, rwr_scores.txt, start_seeds.txt). After this, the user can reuse the results (so make sure you don't delete them!) and test other parameters by providing flag **-ru**.

```bash
# Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -d 0.85 -k 0.85 -fg "K"

# Performing a run reusing network propagation result, testing two different KDE values, and two different gene sets, and export the network in a different format.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ru -k 0.8 -k 0.9 -fg "C" -fg "B" -nf "edgelist"
```

### 4. Batch job submission (LSF cluster).
Each phuego run works with one protein list and one damping factor. If you have multiple protein lists (e.g., from a set of experiment), and/or would like to test multiple damping factors, you could use a .sh script to call phuego multiple times, and submit jobs in batch manner. Below we provide a .sh script for a LSF cluster as an example.

To do so, first create a test_datasets.txt file that store the path to all your protein list files:

```text
path/to/protein_list_1.txt
path/to/protein_list_2.txt
path/to/protein_list_3.txt
```

Then submit your jobs using the following .sh script. The output will be organized into a two-layer folder structure under your specified result_dir. Modify argument value to suit your need.

```bash
#!/bin/bash

# Path to the dataset file
dataset_file="path/to/test_datasets.txt"

# Read dataset names from the file into an array
readarray -t datasets < "$dataset_file"

# Result folder.
result_dir="path/to/result_dir"

# Run phuego.
i=0
for line in "${datasets[@]}"; do
    # Create numerical job name.
    job_name="job_$((i+1))"
    ((i++))

    # Extract the last level from the input.
    experiment=$(echo "$line" | rev | cut -d'/' -f1 | rev)
    
    # Create the experiment dir.
    exp_dir="$pub_dir/$experiment"

    dampings=(0.5 0.7 0.85)
    for damping in "${dampings[@]}"; do
        # Create the damping dir.
        damping_dir="$exp_dir/$damping"
        mkdir -p $damping_dir

        # Run phuego with this damping factor, intact background and two kde_cutoff. 
        # Run with 4 cores to increase the speed.
        bsub -n 4 -M 4096 -R "rusage[mem=4096]" -o log.txt -e err.txt -J "$job_name"\
        phuego -sf "Path/to/support_data/" -rf "$damping_dir" -tpath "$line" \
        -d $damping -k 0.85 -k 0.9 -fg "B" -fg "K" -nf "graphml"
    done
done
```

### 5. Running phuego with removed network nodes.
In a drugging or a knockout experiment, one might want to removed the knocked out targets from the reference network before performing network propagation, assuming that they are no longer functional. To do so, one could specify a .csv file as below, and provide to phuego:

```text
UniprotID_1,UniprotID_2,UniprotID_3
UniprotID_1,UniprotID_2,UniprotID_3
```

Here, row 1 is a list of targets to be removed from the network propagation of upregulated input proteins, and row 2 for downregulated. Normally, one would expect these to be the same. The list can be provided as following:

```bash
# Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ipath "path/to/targets_list.csv" -d 0.85 -k 0.85 -fg "K" 
```

## Using the python package

The functions in the phuego package can be imported and used in your own python scripts. This makes it easier for integrating phuego into your own workflow.

### 1. Downloading supporting dataset

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

### 2. Checking the input data
```python
from phuego import load_test_example

# The function return the absolute path of the test dataset, and the dataset itself as a dataframe.
test_path, test_df = load_test_example()
print(test_df)
```

### 3. Performing a test/mock run with the phuego test dataset

Performing a mock run.

```python
from phuego import phuego_mock

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
convert2folder=True
include_isolated_egos_in_KDE_net=False
net_format="graphml"

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
    convert2folder=convert2folder,
    include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
    net_format=net_format,
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
rwr_threshold = 0.05
kde_cutoff = [0.85]
use_existing_rwr = False
convert2folder=True
include_isolated_egos_in_KDE_net=False
net_format="graphml"

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
    rwr_threshold=rwr_threshold,
    kde_cutoff=kde_cutoff,
    use_existing_rwr=use_existing_rwr,
    convert2folder=convert2folder,
    include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
    net_format=net_format,
    )
```

### 4. Running your own protein list
To run phuego on your own protein list, simply provide the **test_path** to the above code, and remove the **load_test_example()** line. 

To reuse network propagation result and explore different KDE cutoff / genesets, set **use_existing_rwr = True** (also see above).


## Input and results.
### Supporting dataset. 
The supporting datasets contain the base network, randomized networks, semantic similarity, geneset database etc. For details, refer to the associated [publication](#citation).

### Input
phuEGO accept in input a .txt file where first column are uniprot Ids and second columns are LFC or values associated to the importance of the protein in the first column, such as the following. The user can also provide a [.csv](#5-running-phuego-with-removed-network-nodes) file for removed nodes.

```text
P29317	2.1043
P00533	1.36255
P10398	0.7655
Q9UHY1	0.585
P53999	-1.3219
Q96II8	-1.3219
P17096	-1.0
P32519	-1.0
```

### Output, overview:
Phuego output numerous files generated at each steps of the algorithm, organized into a folder structure. If the user prefer to navigate the files without the folder structure, they could set the flag **--dont_convert2folder**.

In brief, phuego processed the user input in two regulation directions: 

    - increased
    - decreased

In each direction, there are four levels of phuego output: 

    - seeds
    - rwr
    - ego
    - modules

#### Output: networks

For most users, the most interesting result would be the **module network files**. So we will explain these first:

#### module_net.graphml
This is a network output of the igraph package, in an user specified format (**-nf**). It contains information for module annotation, and can be readily imported into other software such as Cytoscape.

#### module_net_edgelist.csv

#### module_net_nodes_attribute.csv


### Output: GSEA
#### Fisher test
Inside the fishers folder 8 files are can be found:
Pfisher.txt refers to the enriched terms against Gene Ontology biological process
Ffisher.txt refers to the enriched terms against Gene Ontology functional 
Cfisher.txt refers to the enriched terms against Gene Ontology cellular component
Kfisher.txt refers to the enriched terms against KEGG
Rfisher.txt refers to the enriched terms against Reactome when only the leaves are consider as annotation
RTfisher.txt refers to the enriched terms against Reactome when all the hierarchy is considered
Dfisher.txt refers to the enriched terms against DisGenenet
Bfisher.txt refers to the enriched terms against Bioplanet

### Other output



#### pvalues.txt
This file contains the pvalues for each node of the network. It is divided in seven columns, the first column refers to the uniprot ids.
Columns from 2,3,and 4 refers to the pvalues associated with the increased phosphorylation nodes pvalues, of which:
    -second column refers to the pvalues when increased phosphorilated tyrosine are used as seed nodes 
    -third column refers to the pvalues when all the other increased phosphorilated kinases are used as seed nodes, 
    -fourth column refers to the pvalues when the increased phosphorilated substrates are used as seed nodes.

Columns from 5,6 and 7 refers to the pvalues associated with the decreased phosphorylated nodes pvalues, of which:
    -fifth column refers to the pvalues when decreased phosphorilated tyrosine are used as seed nodes 
    -sixt column refers to the pvalues when all the other decreased phosphorilated kinases are used as seed nodes, 
    -seventh column refers to the pvalues when the decreased phosphorilated substrates are used as seed nodes.

A value greater than 950 indicates a pvalues<0.05 as well as a values greater than 990 indicates a pvalues<0.01
#### rwr_scores.txt
This file has the same format of pvalues.txt with the difference that values indicates rwr scores 

#### start_seeds.txt 
First column is uniprot id, second column is phosphorylation LFC (is it exactly the same as the user input?)

#### KDE_egos.txt
First column is a seed node, the remaining column are the neighbors associated with the seed nodes. 

#### module_egos.txt
Has the same format of KDE_egos.txt but refers to the module specific nodes associated to the seed nodes.



## Development

phuego can be integrated into any phosphoproteomics analysis pipeline. It will soon become available through NF-core/Nextflow for easier pipeline integration.

## Citation

Please cite phuego if you use it in your analysis.

```BibTeX
```

## Contributors

The algorithm and initial scripts of phuego are developed by Girolamo Giudice ([@girolamogiudice](https://github.com/girolamogiudice)) and Evangelia Petsalaki at [@EMBL-EBI] (https://www.ebi.ac.uk/).

The Python package and CLI application are developed by Haoqi Chen ([@haoqichen20] (https://github.com/haoqichen20) at [@EMBL-EBI] (https://www.ebi.ac.uk/)).