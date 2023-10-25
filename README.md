<img src="https://github.com/haoqichen20/phuego/blob/master/docs/images/logo.png" width="400">

### phuEGO: a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets
---

phuEGO is a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets. It combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets.

To run phuEGO, the user can either directly calls the [CLI](#2-using-the-cli) application from command line, or integrate the [python functions](#3-using-the-python-package) into their own script. We recommend the user to use the CLI and submit jobs to a server where possible, to ensure robust execution of the tasks.

## 1. Installation

```bash
# Install phuEGO in your python environment with 3.13 > python > 3.9
pip install phuego
# Check version to ensure successful installation.
phuego --version
# Check all parameters of CLI application.
phuego --help
```

## 2. Using the CLI

The CLI application is an easy way to run phuEGO on single or a whole batch of dataset directly from the command line. 

**Specific note for Windows users**: please provide the paths with forward slash '/' instead of backward slash '\\'.

### 1). Downloading supporting dataset
When using phuEGO **for the first time**, a support dataset that contains three zipped files (https://zenodo.org/record/8094690) need to be downloaded. We recommend using the zenodo-get packages:

```bash
# Install zenodo_get to download the supporting dataset from Zenodo.
pip install zenodo-get
cd path/to/desired/folder
zenodo_get 8094690
```

Decompressed the files into the support datafolder using the following Python code or other tools, such as tar. **Do not modify the folder structure.**

```python
import tarfile

support_data_folder = "path/to/desired/folder" #Change this.

def decompress_tar_gz(file_path, output_dir):
    with tarfile.open(file_path, "r:gz") as tar:
        tar.extractall(path=output_dir)       
        
paths = ["./gic_sim.tar.gz", "./fisher_etc.tar.gz", "./networks.tar.gz"]
for path in paths:
    decompress_tar_gz(path, support_data_folder)
```

### 2). Performing a test/mock run with the phuEGO test dataset
Phuego contains a test dataset, which is a list of differentially phosphorylated proteins between untreated cells and those treated with epidermal growth factor (EGF) (PMID: 19651622). To view the test dataset, see [here](#2-checking-the-input-data).

To understand the [input](#2-input) and [output](#3-output) of phuEGO, the user could either perform a mock run or a test run using the test dataset. A mock run performs 10 propagations on randomized networks, and thus isn't sufficient for reliable statistics, but typically finishs within a few minutes. A test run is a complete run with 1000 propagation on randomized networks, and typically takes > 1hr to finish.

```bash
# Performing a mock run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock

# Performing a mock run, with two kde_cutoff and two geneset database for geneset enrichment analysis.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock -k 0.85 -k 0.9 -fg "K" -fg "B"

# Performing a test run.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_test
```

### 3). Running your own protein list, reusing rwr results
To run phuEGO using your own list of protein, provide the input file path to phuEGO.

**Important:** the network propagation (rwr) step is the most time consuming step of phuEGO. For each combination of protein list and damping factor, a separate propagation would be run and the result will be stored in the provided result folder (pvalues.txt, rwr_scores.txt, start_seeds.txt). After this, the user can reuse the results (so make sure you don't delete them!) and test other parameters by providing flag **-ru**.

```bash
# Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -d 0.85 -k 0.85 -fg "K"

# Performing a run reusing network propagation result, testing two different KDE values, and two different gene sets, and export the network in a different format.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ru -k 0.8 -k 0.9 -fg "C" -fg "B" -nf "edgelist"
```

### 4). Batch job submission (LSF cluster)
Each phuEGO run works with one protein list and one damping factor. If you have multiple protein lists (e.g., from a set of experiment), and/or would like to test multiple damping factors, you could use a .sh script to call phuEGO multiple times, and submit jobs in batch manner. Below we provide a .sh script for a LSF cluster as an example.

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

# Run phuEGO.
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

        # Run phuEGO with this damping factor, intact background and two kde_cutoff. 
        # Run with 4 cores to increase the speed.
        bsub -n 4 -M 4096 -R "rusage[mem=4096]" -o log.txt -e err.txt -J "$job_name"\
        phuego -sf "Path/to/support_data/" -rf "$damping_dir" -tpath "$line" \
        -d $damping -k 0.85 -k 0.9 -fg "B" -fg "K" -nf "graphml"
    done
done
```

### 5). Running phuEGO with removed network nodes
In a drugging or a knockout experiment, one might want to removed the knocked out targets from the reference network before performing network propagation, assuming that they are no longer functional. To do so, one could specify a .csv file as below, and provide to phuEGO:

```text
UniprotID_1,UniprotID_2,UniprotID_3
UniprotID_1,UniprotID_2,UniprotID_3
```

Here, row 1 is a list of targets to be removed from the network propagation of upregulated input proteins, and row 2 for downregulated. Normally, one would expect these to be the same. The list can be provided as following:

```bash
# Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ipath "path/to/targets_list.csv" -d 0.85 -k 0.85 -fg "K" 
```

## 3. Using the python package

The functions in the phuEGO package can be imported and used in your own python scripts. This makes it easier for integrating phuEGO into your own workflow.

### 1). Downloading supporting dataset

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

### 2). Checking the input data
```python
from phuego import load_test_example

# The function return the absolute path of the test dataset, and the dataset itself as a dataframe.
test_path, test_df = load_test_example()
print(test_df)
```

### 3). Performing a test/mock run with the phuEGO test dataset

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

# Run phuEGO mock.
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

# Run phuEGO with test dataset..
print("Run phuEGO with test dataset, whose first few lines are: \n",test_df.head())
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

### 4). Running your own protein list
To run phuEGO on your own protein list, simply provide the **test_path** to the above code, and remove the **load_test_example()** line. 

To reuse network propagation result and explore different KDE cutoff / genesets, set **use_existing_rwr = True** (also see above).


## 4. Data and results
### 1). Supporting dataset. 
The supporting datasets contain the base network, randomized networks, semantic similarity, geneset database etc. For details, refer to the associated [publication](#6-citation).

### 2). Input
phuEGO accept in input a .txt file where first column are uniprot Ids and second columns are LFC or values associated to the importance of the protein in the first column, such as the following. The user can also provide a .csv file for removed nodes, as explained [above](#5-running-phuego-with-removed-network-nodes).

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

### 3). Output - overview
Phuego output numerous files generated at each steps of the algorithm, organized into a folder structure. If the user prefer to navigate the files without the folder structure, they could set the flag **--dont_convert2folder** when using the CLI, or set **convert2folder=False** when using the python functions.

phuEGO processed the user input in two regulation directions: 

    - increased
    - decreased

In each direction, there are four levels of phuEGO output: 

    - seeds
    - rwr
    - ego (subjected to KDE_cutoff)
    - modules (subjected to KDE_cutoff)

Depending on the levels, three types of files could be provided:

    - protein lists
    - networks
    - fisher's exact test results

In brief, a user input protein list (seeds) is divided into two directions using the scores (e.g., log2FCs). On each direction, it is expanded through network propagation into a very large list (rwr). Then this list is narrowed down (depending on where the user set the KDE_cutoff) by extracting only Ego network of seed nodes (ego), and further filtered by keeping inter-linked egos (modules). **The [modules](#4-output---module-networks) are considered the final result of phuEGO.** For details, refer to [publication](#6-citation).

### 4). Output - module networks

Module networks, generated by the Python igraph package, in an user specified format (via CLI parameter **-nf**). The default .graphml format contains information for module annotation, and can be readily imported into other software such as Cytoscape.

    - decreased/KDE_cutoff/networks/module_net.graphml
    - increased/KDE_cutoff/networks/module_net.graphml

Files contains redundant information of the module networks in a more human-readable format.

    - decreased/KDE_cutoff/networks/module_net_edgelist.csv
    - decreased/KDE_cutoff/networks/module_net_nodes_attribute.csv
    - increased/KDE_cutoff/networks/module_net_edgelist.csv
    - increased/KDE_cutoff/networks/module_net_nodes_attribute.csv

The network of the ego is also output at: 

    - decreased/KDE_cutoff/networks/KDE.graphml
    - increased/KDE_cutoff/networks/KDE.graphml

### 5). Output - GSEA
On each level, geneset enrichment analysis (fisher's exact test) are performed on the protein list, against an user specified geneset database. Through this, the user can examine whether phuEGO retrieves a more meaningful protein list from the seeds. The geneset database are part of the supporting database. In total, 8 databases are included, and one or more databases can be specified by providing abbreviation to the CLI arguments **--fg**. 

    -P: enrichment against Gene Ontology biological process. Output file: Pfisher.txt
    -F: enrichment against Gene Ontology functional . Output file: Ffisher.txt
    -C: enrichment against Gene Ontology cellular component. Output file: Cfisher.txt
    -K: enrichment against KEGG. Output file: Kfisher.txt
    -R: enrichment against Reactome when only the leaves are consider as annotation. Output file: Rfisher.txt
    -R: enrichment against Reactome when all the hierarchy is considered. Output file: Rfisher.txt
    -D: enrichment against DisGenenet. Output file: Dfisher.txt
    -B: enrichment against Bioplanet. Output file: Bfisher.txt

They are distributed across the following paths:

    - seeds: decreased/seed_fisher; increased/seed_fisher
    - rwr: decreased/rwr_fisher; increased/rwr_fisher
    - ego: decreased/KDE_cutoff/fisher; increased/KDE_cutoff/fisher
    - modules: decreased/KDE_cutoff/modules/module_N/fisher; increased/KDE_cutoff/modules/module_N/fisher

### 6). Output - protein list
On the ego level, a text file summarize the protein list within ego network of each seed on a separate row. The first column is a seed node, the remaining column are the neighbors associated with the seed nodes. If **-ie** is provided to CLI, then rows with only the seeds will be excluded during network generation.

    - decreased/KDE_cutoff/KDE_egos.txt
    - increased/KDE_cutoff/KDE_egos.txt

On the module level, for each module, a similar text file is generated, containing only module-specific seeds and the associated nodes.

    - decreased/KDE_cutoff/modules/module_N/module_egos.txt
    - increased/KDE_cutoff/modules/module_N/module_egos.txt

### 7). Output - network propagation results
Three files are stored on the top layers of the result folder. These files can be reused by phuEGO to save time, when the user provide the flag **-ru** to CLI.

    - pvalues.txt
    - rwr_scores.txt
    - start_seeds.txt

The pvalues.txt file contains the pvalues for each node of the network. It is divided in seven columns, the first column refers to the uniprot ids.

Columns from 2,3,and 4 refers to the pvalues associated with the increased phosphorylation nodes pvalues, of which:
    -second column refers to the pvalues when increased phosphorilated tyrosine are used as seed nodes 
    -third column refers to the pvalues when all the other increased phosphorilated kinases are used as seed nodes, 
    -fourth column refers to the pvalues when the increased phosphorilated substrates are used as seed nodes.

Columns from 5,6 and 7 refers to the pvalues associated with the decreased phosphorylated nodes pvalues, of which:
    -fifth column refers to the pvalues when decreased phosphorilated tyrosine are used as seed nodes 
    -sixt column refers to the pvalues when all the other decreased phosphorilated kinases are used as seed nodes, 
    -seventh column refers to the pvalues when the decreased phosphorilated substrates are used as seed nodes.

A value greater than 950 indicates a pvalues<0.05 as well as a values greater than 990 indicates a pvalues<0.01. 

The rwr_scores.txt file has the same format of pvalues.txt with the difference that values indicates rwr scores.

The start_seeds.txt is basically the same as the user input.

## 5. Development

### Future development

phuEGO can be integrated into any phosphoproteomics analysis pipeline. It will eventually become available through NF-core/Nextflow for easier pipeline integration.

### Changelog

For a detailed history of changes, see the [Changelog](https://github.com/haoqichen20/phuego/blob/master/CHANGELOG.md).

## 6. Citation

Please cite phuEGO if you use it in your analysis.

> **phuEGO: A network-based method to reconstruct active signalling pathways from phosphoproteomics datasets** <br> _Girolamo Giudice, Haoqi Chen, Evangelia Petsalaki_ <br>
> bioRxiv, 2023 <br>
> doi: [10.1101/2023.08.07.552249](https://doi.org/10.1101/2023.08.07.552249) <br>

```BibTeX
@article{giudice_phuego_2023,
	title = {{phuEGO}: {A} network-based method to reconstruct active signalling pathways from phosphoproteomics datasets},
	url = {https://www.biorxiv.org/content/early/2023/08/07/2023.08.07.552249},
	doi = {10.1101/2023.08.07.552249},
	abstract = {Signalling networks are critical for virtually all cell functions. Our current knowledge of cell signalling has been summarised in signalling pathway databases, which, while useful, are highly biassed towards well-studied processes, and don’t capture context specific network wiring or pathway cross-talk. Mass spectrometry-based phosphoproteomics data can provide a more unbiased view of active cell signalling processes in a given context, however, it suffers from low signal-to-noise ratio and poor reproducibility across experiments. Methods to extract active signalling signatures from such data struggle to produce unbiased and interpretable networks that can be used for hypothesis generation and designing downstream experiments. Here we present phuEGO, which combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets. We applied phuEGO to five phosphoproteomics data sets from cell lines collected upon infection with SARS CoV2. PhuEGO was better able to identify common active functions across datasets and to point to a subnetwork enriched for known COVID-19 targets. Overall, phuEGO provides a tool to the community for the improved functional interpretation of global phosphoproteomics datasets.Competing Interest StatementThe authors have declared no competing interest.},
	journal = {bioRxiv},
	author = {Giudice, Girolamo and Chen, Haoqi and Petsalaki, Evangelia},
	year = {2023},
	note = {Publisher: Cold Spring Harbor Laboratory
\_eprint: https://www.biorxiv.org/content/early/2023/08/07/2023.08.07.552249.full.pdf},
}
```

## 7. Contributors

The algorithm and initial scripts of phuEGO are developed by Girolamo Giudice ([@girolamogiudice](https://github.com/girolamogiudice)) and Evangelia Petsalaki at [EMBL-EBI](https://www.ebi.ac.uk/).

The Python package and CLI application are developed by Haoqi Chen ([@haoqichen20](https://github.com/haoqichen20)) at [EMBL-EBI](https://www.ebi.ac.uk/).