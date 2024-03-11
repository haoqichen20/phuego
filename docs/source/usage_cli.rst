Running phuEGO
==============

The main analysis of phuEGO is implemented by the CLI. 
Three paths should be provided to the CLI, as well as a number of arguments.
These can be either directly provided using command line flags, or aggregated
in a text file (hence **argument file**) for clarity. 
To get an overview of the CLI arguments, use:

.. code-block:: bash

   phuego --help

.. warning::

   For Windows users, please provide the paths with forward slash ‘/’ 
   instead of backward slash ‘\\’.


Running for the first time
~~~~~~~~~~~~~~~~~~~~~~~~~~

To run phuEGO using your own list of protein, provide the input file
path to phuEGO.

.. code:: bash

   # Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -d 0.85 -k 0.85 -fg "K"

   # Performing a run reusing network propagation result, testing two different KDE values, and two different gene sets, and export the network in a different format.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ru -k 0.8 -k 0.9 -fg "C" -fg "B" -nf "edgelist"


Reusing pagerank results
~~~~~~~~~~~~~~~~~~~~~~~~

**Important:** the network propagation (rwr) step is the most time
consuming step of phuEGO. For each combination of protein list and
damping factor, a separate propagation would be run and the result will
be stored in the provided result folder (pvalues.txt, rwr_scores.txt,
start_seeds.txt). After this, the user can reuse the results (so make
sure you don’t delete them!) and test other parameters by providing flag
**-ru**.


Batch job submission(LSF cluster)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each phuEGO run works with one protein list and one damping factor. If
you have multiple protein lists (e.g., from a set of experiment), and/or
would like to test multiple damping factors, you could use a .sh script
to call phuEGO multiple times, and submit jobs in batch manner. Below we
provide a .sh script for a LSF cluster as an example.

To do so, first create a test_datasets.txt file that store the path to
all your protein list files:

.. code:: text

   path/to/protein_list_1.txt
   path/to/protein_list_2.txt
   path/to/protein_list_3.txt

Then submit your jobs using the following .sh script. The output will be
organized into a two-layer folder structure under your specified
result_dir. Modify argument value to suit your need.

.. code:: bash

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


Using 





.. _remove_perturbed_node:

Remove perturbation target nodes from reference network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a drugging or a knockout experiment, one might want to removed the
knocked out targets from the reference network before performing network
propagation, assuming that they are no longer functional. To do so, one
could specify a .csv file as below, and provide to phuEGO:

.. code:: text

   UniprotID_1,UniprotID_2,UniprotID_3
   UniprotID_1,UniprotID_2,UniprotID_3

Here, row 1 is a list of targets to be removed from the network
propagation of upregulated input proteins, and row 2 for downregulated.
Normally, one would expect these to be the same. The list can be
provided as following:

.. code:: bash

   # Performing a run for the first time. Set damping factor to be 0.85, kde_cutoff to be 0.85, and genesets to be 'KEGG'.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" -tpath "path/to/protein_list.txt" -ipath "path/to/targets_list.csv" -d 0.85 -k 0.85 -fg "K" 


Run phuEGO with user defined layers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~