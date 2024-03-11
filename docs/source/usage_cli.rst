Running phuEGO
==============

The main analysis of phuEGO is implemented by the CLI. Paths for the supporting 
dataset, input data, and the desired result folder should be provided, along with
other parameters. These can be either directly provided using command 
line flags, or aggregated in a text file (hence **argument file**) for clarity. 

The flags and arguments for the command line interface can be divided into two sets.
The first set includes parameters of the methods. 
The second set controls the execution of the program.
For an overview of CLI arguments, use:

.. code-block:: bash

   phuego --help


Command line interface
~~~~~~~~~~~~~~~~~~~~~~

Analyzing an input data named protein_list.txt, with damping factor of 0.85, 
kde_cutoff of 0.85, and geneset of "KEGG".

.. code-block:: bash

   phuego\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"\
    -d 0.85\
    -k 0.85\
    -fg "K"

.. warning::

   For Windows users, please provide the paths with forward slash ‘/’ 
   instead of backward slash ‘\\’.


Argument file
~~~~~~~~~~~~~

Example of an argument file:

.. code-block::

   TODO


Using the argument file:

.. code-block:: bash

   TODO


Re-using network propagation results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As the network propagation step is the most time consuming step, phuEGO save
the intermediate result in three files
(start_seeds.txt, rwr_scores.txt, pvalues.txt).
The user can reuse them on the same input data with the same damping factor,
by providing flag **-ru**. This allows testing downstream steps with
different parameters (e.g., different kde_cutoff) conveniently.

.. container::
   
   Analyzing the same input data while re-using network propagation result.
   Testing two different KDE_cutoff values, two different genesets, 
   and export the network in a different format.

   .. code-block:: bash
   
      phuego\
       -sf "path/to/support_data_folder/"\
       -rf "path/to/desired_result_folder/"\
       -tpath "path/to/protein_list.txt"\
       -ru\
       -k 0.8 -k 0.9\
       -fg "C" -fg "B"\
       -nf "edgelist"


User defined layers
~~~~~~~~~~~~~~~~~~~

phuEGO by default separate proteins into three layers: tyrosine kinases, 
serine/threonine kinases, and others, and perform random-walk-with-restart 
(pagerank) within each layers. Users can provide their own definition of layers,
as well as running the method with one or two layers where they see fit. 

The layers can be defined as below:

..code-block::

   TODO

And it can be run with the following:

.. code-block:: bash

   TODO


.. _remove_perturbed_node:

Removing perturbed nodes
~~~~~~~~~~~~~~~~~~~~~~~~

In a drugging or a knockout experiment, one might want to removed the
knocked out targets from the reference network before performing network
propagation, assuming that they are no longer present or functional. 
To do so, one could specify a .csv file as below, and provide to phuEGO:

.. code-block::

   UniprotID_1,UniprotID_2,UniprotID_3
   UniprotID_1,UniprotID_2,UniprotID_3

Here, row 1 is a list of targets to be removed from the network
propagation of upregulated input proteins, and row 2 for downregulated.
Normally, one would expect these to be the same. The list can be
provided as following:

.. code:: bash

   phuego\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"\
    -ipath "path/to/targets_list.csv"\
    -d 0.85\
    -k 0.85\
    -fg "K" 


Batch job submission
~~~~~~~~~~~~~~~~~~~~

Each phuEGO run works with one protein list and one damping factor. If
you have multiple protein lists (e.g., from a set of experiment), and/or
would like to test multiple damping factors, you could use a .sh script
to call phuEGO multiple times, and submit jobs in batch manner. Below we
provide a .sh script for a LSF cluster as an example.

To do so, first create a test_datasets.txt file that store the path to
all your protein list files:

.. code-block:: text

   path/to/protein_list_1.txt
   path/to/protein_list_2.txt
   path/to/protein_list_3.txt

Then submit your jobs using the following .sh script. The output will be
organized into a two-layer folder structure under your specified
result_dir. Modify argument value to suit your need.

.. code-block:: bash

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
           bsub -n 4 -M 4096 -R "rusage[mem=4096]" -o log.txt -e err.txt -J "$job_name" \
           phuego\
            -sf "Path/to/support_data/"\
            -rf "$damping_dir"\
            -tpath "$line"\
            -d $damping\
            -k 0.85 -k 0.9\
            -fg "B" -fg "K"\
            -nf "graphml"
       done
   done