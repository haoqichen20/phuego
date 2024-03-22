Running phuEGO
==============

The main analysis of phuEGO is implemented by the command line interface (CLI). 
Paths for the supporting dataset, input data, and the desired result folder 
should be provided, along with other parameters. Check all parameters of the CLI with:

.. code-block:: bash

   phuego main -h

Most parameters have default values assigned. The parameters are divided into 
several groups for clarity, following the flow of the analysis. 
Below, we showcase the usage of the software with a few examples.

.. _CLI:

Analysis with default parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyzing an input data named protein_list.txt.

.. code-block:: bash

   phuego main\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"


.. _reuse:

Reusing seed propagation results, kde_cutoff
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As the seed_propagation step is the most time consuming step, phuEGO save
the intermediate result in three files
(start_seeds.txt, rwr_scores.txt, pvalues.txt).
The user can reuse them to skip the seed propagation step
by providing flag **-ru**. This allows testing downstream steps with
different parameters conveniently.

The parameter **-k** is a threshold that dictates the signal contraction steps (ego_decomposition), 
and has the most influence on the results. Larger kde_cutoff results in smaller networks. phuEGO allows
testing various cutoff values in one execution.

.. container::
   
   Analyzing the same input data while re-using network propagation result.
   Testing three different KDE_cutoff values, perform over-representation analysis
   against different genesets, and export the network in a different format.

   .. code-block:: bash
   
      phuego main\
       -sf "path/to/support_data_folder/"\
       -rf "path/to/desired_result_folder/"\
       -tpath "path/to/protein_list.txt"\
       -ru\
       -k "[0.5, 0.7, 0.85]"\
       -fg "['K', 'R', 'B']"\
       -nf "ncol"


.. _damping:

Damping factors
~~~~~~~~~~~~~~~

Damping factor dictates the behaviour of the random walk with restart process, 
which is implemented by the 
`personalized_pagerank <https://igraph.org/python/doc/api/igraph._igraph.GraphBase.html#personalized_pagerank>`
algorithm. Specifically, at every step of the random walk, there is a probability of **1 - damping** of restarting.
Since the personalized_pagerank algorithm is used in three steps of phuEGO, there are three parameters for this. 
These parameters have limited influence of the final results. 
We suggest keeping all damping factors at default value of 0.85, following the original pagerank implementation, 
or testing a couple of values for damping_seed only, which has the most influence of the three.
Analysis with different damping factor values need to be submitted as separate jobs.

.. container::

   .. code-block:: bash
   
      phuego main\
       -sf "path/to/support_data_folder/"\
       -rf "path/to/desired_result_folder/"\
       -tpath "path/to/protein_list.txt"\
       -ds 0.7\
       -k "[0.5, 0.7, 0.85]"\

.. _user_defined_layers:

User defined layers
~~~~~~~~~~~~~~~~~~~

phuEGO by default separate proteins into three layers: tyrosine kinases, 
serine/threonine kinases, and others, and perform random-walk-with-restart 
(pagerank) within each layers. Users can provide their own definition of layers,
as well as running the method with one or two layers where they see fit. 

The customized layer definition can be provided as a tab-delimited text file, with the first row as layer1 and 
the second row as layer2, and the layer namesin the first column. Uniprot ID of proteins should be used. 
For three layer classification, define two layers and other proteins will be classify as a third layer. 
For two layer classification, define one layer only. If the user wish to test the list as one layer, this is 
available through the -ld parameter.

Format of layer definition:

.. code-block::

   layer1_name uniprotID_1 uniprotID_2 uniprotID_3 
   layer2_name uniprotID_4 uniprotID_5 uniprotID_5

Example code:

.. code-block:: bash

   # Run with customized layers
   phuego main\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"\
    -ld "custom"\
    -ldpath "path/to/layer_definition.txt"

   # Run as one layer.
   phuego main\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"\
    -ld "one"


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

   phuego main\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -tpath "path/to/protein_list.txt"\
    -ipath "path/to/targets_list.csv"



.. _batch_job:

Batch job submission
~~~~~~~~~~~~~~~~~~~~

Each phuEGO run works with one protein list and one set of damping factors. If
you have multiple protein lists (e.g., from a set of experiment), and/or
would like to test multiple damping factors, you could submit a job batch. Below we
provide a .sh script for a LSF cluster as an example.

To do so, first create a test_datasets.txt file that store the path to
all your protein list files:

.. code-block::

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

           # Use 4 cores to accelerate.
           bsub -n 4 -M 4096 -R "rusage[mem=4096]" -o log.txt -e err.txt -J "$job_name" \
           phuego main\
            -sf "Path/to/support_data/"\
            -rf "$damping_dir"\
            -tpath "$line"\
            -d $damping
       done
   done