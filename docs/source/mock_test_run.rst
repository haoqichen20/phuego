Testing the package
===================

.. _example_input:

The phuEGO package contains a test dataset, which is a list of differentially
phosphorylated proteins (after :ref:`phosphosite aggregation <aggregation>`) 
between untreated cells 
and those treated with epidermal growth factor (EGF) (PMID: 19651622). 
To view the test dataset, run the following Python code:

.. code-block:: python

   from phuego import load_test_example

   # Path and a dataframe of the test dataset.
   test_path, test_df = load_test_example()
   print(test_path) #user can use this to access the file.
   print(test_df)


Mock run
~~~~~~~~

A mock run performs 10 network propagations on the randomized networks using 
the test dataset as input, which isn't sufficient for reliable statistics. But 
as it typically finishes within a few minutes, it can help understanding the package.

**Important note:** carefully provide the kde_cutoff(-k) and fisher genesets(-fg)
as list of float numbers or strings.
Refer to the command help messages for further details.

.. code-block:: bash

   # Folders won't be automatically created if not existing.
   mkdir "path/to/desired_result_folder"

   # Performing a mock run with default parameter values.
   phuego mock\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\

   # Showcasing the usage of parameters.
   phuego mock\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\
    -ds 0.85\
    -de 0.85\
    -dm 0.85\
    -k "[0.85, 0.9]"\
    -fg "['K', 'B']"\

.. warning::

   For Windows users, please provide the paths with forward slash ‘/’ 
   instead of backward slash ‘\\’.

Test run
~~~~~~~~

A test run is similar to the mock run, but performs a complete analysis with 
1000 propagation on randomized networks. It takes longer to finish (typically 10-30 min).

.. code-block:: bash

   # Performing a test run with default parameters.
   phuego test\
    -sf "path/to/support_data_folder/"\
    -rf "path/to/desired_result_folder/"\