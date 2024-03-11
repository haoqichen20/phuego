Testing the package
===================

The phuEGO package contains a test dataset, which is a list of differentially
phosphorylated proteins (after phosphosite aggregation) between untreated cells 
and those treated with epidermal growth factor (EGF) (PMID: 19651622). 
To view the test dataset, open a Jupyter notebook and run the following:

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
as it typically finishes within minutes, it can help understanding the package.

.. code-block:: bash

   # Performing a mock run.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock

   # Performing a mock run, with two kde_cutoff and two geneset database for geneset enrichment analysis.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock -k 0.85 -k 0.9 -fg "K" -fg "B"


Test run
~~~~~~~~

A test run is similar to the mock run, but performs a complete analysis with 
1000 propagation on randomized networks. It takes longer to finish.

.. code-block:: bash

   # Performing a test run.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_test