Testing the package
===================

Mock run
~~~~~~~~

Phuego contains a test dataset, which is a list of differentially
phosphorylated proteins between untreated cells and those treated with
epidermal growth factor (EGF) (PMID: 19651622). To view the test
dataset, see `here <#2-checking-the-input-data>`__.

To understand the `input <#2-input>`__ and `output <#3-output>`__ of
phuEGO, the user could either perform a mock run or a test run using the
test dataset. A mock run performs 10 propagations on randomized
networks, and thus isn’t sufficient for reliable statistics, but
typically finishs within a few minutes. A test run is a complete run
with 1000 propagation on randomized networks, and typically takes > 1hr
to finish.

.. code:: bash

   # Performing a mock run.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock

   # Performing a mock run, with two kde_cutoff and two geneset database for geneset enrichment analysis.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_mock -k 0.85 -k 0.9 -fg "K" -fg "B"

   # Performing a test run.
   phuego -sf "path/to/support_data_folder/" -rf "path/to/desired_result_folder/" --run_test


Test run
~~~~~~~~

