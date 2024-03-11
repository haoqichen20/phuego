Python APIs
===========

The functions in the phuEGO package can be imported and used in your own
python scripts. This makes it easier for integrating phuEGO into your
own workflow.

.. _downloading-supporting-dataset-1:

1). Downloading supporting dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

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

2). Checking the input data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from phuego import load_test_example

   # The function return the absolute path of the test dataset, and the dataset itself as a dataframe.
   test_path, test_df = load_test_example()
   print(test_df)

.. _performing-a-testmock-run-with-the-phuego-test-dataset-1:

3). Performing a test/mock run with the phuEGO test dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Performing a mock run.

.. code:: python

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

Performing a full test run.

.. code:: python

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

4). Running your own protein list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run phuEGO on your own protein list, simply provide the **test_path**
to the above code, and remove the **load_test_example()** line.

To reuse network propagation result and explore different KDE cutoff /
genesets, set **use_existing_rwr = True** (also see above).
