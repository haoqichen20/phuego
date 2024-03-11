Preparing supporting data
=========================

The supporting datasets contain the base network, randomized networks,
semantic similarity, geneset database etc. For details, refer to the
associated `publication <#6-citation>`__.


Download reference network with semsim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using phuEGO **for the first time**, a support dataset that
contains three zipped files (https://zenodo.org/record/8094690) need to
be downloaded. We recommend using the zenodo-get packages:

.. code:: bash

   # Install zenodo_get to download the supporting dataset from Zenodo.
   pip install zenodo-get
   cd path/to/desired/folder
   zenodo_get 8094690

Decompressed the files into the support datafolder using the following
Python code or other tools, such as tar. **Do not modify the folder
structure.**

.. code:: python

   import tarfile

   support_data_folder = "path/to/desired/folder" #Change this.

   def decompress_tar_gz(file_path, output_dir):
       with tarfile.open(file_path, "r:gz") as tar:
           tar.extractall(path=output_dir)       
           
   paths = ["./gic_sim.tar.gz", "./fisher_etc.tar.gz", "./networks.tar.gz"]
   for path in paths:
       decompress_tar_gz(path, support_data_folder)


Creating your own reference network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~