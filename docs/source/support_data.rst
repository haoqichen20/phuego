Preparing supporting data
=========================

After :doc:`installation`, a supporting dataset need to be prepared. This 
contains the base and randomized networks populated with pre-calculated protein-protein 
semantic similarity (semsim), as well as various genesets for enrichment analysis.
For details, refer to the associated :ref:`publication`.

The supporting dataset for human (ID 9606) with a full 
`IntAct <https://www.ebi.ac.uk/intact/home>`__ network can be 
directly :ref:`downloaded <download>`. Alternatively, the users can :ref:`create <create>` 
customized reference networks (e.g., for other species or networks) using 
a `Nextflow pipeline <https://url_to_be_added>`__ provided separately.


.. _download:

Downloading the human IntAct reference network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _create:

Creating your own reference network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

