Preparing supporting dataset
============================

After :doc:`installation <installation>`, a supporting dataset need to be prepared. This 
contains the base and randomized networks populated with pre-calculated protein-protein 
semantic similarity (semsim), as well as various genesets for enrichment analysis.
For details, refer to the associated :ref:`publication <citation>`.

The supporting dataset for human (ID 9606) with a full 
`IntAct <https://www.ebi.ac.uk/intact/home>`__ network can be 
directly :ref:`downloaded <download>`. Alternatively, the users can :ref:`create <create>` 
customized reference networks (e.g., for other species or networks) using 
a `Nextflow pipeline <https://url_to_be_added>`__ provided separately.


.. _download:

Downloading the pre-made dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container::

    The pre-made supporting data is deposited on 
    `Zenodo <https://zenodo.org/record/8094690>`__. We recommend using the 
    zenodo-get packages for downloading. Note that the dataset require 
    ~15G storage.

    .. code-block:: bash

        pip install zenodo-get
        cd path/to/desired/folder #change this
        zenodo_get 8094690

.. container::

    Decompressed the files into the support datafolder using the following
    Python code or other tools, such as tar. **Do not modify the folder
    structure.**

    .. code-block:: python

        import tarfile

        support_data_folder = "path/to/desired/folder" #Change this.

        def decompress_tar_gz(file_path, output_dir):
            with tarfile.open(file_path, "r:gz") as tar:
                tar.extractall(path=output_dir)       
                
        paths = ["./gic_sim.tar.gz", "./fisher_etc.tar.gz", "./networks.tar.gz"]
        for path in paths:
            decompress_tar_gz(path, support_data_folder)


.. _create:

Creating dataset with customized reference network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will be updated once the `Nextflow pipeline <https://url_to_be_added>`__ 
is ready and online.