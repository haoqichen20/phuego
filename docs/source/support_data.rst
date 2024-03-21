Preparing supporting dataset
============================

After :doc:`installation <installation>`, a supporting dataset need to be prepared. This 
contains the base and randomized networks populated with pre-calculated protein-protein 
semantic similarity (semsim) for the network propagation, as well as various genesets for enrichment analysis.
For details, refer to the :ref:`publication <citation>`.

The supporting dataset for human (ID 9606) with a full 
`IntAct <https://www.ebi.ac.uk/intact/home>`__ network can be 
directly :ref:`downloaded <download>`. **Alternatively**, the users can :ref:`create <create>` 
customized reference networks (e.g., for other species or networks) using 
a `Nextflow pipeline <https://github.com/haoqichen20/phuego_support>`__ provided separately.


.. _download:

Downloading the pre-made dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container::

    The pre-made supporting data (slim version containing only "gic" semantic similarity)
    is deposited on `Zenodo <https://zenodo.org/records/10839901>`__. 
    We recommend using the zenodo-get packages for downloading. 
    Note that the dataset require ~15G storage.

    .. code-block:: bash

        pip install zenodo-get
        cd path/to/desired/folder #change this
        zenodo_get 10839901

.. container::

    Decompressed into a desired path.
    Note that the compressed file contain a top-layer folder, 
    to ensure the integrity of folder structure after decompression. **Do not
    modify the folder structure.**

    .. code-block:: bash
        
        cd path/to/desired/folder #Change this.
        tar -xzf support_data_slim.tar.gz

    If tar is not available, the following Python code might be used:

    .. code-block:: python

        import tarfile

        target_folder = "path/to/desired/folder" #Change this.

        def decompress_tar_gz(file_path, output_dir):
            with tarfile.open(file_path, "r:gz") as tar:
                tar.extractall(path=output_dir)       
                
        decompress_tar_gz("./support_data_slim.tar.gz", target_folder)


.. _create:

Creating dataset with customized reference network (Beta)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Nextflow pipeline <https://github.com/haoqichen20/phuego_support>`__ for creating customized
support data is currently being tested. 
User that wishes to use this pipeline can contact the authors for further information.

A user-defined network will be required as input. This can be of human or other model species. 
The pipeline will then calculate the semantic similarity
for all protein pairs, create randomized network and populate the semantic similarity. 
The user can also use different types of semantic similarity metric where they see suitable.

For large networks such as the `IntAct <https://www.ebi.ac.uk/intact/home>`__ used in the publication, 
the pipeline would require large memory and multiple cores to be executed. Thus, the pipeline should be
run on a server rather than personal computers.
